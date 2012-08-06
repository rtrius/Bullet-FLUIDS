/** FluidSortingGrid.h

	ZLib license
	This software is provided 'as-is', without any express or implied
	warranty. In no event will the authors be held liable for any damages
	arising from the use of this software.
	
	Permission is granted to anyone to use this software for any purpose,
	including commercial applications, and to alter it and redistribute it
	freely, subject to the following restrictions:
	
	1. The origin of this software must not be misrepresented; you must not
	   claim that you wrote the original software. If you use this software
	   in a product, an acknowledgment in the product documentation would be
	   appreciated but is not required.
	2. Altered source versions must be plainly marked as such, and must not be
	   misrepresented as being the original software.
	3. This notice may not be removed or altered from any source distribution.
*/
#ifndef FLUID_SORTING_GRID_H
#define FLUID_SORTING_GRID_H

#include "LinearMath/btAlignedObjectArray.h"

#include "FluidGrid.h"

class btVector3;
struct FluidParticles;

typedef unsigned int SortGridValue;	//Range must contain SORT_GRID_INDEX_RANGE^3
typedef char SortGridIndex;
const SortGridValue SORT_GRID_INDEX_RANGE = 256;		//2^( 8*sizeof(SortGridIndex) )


struct SortGridCell
{
	int m_firstIndex;
	int m_lastIndex;
};

struct SortGridIndicies		//Contains a world position in units of 'FluidSortingGrid.m_gridCellSize'
{
	SortGridIndex x;		
	SortGridIndex y;
	SortGridIndex z;
private:
	SortGridIndex padding;
	
public:
	bool operator==(const SortGridIndicies& GI) const { return (x == GI.x && y == GI.y && z == GI.z); }
	bool operator!=(const SortGridIndicies& GI) const { return (x != GI.x || y != GI.y || z != GI.z); }
	
	bool operator>(const SortGridIndicies& GI) const
	{
		if(z != GI.z) return (z > GI.z);
		if(y != GI.y) return (y > GI.z);
		return (x > GI.x);
	}	
	bool operator<(const SortGridIndicies& GI) const
	{
		if(z != GI.z) return (z < GI.z);
		if(y != GI.y) return (y < GI.z);
		return (x < GI.x);
	}
	
	SortGridValue getValue() const 
	{
		//Convert range from [-128, 127] to [0, 255] before combining
		return (  (static_cast<int>(x)+128) 
				+ (static_cast<int>(y)+128)*SORT_GRID_INDEX_RANGE 
				+ (static_cast<int>(z)+128)*SORT_GRID_INDEX_RANGE*SORT_GRID_INDEX_RANGE );
	}
};


///SORT_GRID_INDEX_RANGE^3 sized grid that only stores nonempty cells.
class FluidSortingGrid : public FluidGrid
{
	btScalar m_gridCellSize;

	btAlignedObjectArray<SortGridValue> m_activeCells;	//Stores the value of each nonempty grid cell
	btAlignedObjectArray<SortGridCell> m_cellContents;	//Stores the range of indicies that correspond to the values in m_activeCells
	
	btVector3 m_pointMin;
	btVector3 m_pointMax;
	
public:
	FluidSortingGrid() {}
	FluidSortingGrid(btScalar simScale, btScalar simCellSize) { setup(simScale, simCellSize); }

	void setup(btScalar simScale, btScalar simCellSize) 
	{		
		btScalar worldCellSize = simCellSize / simScale;
		m_gridCellSize = worldCellSize;
	}

	virtual void clear() 
	{ 
		m_activeCells.resize(0);
		m_cellContents.resize(0);
	}
	virtual void insertParticles(FluidParticles *fluids);
	virtual void findCells(const btVector3 &position, btScalar radius, FindCellsResult *out_gridCells) const;
	
	virtual FluidGridIterator getGridCell(int gridCellIndex) const
	{
		return FluidGridIterator(m_cellContents[gridCellIndex].m_firstIndex, m_cellContents[gridCellIndex].m_lastIndex);
	}
	virtual void getGridCellIndiciesInAabb(const btVector3 &min, const btVector3 &max, btAlignedObjectArray<int> *out_indicies) const;
	
	virtual FluidGridType getGridType() const { return FT_IndexRange; }
	virtual btScalar getCellSize() const { return m_gridCellSize; }
	
	virtual int getCombinedPosition(int gridCellIndex) const { return m_activeCells[gridCellIndex]; }
	virtual int getNumGridCells() const { return m_activeCells.size(); }
	virtual void getResolution(int *out_resolutionX, int *out_resolutionY, int *out_resolutionZ) const;
	virtual void removeFirstParticle(int gridCellIndex, const btAlignedObjectArray<int> &nextFluidIndex);
	
private:
	const SortGridCell* getCell(SortGridValue value) const
	{
		int index = m_activeCells.findBinarySearch(value);	//findBinarySearch() returns m_activeCells.size() on failure
		return ( index != m_activeCells.size() ) ? &m_cellContents[index] : 0;
	}

	SortGridIndicies generateIndicies(const btVector3 &position) const;

	void findAdjacentGridCells(SortGridIndicies indicies, FindCellsResult *out_gridCells) const;
};

#endif



