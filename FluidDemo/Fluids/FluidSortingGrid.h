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


#ifdef SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED		//Defined in "FluidGrid.h"
	typedef unsigned long long int SortGridValue;	//Range must contain SORT_GRID_INDEX_RANGE^3
	//typedef short int SortGridIndex;
	//const SortGridValue SORT_GRID_INDEX_RANGE = 65536;		//2^( 8*sizeof(SortGridIndex) )
	typedef int SortGridIndex;
	const SortGridValue SORT_GRID_INDEX_RANGE = 2097152;	//2^21
#else
	typedef unsigned int SortGridValue;				//Range must contain SORT_GRID_INDEX_RANGE^3
	typedef char SortGridIndex;
	const SortGridValue SORT_GRID_INDEX_RANGE = 256;			//2^( 8*sizeof(SortGridIndex) )
#endif

//SortGridIndex_large is a signed type with range including all values in:
//	[-HALVED_SORT_GRID_INDEX_RANGE, (HALVED_SORT_GRID_INDEX_RANGE - 1) + HALVED_SORT_GRID_INDEX_RANGE]
//
//	e.g if SORT_GRID_INDEX_RANGE == 256, HALVED_SORT_GRID_INDEX_RANGE == 128
//	then SortGridIndex_large must contain [-128, 127 + 128] == [-128, 255].
//	It is used to convert from SortGridIndex(signed, small range), to SortGridValue(unsigned, large range).
typedef int SortGridIndex_large;	

const SortGridIndex_large HALVED_SORT_GRID_INDEX_RANGE = SORT_GRID_INDEX_RANGE/2;
//const int BITS_PER_SORT_GRID_INDEX = 8 * sizeof(SortGridIndex);


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
		if(y != GI.y) return (y > GI.y);
		return (x > GI.x);
	}	
	bool operator<(const SortGridIndicies& GI) const
	{
		if(z != GI.z) return (z < GI.z);
		if(y != GI.y) return (y < GI.y);
		return (x < GI.x);
	}
	
	SortGridValue getValue() const
	{
		//Convert range 
		//from [-HALVED_SORT_GRID_INDEX_RANGE, HALVED_SORT_GRID_INDEX_RANGE - 1] 
		//  to [0, SORT_GRID_INDEX_RANGE - 1] before combining
		//e.g. from [-128, 127] to [0, 255] before combining	
		SortGridIndex_large signedX = static_cast<SortGridIndex_large>(x) + HALVED_SORT_GRID_INDEX_RANGE;
		SortGridIndex_large signedY = static_cast<SortGridIndex_large>(y) + HALVED_SORT_GRID_INDEX_RANGE;
		SortGridIndex_large signedZ = static_cast<SortGridIndex_large>(z) + HALVED_SORT_GRID_INDEX_RANGE;
		
		SortGridValue unsignedX = static_cast<SortGridValue>(signedX);
		SortGridValue unsignedY = static_cast<SortGridValue>(signedY) * SORT_GRID_INDEX_RANGE;
		SortGridValue unsignedZ = static_cast<SortGridValue>(signedZ) * SORT_GRID_INDEX_RANGE * SORT_GRID_INDEX_RANGE;
		
		return unsignedX + unsignedY + unsignedZ;
	}
};

///SORT_GRID_INDEX_RANGE^3 sized grid that only stores nonempty cells.
class FluidSortingGrid : public FluidGrid
{
	btScalar m_gridCellSize;

	btAlignedObjectArray<SortGridValue> m_activeCells;	//Stores the value of each nonempty grid cell
	btAlignedObjectArray<SortGridCell> m_cellContents;	//Stores the range of indicies that correspond to the values in m_activeCells
	
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
	
	virtual void getIndiciesReduce(int gridCellIndex, int *out_x, int *out_y, int *out_z) const
	{
		splitIndex(SORT_GRID_INDEX_RANGE, SORT_GRID_INDEX_RANGE, m_activeCells[gridCellIndex], out_x, out_y, out_z);
	}
	virtual int getNumGridCells() const { return m_activeCells.size(); }
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



