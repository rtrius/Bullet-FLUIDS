/** hashgrid.h

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
#ifndef HASHGRID_H_INCLUDED
#define HASHGRID_H_INCLUDED

#include "LinearMath/btAlignedObjectArray.h"

#include "FluidGrid.h"

class btVector3;
struct FluidParticles;

typedef unsigned int GridHash;	//Range must contain HASH_GRID_INDEX_RANGE^3
typedef char HashGridIndex;
const GridHash HASH_GRID_INDEX_RANGE = 256;		//2^( 8*sizeof(HashGridIndex) )


struct HashGridCell
{
	int m_firstIndex;
	int m_lastIndex;
};
struct HashGridQueryResult { HashGridCell *m_cells[8]; };

struct HashGridIndicies		//Contains a world position in units of 'HashGrid.m_gridCellSize'
{
	HashGridIndex x;		
	HashGridIndex y;
	HashGridIndex z;
private:
	HashGridIndex padding;
	
public:
	bool operator==(const HashGridIndicies& GI) const { return (x == GI.x && y == GI.y && z == GI.z); }
	bool operator!=(const HashGridIndicies& GI) const { return (x != GI.x || y != GI.y || z != GI.z); }
	
	bool operator>(const HashGridIndicies& GI) const
	{
		if(z != GI.z) return (z > GI.z);
		if(y != GI.y) return (y > GI.z);
		return (x > GI.x);
	}	
	bool operator<(const HashGridIndicies& GI) const
	{
		if(z != GI.z) return (z < GI.z);
		if(y != GI.y) return (y < GI.z);
		return (x < GI.x);
	}
	
	GridHash getHash() const 
	{
		//Convert range from [-128, 127] to [0, 255] before combining
		return (  (static_cast<int>(x)+128) 
				+ (static_cast<int>(y)+128)*HASH_GRID_INDEX_RANGE 
				+ (static_cast<int>(z)+128)*HASH_GRID_INDEX_RANGE*HASH_GRID_INDEX_RANGE );
	}
};


///HASH_GRID_INDEX_RANGE^3 sized grid that only stores nonempty cells.
class HashGrid : public FluidGrid
{
	btScalar m_gridCellSize;

	btAlignedObjectArray<GridHash> m_activeCells;		//Stores the hash of each nonempty grid cell
	btAlignedObjectArray<HashGridCell> m_cellContents;	//Stores the range of indicies that correspond to the hashes in m_activeCells
	
	btVector3 m_pointMin;
	btVector3 m_pointMax;
	
public:
	HashGrid() {}
	HashGrid(btScalar simScale, btScalar simCellSize) { setup(simScale, simCellSize); }

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
	const HashGridCell* getCell(GridHash hash) const
	{
		int index = m_activeCells.findBinarySearch(hash);	//findBinarySearch() returns m_activeCells.size() on failure
		return ( index != m_activeCells.size() ) ? &m_cellContents[index] : 0;
	}

	HashGridIndicies generateIndicies(const btVector3 &position) const;

	void findAdjacentGridCells(HashGridIndicies indicies, FindCellsResult *out_gridCells) const;
};

#endif



