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


#include "LinearMath/btQuickProf.h"
#include "LinearMath/btAlignedObjectArray.h"

#include "fluid.h"

class btVector3;

typedef unsigned int GridHash;	//Range must contain HASH_GRID_INDEX_RANGE^3
typedef char HashGridIndex;
const GridHash HASH_GRID_INDEX_RANGE = 255;




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
		return (  static_cast<GridHash>(x) 
				+ static_cast<GridHash>(y)*HASH_GRID_INDEX_RANGE 
				+ static_cast<GridHash>(z)*HASH_GRID_INDEX_RANGE*HASH_GRID_INDEX_RANGE );
	}
};


///Place must be within [0,9]
///returns the digit at 10^place (e.g. if place == 2, digit at 'hundreds' is returned)
inline unsigned int getDigit(unsigned int value, unsigned int place)
{
	static const unsigned int POWERS_OF_TEN[] = 
	{
		1, 10, 100,
		1000, 10000, 100000,
		1000000, 10000000, 100000000,
		1000000000		//Billion; assuming unsigned int == 32 bits, max value is ~4 billion
	};
		
	//Truncate digits with lesser significance than place
	value /= POWERS_OF_TEN[place];
	
	//Truncate digits with greater significance than place
	return value - (value/10)*10;
}

class HashGridSort
{
	Fluids *m_fluids;
	btAlignedObjectArray<GridHash> *m_hashes;
	
public:
	HashGridSort(Fluids *fluids, btAlignedObjectArray<GridHash> *hashes) : m_fluids(fluids), m_hashes(hashes) {}
	
	void quickSort() { if( size() > 1 ) quickSortInternal(0, size() - 1); }
	
	void radixSort()
	{
		const int LAST_PLACE = 8;	//Orders of magnitude of 2^32(~4 billion) == 10 digits == places[0,9] == 10^8
		for(int i = 0; i <= LAST_PLACE; ++i) countingSort(i);
	}
	
	int size() const { return m_fluids->size(); }
	
	//Comparison sort
	bool compare(int index0, int index1) const { return (*m_hashes)[index0] < (*m_hashes)[index1]; }
	void swap(int index0, int index1)
	{
		m_fluids->m_pos.swap(index0, index1);
		m_fluids->m_vel.swap(index0, index1);
		m_fluids->m_vel_eval.swap(index0, index1);
		m_fluids->m_sph_force.swap(index0, index1);
		m_fluids->m_externalAcceleration.swap(index0, index1);
		m_fluids->m_prev_pos.swap(index0, index1);
		m_fluids->m_pressure.swap(index0, index1);
		m_fluids->m_density.swap(index0, index1);
		m_fluids->m_nextFluidIndex.swap(index0, index1);
		//m_fluids->m_neighborTable.swap(index0, index1);
		
		m_hashes->swap(index0, index1);
	}
	
	//Radix	sort
	GridHash getValue(int index) const { return (*m_hashes)[index]; }
	
private:
	void quickSortInternal(int lo, int hi);
	
	void countingSort(unsigned int place);	//Sorts by the digit at 10^place
};



///HASH_GRID_INDEX_RANGE^3 sized grid that only stores nonempty cells.
class HashGrid
{
	float m_gridCellSize;

	btAlignedObjectArray<GridHash> m_activeCells;		//Stores the hash of each nonempty grid cell
	btAlignedObjectArray<HashGridCell> m_cellContents;	//Stores the range of indicies that correspond to the hashes in m_activeCells
	
public:
	void setup(float simScale, float simCellSize) 
	{		
		float worldCellSize = simCellSize / simScale;
		m_gridCellSize = worldCellSize;
	}

	void clear() 
	{ 
		m_activeCells.resize(0);
		m_cellContents.resize(0);
	}
	
	void insertParticles(Fluids *fluids);
	
	void findCells(const btVector3 &position, float radius, HashGridQueryResult *out_gridCells);
	HashGridCell* getCell(const GridHash &hash)
	{
		int index = m_activeCells.findBinarySearch(hash);	//findBinarySearch() returns m_activeCells.size() on failure
		return ( index != m_activeCells.size() ) ? &m_cellContents[index] : 0;
	}
	
	HashGridIndicies generateIndicies(const btVector3 &position) const;
};

/*
class FluidIterator
{
public:
	virtual ~FluidIterator() {}

	virtual int getFirstIndex() = 0;
	virtual int nextIndex() = 0;
	
	virtual bool isIndexValid(int index) const = 0;	
};

class HashGridIterator : public FluidIterator
{
	HashGridCell m_cell;
	
	int m_currentIndex;
	
public:
	HashGridIterator(const HashGridCell &C) : m_cell(C) {}

	virtual int getFirstIndex() 
	{
		m_currentIndex = m_cell.m_firstIndex;
		return m_cell.m_firstIndex; 
	}
	virtual int nextIndex() { return ++m_currentIndex; }
	virtual bool isIndexValid(int index) const { return (index <= m_cell.m_lastIndex); }
};
class StaticGridIterator : public FluidIterator
{
	btAlignedObjectArray<Fluid> *m_fluids;
	int m_firstIndex;
	
	int m_currentIndex;
	
public:
	StaticGridIterator(btAlignedObjectArray<Fluid> *fluids, int firstIndex) : m_fluids(fluids), m_firstIndex(firstIndex) {}

	virtual int getFirstIndex()
	{
		m_currentIndex = m_firstIndex;
		return m_firstIndex;
	}
	virtual int nextIndex()
	{ 
		m_currentIndex = (*m_fluids)[m_currentIndex].nextFluidIndex;
		return m_currentIndex; 
	}
	virtual bool isIndexValid(int index) const { return (index != INVALID_PARTICLE_INDEX); }
};

///FluidIterator *FI;
///for( int i = FI->getFirstIndex(); FI->isIndexValid(i); i = FI->nextIndex() )
///{
///	Fluid *f = fluidSystem->getFluid(i);
///}
*/

#endif



