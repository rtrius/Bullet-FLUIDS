/** hashgrid.cpp

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

#include "hashgrid.h"

#include <cstring> //memcpy()

#include "LinearMath/btVector3.h"
#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

#include "FluidParticles.h"

////////////////////////////////////////////////////////////////////////////////
/// class HashGrid
////////////////////////////////////////////////////////////////////////////////
struct HashIndexPair
{
	GridHash m_hash;
	int m_index;
	
	HashIndexPair() {}
	HashIndexPair(GridHash hash, int index) : m_hash(hash), m_index(index) {}
};
struct HashIndexPair_SortPredicate 
{ 
	inline bool operator() (const HashIndexPair &a, const HashIndexPair &b) const 
	{
		return (a.m_hash < b.m_hash);
	}
};

void rearrangeToMatchSortedHashes(const btAlignedObjectArray<HashIndexPair> &sortedHashes, btAlignedObjectArray<btVector3> &out_rearranged)
{
	static btAlignedObjectArray<btVector3> result;
	result.resize( sortedHashes.size() );
	
	for(int i = 0; i < sortedHashes.size(); ++i)
	{
		int oldIndex = sortedHashes[i].m_index;
		int newIndex = i;
			
		result[newIndex] = out_rearranged[oldIndex];
	}
	
	out_rearranged = result;
}
void sortHashGrid(FluidParticles *fluids, btAlignedObjectArray<HashIndexPair> *hashes)
{
	{
		BT_PROFILE("sortHashGrid() - quickSort");
		hashes->quickSort( HashIndexPair_SortPredicate() );
	}
	
	{
		BT_PROFILE("sortHashGrid() - move data");
		
		//Other arrays in fluids are discarded and recalculated when FluidSystem::stepSimulation() is called
		rearrangeToMatchSortedHashes(*hashes, fluids->m_pos);
		rearrangeToMatchSortedHashes(*hashes, fluids->m_vel);
		rearrangeToMatchSortedHashes(*hashes, fluids->m_vel_eval);
		rearrangeToMatchSortedHashes(*hashes, fluids->m_externalAcceleration);
	}
}


void HashGrid::insertParticles(FluidParticles *fluids)
{
	static btAlignedObjectArray<HashIndexPair> hashes;
	{
		BT_PROFILE("hashgrid() - generate");
		hashes.resize( fluids->size() );
		
		resetPointAabb();
		for(int i = 0; i < fluids->size(); ++i) 
		{
			updatePointAabb(fluids->m_pos[i]);
		
			HashGridIndicies indicies = generateIndicies(fluids->m_pos[i]);
			
			hashes[i] = HashIndexPair( indicies.getHash(), i );
		}
	}
	
	
	//Sort fluidSystem and hashes by hashes
	{
		BT_PROFILE("hashgrid() - sort");
		sortHashGrid(fluids, &hashes);
	}
	
	m_activeCells.resize(0);
	m_cellContents.resize(0);
	{
		BT_PROFILE("hashgrid() - find unique");
		
		//Iterate through hashes to find the unique GridHashes,
		//and the index ranges(fluids[] index) at which each hash appears
		if( hashes.size() ) 
		{
			m_activeCells.push_back( hashes[0].m_hash );
			m_cellContents.push_back( HashGridCell() );
			m_cellContents[0].m_firstIndex = 0;
			m_cellContents[0].m_lastIndex = 0;
			
			for(int i = 1; i < hashes.size(); ++i)
			{
				if( hashes[i].m_hash != hashes[i - 1].m_hash )
				{
					m_activeCells.push_back( hashes[i].m_hash );
					m_cellContents.push_back( HashGridCell() );
					
					int lastIndex = m_cellContents.size() - 1;
					m_cellContents[lastIndex].m_firstIndex = i;
					m_cellContents[lastIndex].m_lastIndex = i;
					
					//
					m_cellContents[lastIndex - 1].m_lastIndex = i - 1;
				}
			}
			
			int hashesLastIndex = hashes.size() - 1;
			if( hashes[hashesLastIndex].m_hash == hashes[hashesLastIndex - 1].m_hash )
			{
				int uniqueLastIndex = m_cellContents.size() - 1;
				m_cellContents[uniqueLastIndex].m_lastIndex = hashesLastIndex;
			}
		}
	}
	
	generateCellProcessingGroups();
}

void HashGrid::findCells(const btVector3 &position, btScalar radius, FindCellsResult *out_gridCells) const
{
	btVector3 sphereMin( position.x() - radius, position.y() - radius, position.z() - radius );
	findAdjacentGridCells( generateIndicies(sphereMin), out_gridCells );
}


void HashGrid::getGridCellIndiciesInAabb(const btVector3 &min, const btVector3 &max, btAlignedObjectArray<int> *out_indicies) const
{
	HashGridIndicies minIndicies = generateIndicies(min);
	HashGridIndicies maxIndicies = generateIndicies(max);

	for(HashGridIndex z = minIndicies.z; z <= maxIndicies.z; ++z)
		for(HashGridIndex y = minIndicies.y; y <= maxIndicies.y; ++y)
			for(HashGridIndex x = minIndicies.x; x <= maxIndicies.x; ++x)
			{
				HashGridIndicies current;
				current.x = x;
				current.y = y;
				current.z = z;
			
				//findBinarySearch() returns m_activeCells.size() on failure
				int gridCellIndex = m_activeCells.findBinarySearch( current.getHash() );
				if( gridCellIndex != m_activeCells.size() ) out_indicies->push_back(gridCellIndex);
			}
}

void HashGrid::getResolution(int *out_resolutionX, int *out_resolutionY, int *out_resolutionZ) const
{
	*out_resolutionX = HASH_GRID_INDEX_RANGE;
	*out_resolutionY = HASH_GRID_INDEX_RANGE;
	*out_resolutionZ = HASH_GRID_INDEX_RANGE;
}
void HashGrid::removeFirstParticle(int gridCellIndex, const btAlignedObjectArray<int> &nextFluidIndex)
{
	//Effect of removeFirstParticle() if there is only 1 particle in the cell:
	//Since the iteration loop for a HashGridCell is effectively:
	//	HashGridCell cell;
	//	for(int n = cell.m_firstIndex; n <= cell.m_lastIndex; ++n)
	//
	//The loop will not execute when (cell.m_firstIndex > cell.m_lastIndex),
	//which has the same effect as removing the HashGridCell from m_cellContents.
	++m_cellContents[gridCellIndex].m_firstIndex;
}

HashGridIndicies HashGrid::generateIndicies(const btVector3 &position) const
{
	//Using only 'result.x = static_cast<HashGridIndex>( position.x() / m_gridCellSize )'
	//would cause positions in (-m_gridCellSize, m_gridCellSize) to hash to 0;
	//that is, cell (0,0,0) would be twice as large as desired.
	//
	//To resolve this, define the indicies such that:
	//[0, m_gridCellSize) hashes to 0
	//[-m_gridCellSize, 0) hashes to -1

	HashGridIndicies result;
	
	btVector3 discretePosition = position / m_gridCellSize;
	result.x = static_cast<HashGridIndex>( (position.x() >= 0.0f) ? discretePosition.x() : floor(discretePosition.x()) );
	result.y = static_cast<HashGridIndex>( (position.y() >= 0.0f) ? discretePosition.y() : floor(discretePosition.y()) );
	result.z = static_cast<HashGridIndex>( (position.z() >= 0.0f) ? discretePosition.z() : floor(discretePosition.z()) );
	
	return result;
}

void HashGrid::findAdjacentGridCells(HashGridIndicies indicies, FindCellsResult *out_gridCells) const
{	
	HashGridIndicies cellIndicies[8];
	cellIndicies[0] = indicies;
	
	for(int i = 1; i < 4; ++i) cellIndicies[i] = cellIndicies[0];
	cellIndicies[1].x++;
	cellIndicies[2].y++;
	cellIndicies[3].x++;
	cellIndicies[3].y++;
	
	for(int i = 0; i < 4; ++i) 
	{
		cellIndicies[i+4] = cellIndicies[i];
		cellIndicies[i+4].z++;
	}
	
	for(int i = 0; i < 8; ++i) 
	{
		const HashGridCell *cell = getCell( cellIndicies[i].getHash() );
		
		out_gridCells->m_iterators[i].m_firstIndex = (cell) ? cell->m_firstIndex : INVALID_PARTICLE_INDEX;
		out_gridCells->m_iterators[i].m_lastIndex = (cell) ? cell->m_lastIndex : INVALID_PARTICLE_INDEX;
	}
}

