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

#include "LinearMath/btVector3.h"

#include <cstring> //memcpy()


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
void sortHashGrid(Fluids *fluids, btAlignedObjectArray<HashIndexPair> *hashes)
{
	hashes->quickSort( HashIndexPair_SortPredicate() );
	
	static Fluids fluidResult;
	fluidResult.resize( hashes->size() );
	
	for(int i = 0; i < hashes->size(); ++i)
	{
		int oldIndex = (*hashes)[i].m_index;
		int newIndex = i;
		
		//Commented out fields are discarded and recalculated when FluidSystem::stepSimulation() is called
		fluidResult.m_pos[newIndex] = fluids->m_pos[oldIndex];
		fluidResult.m_vel[newIndex] = fluids->m_vel[oldIndex];
		fluidResult.m_vel_eval[newIndex] = fluids->m_vel_eval[oldIndex];
		//fluidResult.m_sph_force[newIndex] = fluids->m_sph_force[oldIndex];
		fluidResult.m_externalAcceleration[newIndex] = fluids->m_externalAcceleration[oldIndex];
		//fluidResult.m_prev_pos[newIndex] = fluids->m_prev_pos[oldIndex];
		//fluidResult.m_pressure[newIndex] = fluids->m_pressure[oldIndex];
		//fluidResult.m_density[newIndex] = fluids->m_density[oldIndex];
		//fluidResult.m_nextFluidIndex[newIndex] = fluids->m_nextFluidIndex[oldIndex];
	}
	
	//	portability issues with memcpy()?
	//Swap to avoid copying into input arrays
	const int FLUID_ARRAY_SIZE = sizeof(Fluids);
	char fluidArrayBuffer[FLUID_ARRAY_SIZE];
	memcpy(fluidArrayBuffer, fluids, FLUID_ARRAY_SIZE);
	memcpy(fluids, &fluidResult, FLUID_ARRAY_SIZE);
	memcpy(&fluidResult, fluidArrayBuffer, FLUID_ARRAY_SIZE);
}

void HashGrid::insertParticles(Fluids *fluids)
{
	static btAlignedObjectArray<HashIndexPair> hashes;
	{
		BT_PROFILE("hashgrid() - generate");
		hashes.resize( fluids->size() );
		for(int i = 0; i < fluids->size(); ++i) 
		{
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
}

void HashGrid::findCells(const btVector3 &position, btScalar radius, HashGridQueryResult *out_gridCells)
{
	btVector3 sphereMin( position.x() - radius, position.y() - radius, position.z() - radius );
	findAdjacentGridCells( generateIndicies(sphereMin), out_gridCells );
}

void HashGrid::findAdjacentGridCells(HashGridIndicies indicies, HashGridQueryResult *out_gridCells)
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
	
	for(int i = 0; i < 8; ++i) out_gridCells->m_cells[i] = getCell( cellIndicies[i].getHash() );
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
