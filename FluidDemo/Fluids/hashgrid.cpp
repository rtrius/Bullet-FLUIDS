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
/// class HashGridSort
////////////////////////////////////////////////////////////////////////////////
void HashGridSort::quickSortInternal(int lo, int hi)	
{
	//	from btAlignedObjectArray

	//  lo is the lower index, hi is the upper index
	//  of the region of array a that is to be sorted
	int i = lo;
	int j = hi;
	int pivot = (lo+hi)/2;

	//  partition
	do
	{
		while ( compare(i, pivot) ) i++;
		while ( compare(pivot, j) ) j--;
		if(i <= j)
		{
			swap(i,j);
			i++; 
			j--;
		}
		
	} while (i <= j);

	//  recursion
	if(lo < j) quickSortInternal(lo, j);
	if(i < hi) quickSortInternal(i, hi);
}



	
void HashGridSort::countingSort(unsigned int byteSignificance)
{
	//Base 256 counting sort; 1 array for each digit
	int count[256];	
	for(int i = 0; i < 256; ++i) count[i] = 0;
	
	for(int i = 0; i < size(); ++i) ++count[ getByte( getValue(i), byteSignificance ) ];
	
	{
		int total = 0;
		for(int i = 0; i < 256; ++i)
		{
			int currentCount = count[i];
			count[i] = total;
			total += currentCount;
		}
	}

	static Fluids fluidResult;
	static btAlignedObjectArray<GridHash> hashResult;
	fluidResult.resize( size() );
	hashResult.resize( size() );
	
	for(int i = 0; i < size(); ++i)
	{
		int key = getByte( getValue(i), byteSignificance );
		
		int index = count[key];
		
		//Commented out fields are discarded and recalculated when FluidSystem::stepSimulation() is called
		fluidResult.m_pos[index] = m_fluids->m_pos[i];
		fluidResult.m_vel[index] = m_fluids->m_vel[i];
		fluidResult.m_vel_eval[index] = m_fluids->m_vel_eval[i];
		//fluidResult.m_sph_force[index] = m_fluids->m_sph_force[i];
		fluidResult.m_externalAcceleration[index] = m_fluids->m_externalAcceleration[i];
		//fluidResult.m_prev_pos[index] = m_fluids->m_prev_pos[i];
		//fluidResult.m_pressure[index] = m_fluids->m_pressure[i];
		//fluidResult.m_density[index] = m_fluids->m_density[i];
		//fluidResult.m_nextFluidIndex[index] = m_fluids->m_nextFluidIndex[i];
		
		hashResult[index] = (*m_hashes)[i];
		
		++count[key];
	}
	
	//	portability issues with memcpy()?
	
	//Swap to avoid copying into input arrays
	const int FLUID_ARRAY_SIZE = sizeof(Fluids);
	char fluidArrayBuffer[FLUID_ARRAY_SIZE];
	memcpy(fluidArrayBuffer, m_fluids, FLUID_ARRAY_SIZE);
	memcpy(m_fluids, &fluidResult, FLUID_ARRAY_SIZE);
	memcpy(&fluidResult, fluidArrayBuffer, FLUID_ARRAY_SIZE);
	
	const int HASH_ARRAY_SIZE = sizeof(btAlignedObjectArray<GridHash>);
	char hashArrayBuffer[HASH_ARRAY_SIZE];
	memcpy(hashArrayBuffer, m_hashes, HASH_ARRAY_SIZE);
	memcpy(m_hashes, &hashResult, HASH_ARRAY_SIZE);
	memcpy(&hashResult, hashArrayBuffer, HASH_ARRAY_SIZE);
}


////////////////////////////////////////////////////////////////////////////////
/// class HashGrid
////////////////////////////////////////////////////////////////////////////////
void HashGrid::insertParticles(Fluids *fluids)
{
	static btAlignedObjectArray<GridHash> hashes;
	{
		BT_PROFILE("hashgrid() - generate");
		hashes.resize( fluids->size() );
		for(int i = 0; i < fluids->size(); ++i) hashes[i] = generateIndicies(fluids->m_pos[i]).getHash();
	}
	
	//Sort fluidSystem and hashes by hashes
	{
		BT_PROFILE("hashgrid() - sort");
		HashGridSort HGS(fluids, &hashes);
		//HGS.quickSort();		//	not functional; does not completely sort arrays
		HGS.radixSort();
	}
	
	m_activeCells.resize(0);
	m_cellContents.resize(0);
	{
		BT_PROFILE("hashgrid() - find unique");
		
		//Iterate through hashes to find the unique GridHashes,
		//and the index ranges(fluids[] index) at which each hash appears
		if( hashes.size() ) 
		{
			m_activeCells.push_back( hashes[0] );
			m_cellContents.push_back( HashGridCell() );
			m_cellContents[0].m_firstIndex = 0;
			m_cellContents[0].m_lastIndex = 0;
			
			for(int i = 1; i < hashes.size(); ++i)
			{
				if( hashes[i] != hashes[i - 1] )
				{
					m_activeCells.push_back( hashes[i] );
					m_cellContents.push_back( HashGridCell() );
					
					int lastIndex = m_cellContents.size() - 1;
					m_cellContents[lastIndex].m_firstIndex = i;
					m_cellContents[lastIndex].m_lastIndex = i;
					
					//
					m_cellContents[lastIndex - 1].m_lastIndex = i - 1;
				}
			}
			
			int hashesLastIndex = hashes.size() - 1;
			if( hashes[hashesLastIndex] == hashes[hashesLastIndex - 1] )
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

	HashGridIndicies cellIndicies[8];
	cellIndicies[0] = generateIndicies(sphereMin);

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
	//Using only 'result.x = static_cast<char>( position.x() / m_gridCellSize )'
	//would cause positions in (-m_gridCellSize, m_gridCellSize) to hash to 0;
	//that is, cell (0,0,0) would be twice as large as desired.
	//
	//To resolve this, define the indicies such that:
	//[0, m_gridCellSize) hashes to 0
	//[-m_gridCellSize, 0) hashes to -1

	HashGridIndicies result;
	
	if( position.x() >= 0.0 ) result.x = static_cast<char>( position.x() / m_gridCellSize );
	else result.x = static_cast<char>( ceil(position.x() / m_gridCellSize) );
	if( position.y() >= 0.0 ) result.y = static_cast<char>( position.y() / m_gridCellSize );
	else result.y = static_cast<char>( ceil(position.y() / m_gridCellSize) );
	if( position.z() >= 0.0 ) result.z = static_cast<char>( position.z() / m_gridCellSize );
	else result.z = static_cast<char>( ceil(position.z() / m_gridCellSize) );

	return result;
}
