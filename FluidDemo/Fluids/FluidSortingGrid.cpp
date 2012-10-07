/* FluidSortingGrid.cpp

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

#include "FluidSortingGrid.h"

#include "LinearMath/btVector3.h"
#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

#include "FluidParticles.h"


struct ValueIndexPair_SortPredicate 
{
	inline bool operator() (const ValueIndexPair &a, const ValueIndexPair &b) const 
	{
		return (a.m_value < b.m_value);
	}
};

void rearrangeToMatchSortedValues(const btAlignedObjectArray<ValueIndexPair> &sortedValues, btAlignedObjectArray<btVector3> &out_rearranged)
{
	static btAlignedObjectArray<btVector3> result;
	result.resize( sortedValues.size() );
	
	for(int i = 0; i < sortedValues.size(); ++i)
	{
		int oldIndex = sortedValues[i].m_index;
		int newIndex = i;
			
		result[newIndex] = out_rearranged[oldIndex];
	}
	
	out_rearranged = result;
}
void sortParticlesByValues(FluidParticles *fluids, btAlignedObjectArray<ValueIndexPair> *values)
{
	{
		BT_PROFILE("sortParticlesByValues() - quickSort");
		values->quickSort( ValueIndexPair_SortPredicate() );
	}
	
	{
		BT_PROFILE("sortParticlesByValues() - move data");
		
		//Other arrays in fluids are discarded and recalculated when FluidWorld::stepSimulation() is called
		rearrangeToMatchSortedValues(*values, fluids->m_pos);
		rearrangeToMatchSortedValues(*values, fluids->m_vel);
		rearrangeToMatchSortedValues(*values, fluids->m_vel_eval);
		rearrangeToMatchSortedValues(*values, fluids->m_externalAcceleration);
	}
}


void FluidSortingGrid::insertParticles(FluidParticles *fluids)
{
	static btAlignedObjectArray<ValueIndexPair> pairs;
	{
		BT_PROFILE("FluidSortingGrid() - generate");
		pairs.resize( fluids->size() );
		
		resetPointAabb();
		for(int i = 0; i < fluids->size(); ++i) 
		{
			updatePointAabb(fluids->m_pos[i]);
		
			SortGridIndicies indicies = generateIndicies(fluids->m_pos[i]);
			
			pairs[i] = ValueIndexPair( indicies.getValue(), i );
		}
	}
	
	
	//Sort fluidSystem and values by m_value(s) in pairs
	{
		BT_PROFILE("FluidSortingGrid() - sort");
		sortParticlesByValues(fluids, &pairs);
	}
	
	m_activeCells.resize(0);
	m_cellContents.resize(0);
	{
		BT_PROFILE("FluidSortingGrid() - find unique");
		
		//Iterate through pairs to find the unique SortGridValue(s),
		//and the index ranges(fluids[] index) at which each value appears
		if( pairs.size() ) 
		{
			m_activeCells.push_back( pairs[0].m_value );
			m_cellContents.push_back( FluidGridIterator() );
			m_cellContents[0].m_firstIndex = 0;
			m_cellContents[0].m_lastIndex = 0;
			
			for(int i = 1; i < pairs.size(); ++i)
			{
				if( pairs[i].m_value != pairs[i - 1].m_value )
				{
					m_activeCells.push_back( pairs[i].m_value );
					m_cellContents.push_back( FluidGridIterator() );
					
					int lastIndex = m_cellContents.size() - 1;
					m_cellContents[lastIndex].m_firstIndex = i;
					m_cellContents[lastIndex].m_lastIndex = i;
					
					//
					m_cellContents[lastIndex - 1].m_lastIndex = i - 1;
				}
			}
			
			int valuesLastIndex = pairs.size() - 1;
			if( pairs[valuesLastIndex].m_value == pairs[valuesLastIndex - 1].m_value )
			{
				int uniqueLastIndex = m_cellContents.size() - 1;
				m_cellContents[uniqueLastIndex].m_lastIndex = valuesLastIndex;
			}
		}
	}
	
	generateCellProcessingGroups();
}

void FluidSortingGrid::findCells(const btVector3 &position, btScalar radius, FluidSortingGrid::FoundCells *out_gridCells) const
{
	btVector3 sphereMin( position.x() - radius, position.y() - radius, position.z() - radius );
	findAdjacentGridCells( generateIndicies(sphereMin), out_gridCells );
}


void FluidSortingGrid::getGridCellIndiciesInAabb(const btVector3 &min, const btVector3 &max, btAlignedObjectArray<int> *out_indicies) const
{
	SortGridIndicies minIndicies = generateIndicies(min);
	SortGridIndicies maxIndicies = generateIndicies(max);

	for(SortGridIndex z = minIndicies.z; z <= maxIndicies.z; ++z)
		for(SortGridIndex y = minIndicies.y; y <= maxIndicies.y; ++y)
			for(SortGridIndex x = minIndicies.x; x <= maxIndicies.x; ++x)
			{
				SortGridIndicies current;
				current.x = x;
				current.y = y;
				current.z = z;
			
				//findBinarySearch() returns m_activeCells.size() on failure
				int gridCellIndex = m_activeCells.findBinarySearch( current.getValue() );
				if( gridCellIndex != m_activeCells.size() ) out_indicies->push_back(gridCellIndex);
			}
}

void FluidSortingGrid::internalRemoveFirstParticle(int gridCellIndex)
{
	//Effect of internalRemoveFirstParticle() if there is only 1 particle in the cell:
	//Since the iteration loop for a FluidGridIterator, 
	//when using FluidSortingGrid, is effectively:
	//	FluidGridIterator cell;
	//	for(int n = cell.m_firstIndex; n <= cell.m_lastIndex; ++n)
	//
	//The loop will not execute when (cell.m_firstIndex > cell.m_lastIndex),
	//which has the same effect as removing the FluidGridIterator from m_cellContents.
	++m_cellContents[gridCellIndex].m_firstIndex;
}

SortGridIndicies FluidSortingGrid::generateIndicies(const btVector3 &position) const
{
	//Using only 'result.x = static_cast<SortGridIndex>( position.x() / m_gridCellSize )'
	//would cause positions in (-m_gridCellSize, m_gridCellSize) to convert to 0;
	//that is, cell (0,0,0) would be twice as large as desired.
	//
	//To resolve this, define the indicies such that:
	//[0, m_gridCellSize) converts to 0
	//[-m_gridCellSize, 0) converts to -1

	SortGridIndicies result;
	
	btVector3 discretePosition = position / m_gridCellSize;
	result.x = static_cast<SortGridIndex>( (position.x() >= 0.0f) ? discretePosition.x() : floor(discretePosition.x()) );
	result.y = static_cast<SortGridIndex>( (position.y() >= 0.0f) ? discretePosition.y() : floor(discretePosition.y()) );
	result.z = static_cast<SortGridIndex>( (position.z() >= 0.0f) ? discretePosition.z() : floor(discretePosition.z()) );
	
	btAssert(-HALVED_SORT_GRID_INDEX_RANGE <= result.x && result.x <= HALVED_SORT_GRID_INDEX_RANGE - 1);
	btAssert(-HALVED_SORT_GRID_INDEX_RANGE <= result.y && result.y <= HALVED_SORT_GRID_INDEX_RANGE - 1);
	btAssert(-HALVED_SORT_GRID_INDEX_RANGE <= result.z && result.z <= HALVED_SORT_GRID_INDEX_RANGE - 1);
	
	return result;
}

void FluidSortingGrid::findAdjacentGridCells(SortGridIndicies indicies, FluidSortingGrid::FoundCells *out_gridCells) const
{	
	//INVALID_LAST_INDEX must be lower than INVALID_FIRST_INDEX,
	//such that the below loop will not execute.
	//		for(int i = FI.m_firstIndex; i <= FI.m_lastIndex; ++i)
	const int INVALID_FIRST_INDEX = -1;
	const int INVALID_LAST_INDEX = INVALID_FIRST_INDEX - 1;

	SortGridIndicies cellIndicies[8];
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
		const FluidGridIterator *cell = getCell( cellIndicies[i].getValue() );
		
		out_gridCells->m_iterators[i].m_firstIndex = (cell) ? cell->m_firstIndex : INVALID_FIRST_INDEX;
		out_gridCells->m_iterators[i].m_lastIndex = (cell) ? cell->m_lastIndex : INVALID_LAST_INDEX;
	}
}

void FluidSortingGrid::generateCellProcessingGroups()
{
	//Although a single particle only accesses 2^3 grid cells,
	//the particles within a single grid cell may access up to 
	//3^3 grid cells overall.
	//
	//In order to simultaneously calculate pressure on a per grid 
	//cell basis, with particles being removed from a grid cell
	//as it is processed, the cells must be split into 27
	//sequentially processed groups.

	BT_PROFILE("generateCellProcessingGroups()");
	
	for(int i = 0; i < FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS; ++i) m_cellProcessingGroups[i].resize(0);
		
	for(int cell = 0; cell < getNumGridCells(); ++cell)
	{
		FluidGridIterator FI = getGridCell(cell);
		if( !(FI.m_firstIndex <= FI.m_lastIndex) ) continue;
	
		int index_x, index_y, index_z;
		internalGetIndiciesReduce(cell, &index_x, &index_y, &index_z);

		//For FluidSortingGrid: Convert range from [0, 255] to [1, 256]
		//(SortGridIndicies::getValue() already converts from [-128, 127], to [0, 255])
		index_x += 1;
		index_y += 1;
		index_z += 1;
		
		//For each dimension, place indicies into one of 3 categories such that
		//indicies (1, 2, 3, 4, 5, 6, ...) correspond to categories (1, 2, 3, 1, 2, 3, ...)
		int group = 0;
		
		if(index_x % 3 == 0) group += 0;
		else if(index_x % 2 == 0) group += 1;
		else group += 2;
		
		if(index_y % 3 == 0) group += 0;
		else if(index_y % 2 == 0) group += 3;
		else group += 6;
		
		if(index_z % 3 == 0) group += 0;
		else if(index_z % 2 == 0) group += 9;
		else group += 18;
		
		m_cellProcessingGroups[group].push_back(cell);
	}
}
