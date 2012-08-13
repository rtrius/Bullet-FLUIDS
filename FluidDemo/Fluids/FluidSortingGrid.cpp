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


struct ValueIndexPair
{
	SortGridValue m_value;
	int m_index;
	
	ValueIndexPair() {}
	ValueIndexPair(SortGridValue value, int index) : m_value(value), m_index(index) {}
};
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

void FluidSortingGrid::findCells(const btVector3 &position, btScalar radius, FluidGrid::FoundCells *out_gridCells) const
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

void FluidSortingGrid::internalRemoveFirstParticle(int gridCellIndex, const btAlignedObjectArray<int> &nextFluidIndex)
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

void FluidSortingGrid::findAdjacentGridCells(SortGridIndicies indicies, FluidGrid::FoundCells *out_gridCells) const
{	
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
		
		out_gridCells->m_iterators[i].m_firstIndex = (cell) ? cell->m_firstIndex : INVALID_PARTICLE_INDEX;
		out_gridCells->m_iterators[i].m_lastIndex = (cell) ? cell->m_lastIndex : INVALID_PARTICLE_INDEX_MINUS_ONE;
	}
}

