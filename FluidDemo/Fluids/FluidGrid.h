/** FluidGrid.h
	Copyright (C) 2012 Jackson Lee

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
#ifndef FLUID_GRID_ITERATOR_H
#define FLUID_GRID_ITERATOR_H

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"

#include "LinearMath/btQuickProf.h"


#include "FluidParticles.h"

enum FluidGridType
{
	FT_LinkedList,
	FT_IndexRange
};
struct FluidGridIterator
{
	int m_firstIndex;
	int m_lastIndex;	//For FT_LinkedList, contains a value greater than the highest fluid particle index
	
	FluidGridIterator() {}
	FluidGridIterator(int firstIndex, int lastIndex) : m_firstIndex(firstIndex), m_lastIndex(lastIndex) {}
	
	static inline bool isIndexValid(int index, int lastIndex) 
	{
		//Valid index condition for Linked List Grids: (index != INVALID_PARTICLE_INDEX)
		//Valid index condition for Index Range Grids: (index <= lastIndex)
	
		return (index != INVALID_PARTICLE_INDEX && index <= lastIndex); 
	}
	static inline int getNextIndex(int index, const bool isLinkedList, const btAlignedObjectArray<int> &nextFluidIndex)
	{
		return (isLinkedList) ? nextFluidIndex[index] : index + 1;
	}
};

const int RESULTS_PER_GRID_SEARCH = 8;		//Number of grid cells returned from FluidGrid::findCells()
struct FindCellsResult { FluidGridIterator m_iterators[RESULTS_PER_GRID_SEARCH]; };

const int NUM_CELL_PROCESSING_GROUPS = 27;
class FluidGrid
{
protected:
	///AABB calculated from the center of fluid particles, without considering particle radius
	btVector3 m_pointMin;
	btVector3 m_pointMax;
	
	btAlignedObjectArray<int> m_cellProcessingGroups[NUM_CELL_PROCESSING_GROUPS];
	
public:
	virtual ~FluidGrid() {}

	virtual void clear() = 0;
	virtual void insertParticles(FluidParticles *fluids) = 0;
	
	virtual void findCells(const btVector3 &position, btScalar radius, FindCellsResult *out_gridCells) const = 0;
	
	virtual FluidGridIterator getGridCell(int gridCellIndex) const = 0;
	virtual void getGridCellIndiciesInAabb(const btVector3 &min, const btVector3 &max, btAlignedObjectArray<int> *out_indicies) const = 0;
	
	virtual FluidGridType getGridType() const = 0;
	virtual btScalar getCellSize() const = 0;
	
	//for 'reduced' processing
	virtual void removeFirstParticle(int gridCellIndex, const btAlignedObjectArray<int> &nextFluidIndex) = 0;
	
	void getPointAabb(btVector3 *out_pointMin, btVector3 *out_pointMax) const
	{
		*out_pointMin = m_pointMin;
		*out_pointMax = m_pointMax;
	}

	const btAlignedObjectArray<int>& getCellProcessingGroup(int index) const { return m_cellProcessingGroups[index]; }
	
protected:	
	void resetPointAabb()
	{
		m_pointMin.setValue(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);
		m_pointMax.setValue(-BT_LARGE_FLOAT, -BT_LARGE_FLOAT, -BT_LARGE_FLOAT);
	}
	
	void updatePointAabb(const btVector3 &position)
	{
		for(int i = 0; i < 3; ++i)
		{
			if(position.m_floats[i] < m_pointMin.m_floats[i]) m_pointMin.m_floats[i] = position.m_floats[i];
			if(position.m_floats[i] > m_pointMax.m_floats[i]) m_pointMax.m_floats[i] = position.m_floats[i];
		}
	}
	
	//for 'reduced' processing
	virtual int getNumGridCells() const = 0;
	virtual void getIndiciesReduce(int gridCellIndex, int *out_x, int *out_y, int *out_z) const = 0;
	
	void generateCellProcessingGroups()
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
		
		for(int i = 0; i < NUM_CELL_PROCESSING_GROUPS; ++i) m_cellProcessingGroups[i].resize(0);
			
		for(int cell = 0; cell < getNumGridCells(); ++cell)
		{
			FluidGridIterator FI = getGridCell(cell);
			if( !FluidGridIterator::isIndexValid(FI.m_firstIndex, FI.m_lastIndex) ) continue;
		
			int index_x, index_y, index_z;
			getIndiciesReduce(cell, &index_x, &index_y, &index_z);

			//For FluidStaticGrid: Convert range from [0, n] to [1, n+1]
			//
			//For FluidSortingGrid: Convert range from [0, 255] to [1, 256]
			//(SortGridValue::setValue() already converts from [-128, 127], to [0, 255])
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

#define SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED
#ifndef SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED
	//Extracts the (x,y,z) indicies from combinedIndex, where
	//combinedIndex == x + y*resolutionX + z*resolutionX*resolutionY
	static void splitIndex(int resolutionX, int resolutionY, int combinedIndex, int *out_x, int *out_y, int *out_z)
	{
		int x = combinedIndex % resolutionX;
		int z = combinedIndex / (resolutionX*resolutionY);
		int y = (combinedIndex - z*resolutionX*resolutionY) / resolutionX;
				
		*out_x = x;
		*out_z = z;
		*out_y = y;
	}
#else
	static void splitIndex(unsigned long long int resolutionX, unsigned long long int resolutionY, 
						   unsigned long long int combinedIndex, int *out_x, int *out_y, int *out_z)
	{
		unsigned long long int cellsPerLine = resolutionX;
		unsigned long long int cellsPerPlane = resolutionX * resolutionY;
											 
		unsigned long long int x = combinedIndex % cellsPerLine;
		unsigned long long int z = combinedIndex / cellsPerPlane;
		unsigned long long int y = (combinedIndex - z*cellsPerPlane) / cellsPerLine;
				
		*out_x = static_cast<int>(x);
		*out_z = static_cast<int>(z);
		*out_y = static_cast<int>(y);
	}
#endif
	
};


#endif

