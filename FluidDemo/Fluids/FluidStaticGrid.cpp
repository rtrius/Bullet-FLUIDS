/*
  FLUIDS v.1 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

  ZLib license
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
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

#include "FluidStaticGrid.h"

#include "FluidParticles.h"	//for INVALID_PARTICLE_INDEX

void FluidStaticGrid::setup(const btVector3 &min, const btVector3 &max, btScalar simScale, btScalar simCellSize, btScalar border)
{
	btScalar worldCellSize = simCellSize / simScale;
	m_params.m_gridCellSize = worldCellSize;

	m_params.m_min.setValue( min.x() - border, min.y() - border, min.z() - border );
	m_params.m_max.setValue( max.x() + border, max.y() + border, max.z() + border );

	m_params.m_resolutionX = static_cast<int>( ceil((max.x() - min.x()) / worldCellSize) );
	m_params.m_resolutionY = static_cast<int>( ceil((max.y() - min.y()) / worldCellSize) );
	m_params.m_resolutionZ = static_cast<int>( ceil((max.z() - min.z()) / worldCellSize) );
	
	m_params.m_numCells = m_params.m_resolutionX * m_params.m_resolutionY * m_params.m_resolutionZ;

	m_grid.resize(m_params.m_numCells);
	m_gridNumFluids.resize(m_params.m_numCells);
	for(int n = 0; n < m_params.m_numCells; n++)
	{
		m_grid[n] = INVALID_PARTICLE_INDEX;
		m_gridNumFluids[n] = 0;
	}
}

void FluidStaticGrid::clear()
{
	for(int n = 0; n < m_params.m_numCells; n++)
	{
		m_grid[n] = INVALID_PARTICLE_INDEX;
		m_gridNumFluids[n] = 0;
	}
}

void FluidStaticGrid::insertParticles(FluidParticles *fluids)
{
	resetPointAabb();
	for(int i = 0; i < fluids->size(); ++i)
	{
		updatePointAabb(fluids->m_pos[i]);
		insertParticle(fluids->m_pos[i], i, &fluids->m_nextFluidIndex[i]);
	}
		
	generateCellProcessingGroups();
}

void FluidStaticGrid::findCells(const btVector3 &position, btScalar radius, FluidGrid::FoundCells *out_gridCells) const
{
	//Store a 2x2x2 grid cell query result in out_gridCells,
	//where out_gridCells.m_iterators[0], the cell with the lowest index,
	//corresponds to the minimum point of the sphere's AABB
	
	//Determine the grid cell index at the minimum point of the particle's AABB
	int index_x = static_cast<int>( (-radius + position.x() - m_params.m_min.x()) / m_params.m_gridCellSize );
	int index_y = static_cast<int>( (-radius + position.y() - m_params.m_min.y()) / m_params.m_gridCellSize );
	int index_z = static_cast<int>( (-radius + position.z() - m_params.m_min.z()) / m_params.m_gridCellSize );
	
	//Clamp index to grid bounds
	if(index_x < 0) index_x = 0;
	if(index_y < 0) index_y = 0;
	if(index_z < 0) index_z = 0;
	
		//Since a 2x2x2 volume is accessed, subtract 2 from the upper index bounds
			//Subtract 1 as a 2x2x2 volume is accessed, and the index we want is the 'min' index
			//Subtract 1 again as indicies start from 0 (m_params.m_resolutionX/Y/Z is out of bounds)
	if(index_x >= m_params.m_resolutionX - 2) index_x = m_params.m_resolutionX - 2;
	if(index_y >= m_params.m_resolutionY - 2) index_y = m_params.m_resolutionY - 2;
	if(index_z >= m_params.m_resolutionZ - 2) index_z = m_params.m_resolutionZ - 2;
	
	//Calculate grid cell indicies
	const int stride_x = 1;
	const int stride_y = m_params.m_resolutionX;
	const int stride_z = m_params.m_resolutionX*m_params.m_resolutionY;
	
	int gridCellIndicies[8];
	gridCellIndicies[0] = (index_z * m_params.m_resolutionY + index_y) * m_params.m_resolutionX + index_x;
	gridCellIndicies[1] = gridCellIndicies[0] + stride_x;
	gridCellIndicies[2] = gridCellIndicies[0] + stride_y;
	gridCellIndicies[3] = gridCellIndicies[0] + stride_y + stride_x;

	gridCellIndicies[4] = gridCellIndicies[0] + stride_z;
	gridCellIndicies[5] = gridCellIndicies[1] + stride_z;
	gridCellIndicies[6] = gridCellIndicies[2] + stride_z;
	gridCellIndicies[7] = gridCellIndicies[3] + stride_z;
	
	for(int i = 0; i < 8; ++i) 
	{
		out_gridCells->m_iterators[i].m_firstIndex = m_grid[ gridCellIndicies[i] ];
		out_gridCells->m_iterators[i].m_lastIndex = LAST_INDEX;
	}
}

void FluidStaticGrid::getGridCellIndiciesInAabb(const btVector3 &min, const btVector3 &max, btAlignedObjectArray<int> *out_indicies) const
{
	int minX, minY, minZ;
	int maxX, maxY, maxZ;
	getIndicies(min, &minX, &minY, &minZ);
	getIndicies(max, &maxX, &maxY, &maxZ);
	
	for(int z = minZ; z <= maxZ; ++z)
		for(int y = minY; y <= maxY; ++y)
			for(int x = minX; x <= maxX; ++x)
			{
				int gridCellIndex = (z*m_params.m_resolutionY + y)*m_params.m_resolutionX + x;
				out_indicies->push_back(gridCellIndex);
			}
}
	
void FluidStaticGrid::getIndicies(const btVector3 &position, int *out_index_x, int *out_index_y, int *out_index_z) const
{
	int index_x = static_cast<int>( (position.x() - m_params.m_min.x()) / m_params.m_gridCellSize );
	int index_y = static_cast<int>( (position.y() - m_params.m_min.y()) / m_params.m_gridCellSize );
	int index_z = static_cast<int>( (position.z() - m_params.m_min.z()) / m_params.m_gridCellSize );
	
	if(index_x < 0) index_x = 0;
	if(index_y < 0) index_y = 0;
	if(index_z < 0) index_z = 0;

	if(index_x > m_params.m_resolutionX - 1) index_x = m_params.m_resolutionX - 1;
	if(index_y > m_params.m_resolutionY - 1) index_y = m_params.m_resolutionY - 1;
	if(index_z > m_params.m_resolutionZ - 1) index_z = m_params.m_resolutionZ - 1;

	*out_index_x = index_x;
	*out_index_y = index_y;
	*out_index_z = index_z;
	
	//unsigned int index_x = static_cast<unsigned int>( (position.x() - m_params.m_min.x()) / m_params.m_gridCellSize );
	//unsigned int index_y = static_cast<unsigned int>( (position.y() - m_params.m_min.y()) / m_params.m_gridCellSize );
	//unsigned int index_z = static_cast<unsigned int>( (position.z() - m_params.m_min.z()) / m_params.m_gridCellSize );
	//*out_index_x = (index_x < m_params.m_resolutionX) ? index_x : m_params.m_resolutionX - 1;
	//*out_index_y = (index_y < m_params.m_resolutionY) ? index_y : m_params.m_resolutionY - 1;
	//*out_index_z = (index_z < m_params.m_resolutionZ) ? index_z : m_params.m_resolutionZ - 1;
}

void FluidStaticGrid::insertParticle(const btVector3 &position, int particleIndex, int *fluidNextIndex)
{
	int index_x = static_cast<int>( (position.x() - m_params.m_min.x()) / m_params.m_gridCellSize );
	int index_y = static_cast<int>( (position.y() - m_params.m_min.y()) / m_params.m_gridCellSize );
	int index_z = static_cast<int>( (position.z() - m_params.m_min.z()) / m_params.m_gridCellSize );
	
	int cellIndex = (index_z*m_params.m_resolutionY + index_y)*m_params.m_resolutionX + index_x;
	
	if(0 <= cellIndex && cellIndex < m_params.m_numCells) 
	{
		//Add particle to linked list
		*fluidNextIndex = m_grid[cellIndex];
		m_grid[cellIndex] = particleIndex;
		
		//
		m_gridNumFluids[cellIndex]++;
	}
}
