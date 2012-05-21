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
#ifndef DEF_GRID_H_INCLUDED
#define DEF_GRID_H_INCLUDED

#include <vector>

#include "vector3df.h"	
#include "LinearMath/btVector3.h"
//#include "LinearMath/btAlignedObjectArray.h"


const int INVALID_PARTICLE_INDEX = -1;

const int RESULTS_PER_GRID_SEARCH = 8;		//Number of grid cell indicies returned from Grid::findCells()
struct GridCellIndicies { int m_indicies[RESULTS_PER_GRID_SEARCH]; };

struct GridParameters
{
	Vector3DF	m_min;						//Volume of grid (may not match domain volume exactly)
	Vector3DF	m_max;
	
	int			m_resolutionX;				//Number of cells per axis
	int			m_resolutionY;
	int			m_resolutionZ;
	
	int			m_numCells;					//Total number of cells
	float		m_gridCellSize;				//Edge length of a cube-shaped cell
};

struct Fluid;
class Grid
{
	std::vector<int>			m_grid;						//Contains the index of the last added particle in a forward linked list
	std::vector<int>			m_gridNumFluids;			//Contains the number of 'struct Fluid'(s) in the cell
	
	GridParameters				m_params;
	
public:
	void setup(const Vector3DF &min, const Vector3DF &max, float simScale, float cellSize, float border);
	
	void clear();
	void insertParticle(Fluid *p, int particleIndex);
	
	void findCells(const Vector3DF &position, float radius, GridCellIndicies *out_findCellsResult) const;
	int getLastParticleIndex(int gridCellIndex) const { return m_grid[gridCellIndex]; }
	
	const GridParameters& getParameters() const { return m_params; } 
	
	//for OpenCL
	int* getCellsPointer() { return &m_grid[0]; }
	int* getCellsNumFluidsPointer() { return &m_gridNumFluids[0]; }
	
	inline void getIndicies(const Vector3DF &position, int *out_index_x, int *out_index_y, int *out_index_z) const
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
		/*
		unsigned int index_x = static_cast<unsigned int>( (position.x() - m_params.m_min.x()) / m_params.m_gridCellSize );
		unsigned int index_y = static_cast<unsigned int>( (position.y() - m_params.m_min.y()) / m_params.m_gridCellSize );
		unsigned int index_z = static_cast<unsigned int>( (position.z() - m_params.m_min.z()) / m_params.m_gridCellSize );
		*out_index_x = (index_x < m_params.m_resolutionX) ? index_x : m_params.m_resolutionX - 1;
		*out_index_y = (index_y < m_params.m_resolutionY) ? index_y : m_params.m_resolutionY - 1;
		*out_index_z = (index_z < m_params.m_resolutionZ) ? index_z : m_params.m_resolutionZ - 1;
		*/
	}
};

#endif
