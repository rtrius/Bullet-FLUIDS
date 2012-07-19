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

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"

const int RESULTS_PER_GRID_SEARCH = 8;		//Number of grid cell indicies returned from Grid::findCells()
struct GridCellIndicies { int m_indicies[RESULTS_PER_GRID_SEARCH]; };

struct GridParameters
{
	btVector3	m_min;						//Volume of grid (may not match domain volume exactly)
	btVector3	m_max;
	btScalar	m_gridCellSize;				//Edge length of a cube-shaped cell
	
	int			m_resolutionX;				//Number of cells per axis
	int			m_resolutionY;
	int			m_resolutionZ;
	
	int			m_numCells;					//Total number of cells
};

class Grid
{
	btAlignedObjectArray<int>	m_grid;				//Contains the index of the last added particle in a forward linked list
	btAlignedObjectArray<int>	m_gridNumFluids;	//Contains the number of 'struct Fluid'(s) in the cell
	
	GridParameters				m_params;
	
public:
	Grid() { m_params.m_numCells = 0; }

	void setup(const btVector3 &min, const btVector3 &max, btScalar simScale, btScalar cellSize, btScalar border);
	
	void clear();
	void insertParticle(const btVector3 &position, int particleIndex, int *fluidNextIndex);
	
	void findCells(const btVector3 &position, btScalar radius, GridCellIndicies *out_findCellsResult) const;
	int getLastParticleIndex(int gridCellIndex) const { return m_grid[gridCellIndex]; }
	void setLastParticleIndex(int gridCellIndex, int cellValue) { m_grid[gridCellIndex] = cellValue; }
	
	const GridParameters& getParameters() const { return m_params; } 
	
	void getIndicies(const btVector3 &position, int *out_index_x, int *out_index_y, int *out_index_z) const;
	
	//for OpenCL
	int* getCellsPointer() { return &m_grid[0]; }
	int* getCellsNumFluidsPointer() { return &m_gridNumFluids[0]; }
};

#endif
