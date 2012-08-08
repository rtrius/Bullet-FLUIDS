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
#ifndef FLUID_STATIC_GRID_H
#define FLUID_STATIC_GRID_H

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"

#include "FluidGrid.h"

struct FluidParticles;

struct FluidStaticGridParameters
{
	btVector3	m_min;						//Volume of grid (may not match domain volume exactly)
	btVector3	m_max;
	btScalar	m_gridCellSize;				//Edge length of a cube-shaped cell
	
	int			m_resolutionX;				//Number of cells per axis
	int			m_resolutionY;
	int			m_resolutionZ;
	
	int			m_numCells;					//Total number of cells
};

class FluidStaticGrid : public FluidGrid
{
	//out_gridCells->m_iterators[i].m_lastIndex = LAST_INDEX -- see FluidGridIterator::isIndexValid()
	static const int LAST_INDEX = 2147483647;	//2^31 - 1; Value greater than the highest fluid particle index

	btAlignedObjectArray<int>	m_grid;				//Contains the index of the last added particle in a forward linked list
	btAlignedObjectArray<int>	m_gridNumFluids;	//Contains the number of 'struct Fluid'(s) in the cell
	
	FluidStaticGridParameters	m_params;
	
public:
	FluidStaticGrid() { m_params.m_numCells = 0; }		//Internal constructor; for FluidSolverReducedGridNeighbor
	FluidStaticGrid(const btVector3 &min, const btVector3 &max, btScalar simScale, btScalar cellSize, btScalar border) 
	{ 
		setup(min, max, simScale, cellSize, border);
	}

	void setup(const btVector3 &min, const btVector3 &max, btScalar simScale, btScalar cellSize, btScalar border);
	
	virtual void clear();	
	virtual void insertParticles(FluidParticles *fluids);
	virtual void findCells(const btVector3 &position, btScalar radius, FindCellsResult *out_gridCells) const;
	
	virtual FluidGridIterator getGridCell(int gridCellIndex) const { return FluidGridIterator(m_grid[gridCellIndex], LAST_INDEX); }
	virtual void getGridCellIndiciesInAabb(const btVector3 &min, const btVector3 &max, btAlignedObjectArray<int> *out_indicies) const;
	
	virtual FluidGridType getGridType() const { return FT_LinkedList; }
	virtual btScalar getCellSize() const { return m_params.m_gridCellSize; }
	
	virtual int getNumGridCells() const { return m_params.m_numCells; }
	virtual void getIndiciesReduce(int gridCellIndex, int *out_x, int *out_y, int *out_z) const
	{
		splitIndex(m_params.m_resolutionX, m_params.m_resolutionY, gridCellIndex, out_x, out_y, out_z);
	}
	virtual void removeFirstParticle(int gridCellIndex, const btAlignedObjectArray<int> &nextFluidIndex)
	{
		if(m_grid[gridCellIndex] != INVALID_PARTICLE_INDEX) m_grid[gridCellIndex] = nextFluidIndex[ m_grid[gridCellIndex] ];
	}
	
	//for OpenCL
	const FluidStaticGridParameters& getParameters() const { return m_params; } 
	int* getCellsPointer() { return &m_grid[0]; }
	int* getCellsNumFluidsPointer() { return &m_gridNumFluids[0]; }
	
private:
	void getIndicies(const btVector3 &position, int *out_index_x, int *out_index_y, int *out_index_z) const;
	void insertParticle(const btVector3 &position, int particleIndex, int *fluidNextIndex);
};

#endif
