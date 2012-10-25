/* FluidSortingGridOpenCL.h
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
#ifndef FLUID_SORTING_GRID_OPENCL_H
#define FLUID_SORTING_GRID_OPENCL_H

#include "LinearMath/btScalar.h"

#include "btExperimentsOpenCL/btRadixSort32CL.h"

#include "opencl_support.h"
#include "FluidOpenCL.h"
#include "../FluidSph.h"
#include "../FluidSortingGrid.h"

class FluidSortingGrid;
class Fluid_OpenCL;

//for FluidSolverOpenCLSymmetric / fluidsSymmetric.cl
#define NUM_FOUND_CELLS_SYMMETRIC_CL 14
typedef struct
{
	FluidGridIterator m_iterators[NUM_FOUND_CELLS_SYMMETRIC_CL];
	
} FoundCellsSymmetric;


///@brief Manages OpenCL buffers corresponding to a FluidSortingGrid.
class FluidSortingGrid_OpenCL
{
	int m_maxActiveCells;

public:	
	OpenCLBuffer m_buffer_numActiveCells;	//int
	OpenCLBuffer m_buffer_activeCells;		//SortGridValue[]
	OpenCLBuffer m_buffer_cellContents;		//FluidGridIterator[]
	
	btOpenCLArray<int> *m_cellProcessingGroups[FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS];
	btOpenCLArray<FoundCellsSymmetric> *m_adjacentCells;
	
	FluidSortingGrid_OpenCL() : m_maxActiveCells(0), m_adjacentCells(0)
	{ 
		for(int i = 0; i < FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS; ++i) m_cellProcessingGroups[i] = 0; 
	}
	~FluidSortingGrid_OpenCL() { deallocate(); }
	
	void writeToOpenCL(cl_context context, cl_command_queue commandQueue, FluidSortingGrid *sortingGrid, bool transferCellProcessingGroups);
	void readFromOpenCL(cl_context context, cl_command_queue commandQueue, FluidSortingGrid *sortingGrid, bool transferCellProcessingGroups);
	
	int getNumActiveCells(cl_command_queue commandQueue);
	
	int getMaxActiveCells() const { return m_maxActiveCells; }
	void resize(cl_context context, cl_command_queue commandQueue, int maxGridCells);
	
private:
	void allocate(cl_context context, cl_command_queue commandQueue, int maxGridCells);
	void deallocate();
};

class FluidSortingGrid_OpenCL_Program
{
	OpenCLBuffer buffer_temp;				//btVector3[] -- used to rearrange fluid particle arrays(position, velocity, etc.)

	cl_program sortingGrid_program;
	cl_kernel kernel_generateValueIndexPairs;
	cl_kernel kernel_rearrangeParticleArrays;
	cl_kernel kernel_generateUniques;

	btRadixSort32CL *m_radixSorter;
	btOpenCLArray<btSortData> *m_valueIndexPairs;
	btAlignedObjectArray<btSortData> m_valueIndexPairsHost;
	
public:
	FluidSortingGrid_OpenCL_Program();
	~FluidSortingGrid_OpenCL_Program() { deactivate(); }

	void initialize(cl_context context, cl_device_id gpu_device, cl_command_queue queue);
	void deactivate();

	void insertParticlesIntoGrid(cl_context context, cl_command_queue commandQueue,
								 FluidSph *fluid, Fluid_OpenCL *fluidData, FluidSortingGrid_OpenCL *gridData);
private:
	
	void generateValueIndexPairs(cl_command_queue commandQueue, int numFluidParticles, btScalar cellSize, cl_mem fluidPositionsBuffer);
	void rearrangeParticleArrays(cl_command_queue commandQueue, int numFluidParticles, cl_mem fluidBuffer);
	void generateUniques(cl_command_queue commandQueue, int numFluidParticles, FluidSortingGrid_OpenCL *gridData);
};

#endif