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
#include "../FluidSph.h"
#include "../FluidSortingGrid.h"

class FluidSortingGrid;
class FluidSphOpenCL;

///Manages OpenCL buffers corresponding to a FluidSortingGrid.
class FluidSortingGridOpenCL
{
public:	
	btOpenCLArray<int> m_numActiveCells;
	
	btOpenCLArray<SortGridValue> m_activeCells;
	btOpenCLArray<FluidGridIterator> m_cellContents;

	FluidSortingGridOpenCL(cl_context context, cl_command_queue queue) 
		: m_numActiveCells(context, queue), m_activeCells(context, queue), m_cellContents(context, queue) {}
		
	void writeToOpenCL(cl_command_queue queue, FluidSortingGrid *sortingGrid);
	void readFromOpenCL(cl_command_queue queue, FluidSortingGrid *sortingGrid);
	
	int getNumActiveCells() const;
};

///Implements FluidSortingGrid::insertParticles() for OpenCL.
class FluidSortingGridOpenCLProgram
{
	cl_program sortingGrid_program;
	cl_kernel kernel_generateValueIndexPairs;
	cl_kernel kernel_rearrangeParticleArrays;
	cl_kernel kernel_generateUniques;

	btOpenCLArray<btVector3> m_tempBuffer;		//Used to rearrange fluid particle arrays(position, velocity, etc.)
	
	btRadixSort32CL m_radixSorter;
	btOpenCLArray<btSortData> m_valueIndexPairs;
	btAlignedObjectArray<btSortData> m_valueIndexPairsHost;
	
public:
	FluidSortingGridOpenCLProgram(cl_context context, cl_command_queue queue, cl_device_id device);
	~FluidSortingGridOpenCLProgram();
	
	void insertParticlesIntoGrid(cl_context context, cl_command_queue commandQueue,
								 FluidSph *fluid, FluidSphOpenCL *fluidData, FluidSortingGridOpenCL *gridData);
								 
private:
	void generateValueIndexPairs(cl_command_queue commandQueue, int numFluidParticles, btScalar cellSize, cl_mem fluidPositionsBuffer);
	void rearrangeParticleArrays(cl_command_queue commandQueue, int numFluidParticles, cl_mem fluidBuffer);
	void generateUniques(cl_command_queue commandQueue, int numFluidParticles, FluidSortingGridOpenCL *gridData);
};

#endif