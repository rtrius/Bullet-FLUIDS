/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. 
   If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_FLUID_SORTING_GRID_OPENCL_H
#define BT_FLUID_SORTING_GRID_OPENCL_H

#include "LinearMath/btScalar.h"

#include "btExperimentsOpenCL/btRadixSort32CL.h"

#include "opencl_support.h"
#include "../btFluidSph.h"
#include "../btFluidSortingGrid.h"

class btFluidSortingGrid;
class btFluidSphOpenCL;

///Manages OpenCL buffers corresponding to a btFluidSortingGrid.
class btFluidSortingGridOpenCL
{
public:	
	btOpenCLArray<int> m_numActiveCells;
	
	btOpenCLArray<btSortGridValue> m_activeCells;
	btOpenCLArray<btFluidGridIterator> m_cellContents;

	btFluidSortingGridOpenCL(cl_context context, cl_command_queue queue) 
		: m_numActiveCells(context, queue), m_activeCells(context, queue), m_cellContents(context, queue) {}
		
	void writeToOpenCL(cl_command_queue queue, btFluidSortingGrid& sortingGrid);
	void readFromOpenCL(cl_command_queue queue, btFluidSortingGrid& sortingGrid);
	
	int getNumActiveCells() const;
};

///Implements btFluidSortingGrid::insertParticles() for OpenCL.
class btFluidSortingGridOpenCLProgram
{
	cl_program sortingGrid_program;
	cl_kernel kernel_generateValueIndexPairs;
	cl_kernel kernel_rearrangeParticleArrays;
	cl_kernel kernel_generateUniques;

	btOpenCLArray<btVector3> m_tempBufferCL;		//Used to rearrange fluid particle arrays(position, velocity, etc.)
	btAlignedObjectArray<btVector3> m_tempBuffer;
	
	btRadixSort32CL m_radixSorter;
	btOpenCLArray<btSortData> m_valueIndexPairs;
	btAlignedObjectArray<btSortData> m_valueIndexPairsHost;
	
public:
	btFluidSortingGridOpenCLProgram(cl_context context, cl_command_queue queue, cl_device_id device);
	~btFluidSortingGridOpenCLProgram();
	
	void insertParticlesIntoGrid(cl_context context, cl_command_queue commandQueue,
								 btFluidSph* fluid, btFluidSphOpenCL* fluidData, btFluidSortingGridOpenCL* gridData);
								 
private:
	void generateValueIndexPairs(cl_command_queue commandQueue, int numFluidParticles, btScalar cellSize, cl_mem fluidPositionsBuffer);
	void rearrangeParticleArrays(cl_command_queue commandQueue, int numFluidParticles, cl_mem fluidBuffer);
	void generateUniques(cl_command_queue commandQueue, int numFluidParticles, btFluidSortingGridOpenCL* gridData);
};

#endif