/* FluidSortingGridOpenCL.cpp
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
#include "FluidSortingGridOpenCL.h"

#include "../FluidSortingGrid.h"


// /////////////////////////////////////////////////////////////////////////////
//class FluidSortingGrid_OpenCL
// /////////////////////////////////////////////////////////////////////////////
void FluidSortingGrid_OpenCL::writeToOpenCL(cl_context context, cl_command_queue commandQueue, 
											FluidSortingGrid *sortingGrid, bool transferCellProcessingGroups)
{
	int numActiveCells = sortingGrid->getNumGridCells();
	
	if(m_maxActiveCells < numActiveCells)
	{
		deallocate();
		allocate(context, commandQueue, numActiveCells);
	}
	
	btAlignedObjectArray<SortGridValue> &activeCells = sortingGrid->internalGetActiveCells();
	btAlignedObjectArray<FluidGridIterator>&cellContents = sortingGrid->internalGetCellContents();
	
	m_buffer_numActiveCells.writeToBuffer( commandQueue, &numActiveCells, sizeof(int) );
	m_buffer_activeCells.writeToBuffer( commandQueue, &activeCells[0], sizeof(SortGridValue)*numActiveCells );
	m_buffer_cellContents.writeToBuffer( commandQueue, &cellContents[0], sizeof(FluidGridIterator)*numActiveCells );
	
	if(transferCellProcessingGroups)
	{
		for(int i = 0; i < FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS; ++i)
		{
			if(m_cellProcessingGroups[i])
			{
				const btAlignedObjectArray<int> &group = sortingGrid->internalGetCellProcessingGroup(i);
				m_cellProcessingGroups[i]->copyFromHost(group, false);
			}
		}
	}
	
	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}
void FluidSortingGrid_OpenCL::readFromOpenCL(cl_context context, cl_command_queue commandQueue, 
											FluidSortingGrid *sortingGrid, bool transferCellProcessingGroups)
{
	btAssert(m_maxActiveCells != 0);
	
	btAlignedObjectArray<SortGridValue> &activeCells = sortingGrid->internalGetActiveCells();
	btAlignedObjectArray<FluidGridIterator>&cellContents = sortingGrid->internalGetCellContents();
	
	int numActiveCells = getNumActiveCells(commandQueue);
	activeCells.resize(numActiveCells);
	cellContents.resize(numActiveCells);
	
	if(numActiveCells)
	{
		m_buffer_activeCells.readFromBuffer( commandQueue, &activeCells[0], sizeof(SortGridValue)*numActiveCells );
		m_buffer_cellContents.readFromBuffer( commandQueue, &cellContents[0], sizeof(FluidGridIterator)*numActiveCells );	
	
		if(transferCellProcessingGroups)
		{
			for(int i = 0; i < FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS; ++i)
			{
				if(m_cellProcessingGroups[i])
				{
					btAlignedObjectArray<int> &group = sortingGrid->internalGetCellProcessingGroup(i);
					m_cellProcessingGroups[i]->copyToHost(group, false);
				}
			}
		}
		else 
		{
			for(int i = 0; i < FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS; ++i) sortingGrid->internalGetCellProcessingGroup(i).resize(0);
		}
		
		cl_int error_code = clFinish(commandQueue);
		CHECK_CL_ERROR(error_code);
	}
	else 
	{
		for(int i = 0; i < FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS; ++i) sortingGrid->internalGetCellProcessingGroup(i).resize(0);
	}
}

FluidSortingGrid_OpenCLPointers FluidSortingGrid_OpenCL::getPointers()
{
	FluidSortingGrid_OpenCLPointers pointers;
	pointers.m_buffer_numActiveCells = m_buffer_numActiveCells.getAddress();
	pointers.m_buffer_activeCells = m_buffer_activeCells.getAddress();
	pointers.m_buffer_cellContents = m_buffer_cellContents.getAddress();
	
	for(int i = 0; i < FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS; ++i)
	{
		pointers.m_numCellsInGroups[i] = m_cellProcessingGroups[i]->size();
		pointers.m_cellProcessingGroups[i] = m_cellProcessingGroups[i]->getBufferCL();
	}
	
	return pointers;
}

int FluidSortingGrid_OpenCL::getNumActiveCells(cl_command_queue commandQueue)
{
	if(!m_maxActiveCells) return 0;

	int numActiveCells;
	m_buffer_numActiveCells.readFromBuffer( commandQueue, &numActiveCells, sizeof(int) );
	
	cl_int errorCode = clFinish(commandQueue);
	CHECK_CL_ERROR(errorCode);
	
	return numActiveCells;
}

void FluidSortingGrid_OpenCL::resize(cl_context context, cl_command_queue commandQueue, int maxGridCells)
{
	deallocate();
	allocate(context, commandQueue, maxGridCells);
}

void FluidSortingGrid_OpenCL::allocate(cl_context context, cl_command_queue commandQueue, int maxGridCells)
{
	m_maxActiveCells = maxGridCells;
	
	m_buffer_numActiveCells.allocate( context, sizeof(int) );
	m_buffer_activeCells.allocate( context, sizeof(SortGridValue)*maxGridCells );
	m_buffer_cellContents.allocate( context, sizeof(FluidGridIterator)*maxGridCells );
	
	for(int i = 0; i < FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS; ++i)
	{
		if(!m_cellProcessingGroups[i])
		{
			void *ptr = btAlignedAlloc( sizeof(btOpenCLArray<int>), 16 );
			m_cellProcessingGroups[i] = new(ptr) btOpenCLArray<int>(context, commandQueue);
		}
	}
	
	const int NUM_ACTIVE_CELLS = 0;
	m_buffer_numActiveCells.writeToBuffer( commandQueue, &NUM_ACTIVE_CELLS, sizeof(int) );
	
	cl_int errorCode = clFinish(commandQueue);
	CHECK_CL_ERROR(errorCode);
}
void FluidSortingGrid_OpenCL::deallocate()
{
	m_maxActiveCells = 0;
	
	m_buffer_numActiveCells.deallocate();
	m_buffer_activeCells.deallocate();
	m_buffer_cellContents.deallocate();
	
	for(int i = 0; i < FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS; ++i)
	{
		if(m_cellProcessingGroups[i])
		{
			m_cellProcessingGroups[i]->~btOpenCLArray<int>();
			btAlignedFree(m_cellProcessingGroups[i]);
			m_cellProcessingGroups[i] = 0;
		}
	}
}


// /////////////////////////////////////////////////////////////////////////////
//class FluidSortingGrid_OpenCL_Program
// /////////////////////////////////////////////////////////////////////////////
FluidSortingGrid_OpenCL_Program::FluidSortingGrid_OpenCL_Program()
{
	sortingGrid_program = INVALID_PROGRAM;
	kernel_generateValueIndexPairs = INVALID_KERNEL;
	kernel_rearrangeParticleArrays = INVALID_KERNEL;
	kernel_generateUniques = INVALID_KERNEL;
	
	m_radixSorter = 0;
	m_valueIndexPairs = 0;
}

void FluidSortingGrid_OpenCL_Program::initialize(cl_context context, cl_device_id gpu_device, cl_command_queue queue)
{
	cl_int error_code;
	
	//Program
	const char CL_SORTING_GRID_PROGRAM_PATH[] = "./Demos/FluidDemo/Fluids/OpenCL_support/fluids.cl";
	sortingGrid_program = compileProgramOpenCL(context, gpu_device, CL_SORTING_GRID_PROGRAM_PATH);

	//Kernels
	kernel_generateValueIndexPairs = clCreateKernel(sortingGrid_program, "generateValueIndexPairs", &error_code);
	CHECK_CL_ERROR(error_code);
	kernel_rearrangeParticleArrays = clCreateKernel(sortingGrid_program, "rearrangeParticleArrays", &error_code);
	CHECK_CL_ERROR(error_code);
	kernel_generateUniques = clCreateKernel(sortingGrid_program, "generateUniques", &error_code);
	CHECK_CL_ERROR(error_code);
	
	//
	if(!m_radixSorter)
	{
		void *ptr = btAlignedAlloc( sizeof(btRadixSort32CL), 16 );
		m_radixSorter = new(ptr) btRadixSort32CL(context, gpu_device, queue);
		
		ptr = btAlignedAlloc( sizeof(btOpenCLArray<btSortData>), 16 );
		m_valueIndexPairs = new(ptr) btOpenCLArray<btSortData>(context, queue);
	}
}
void FluidSortingGrid_OpenCL_Program::deactivate()
{
	cl_int error_code;

	//Kernels
	error_code = clReleaseKernel(kernel_generateValueIndexPairs);
	CHECK_CL_ERROR(error_code);
	error_code = clReleaseKernel(kernel_rearrangeParticleArrays);
	CHECK_CL_ERROR(error_code);
	error_code = clReleaseKernel(kernel_generateUniques);
	CHECK_CL_ERROR(error_code);
		
	//Program
	error_code = clReleaseProgram(sortingGrid_program);
	CHECK_CL_ERROR(error_code);
	
	//Buffers
	buffer_temp.deallocate();
	
	//
	sortingGrid_program = INVALID_PROGRAM;
	kernel_generateValueIndexPairs = INVALID_KERNEL;
	kernel_rearrangeParticleArrays = INVALID_KERNEL;
	kernel_generateUniques = INVALID_KERNEL;
	
	//
	if(m_radixSorter)
	{
		m_radixSorter->~btRadixSort32CL();
		btAlignedFree(m_radixSorter);
		m_radixSorter = 0;
		
		m_valueIndexPairs->~btOpenCLArray<btSortData>();
		btAlignedFree(m_valueIndexPairs);
		m_valueIndexPairs = 0;
	}
}


void rearrangeToMatchSortedValues2(const btAlignedObjectArray<btSortData> &sortedValues, btAlignedObjectArray<btVector3> &out_rearranged)
{
	static btAlignedObjectArray<btVector3> result;
	result.resize( out_rearranged.size() );
	
	for(int i = 0; i < out_rearranged.size(); ++i)
	{
		int oldIndex = sortedValues[i].m_value;
		int newIndex = i;
			
		result[newIndex] = out_rearranged[oldIndex];
	}
	
	out_rearranged = result;
}
void FluidSortingGrid_OpenCL_Program::insertParticlesIntoGrid(cl_context context, cl_command_queue commandQueue,
															  FluidSph *fluid, Fluid_OpenCL *fluidData, FluidSortingGrid_OpenCL *gridData)
{
	BT_PROFILE("FluidSortingGrid_OpenCL_Program::insertParticlesIntoGrid()");

	btAssert(m_radixSorter);
	btAssert(m_valueIndexPairs);
	cl_int error_code;
	
	int numFluidParticles = fluid->numParticles();
	
	//
	unsigned int fluidParticlesSize = sizeof(btVector3)*numFluidParticles;
	if( buffer_temp.getSize() < fluidParticlesSize )
	{
		buffer_temp.deallocate();
		buffer_temp.allocate(context, fluidParticlesSize);
	}
	
	if( gridData->getMaxActiveCells() < numFluidParticles ) gridData->resize(context, commandQueue, numFluidParticles);
	
	if( m_valueIndexPairs->size() != numFluidParticles ) m_valueIndexPairs->resize(numFluidParticles, true);
	
	//
	Fluid_OpenCLPointers fluidPointers = fluidData->getPointers();
	FluidSortingGrid_OpenCLPointers gridPointers = gridData->getPointers();
	
	//
	{
		BT_PROFILE("generateValueIndexPairs()");
		generateValueIndexPairs( commandQueue, numFluidParticles, fluid->getGrid().getCellSize(), fluidPointers.m_buffer_pos );
		
		error_code = clFinish(commandQueue);
		CHECK_CL_ERROR(error_code);
	}
	
	//Note that btRadixSort32CL uses btSortData, while FluidSortingGrid uses ValueIndexPair.
	//btSortData.m_key == ValueIndexPair.m_value (value to sort by)
	//btSortData.m_value == ValueIndexPair.m_index (fluid particle index)
	{
		BT_PROFILE("radix sort");
		m_radixSorter->execute(*m_valueIndexPairs, 32);
	}
	
	//
	{
		BT_PROFILE("rearrange device");
		
		rearrangeParticleArrays(commandQueue, numFluidParticles, fluidPointers.m_buffer_pos);
		error_code = clEnqueueCopyBuffer( commandQueue, *static_cast<cl_mem*>(buffer_temp.getAddress()), 
											*static_cast<cl_mem*>(fluidPointers.m_buffer_pos), 
										  0, 0, sizeof(btVector3)*numFluidParticles, 0, 0, 0 );
		CHECK_CL_ERROR(error_code);
		
		rearrangeParticleArrays(commandQueue, numFluidParticles, fluidPointers.m_buffer_vel_eval);
		error_code = clEnqueueCopyBuffer( commandQueue, *static_cast<cl_mem*>(buffer_temp.getAddress()), 
											*static_cast<cl_mem*>(fluidPointers.m_buffer_vel_eval), 
											0, 0, sizeof(btVector3)*numFluidParticles, 0, 0, 0 );
		CHECK_CL_ERROR(error_code);
		
		error_code = clFinish(commandQueue);
		CHECK_CL_ERROR(error_code);
	}
	
	{
		BT_PROFILE("rearrange host");
		m_valueIndexPairs->copyToHost(m_valueIndexPairsHost, true);
		
		FluidParticles &particles = fluid->internalGetFluidParticles();
		rearrangeToMatchSortedValues2(m_valueIndexPairsHost, particles.m_pos);
		rearrangeToMatchSortedValues2(m_valueIndexPairsHost, particles.m_vel);
		rearrangeToMatchSortedValues2(m_valueIndexPairsHost, particles.m_vel_eval);
		rearrangeToMatchSortedValues2(m_valueIndexPairsHost, particles.m_externalAcceleration);
	}
	
	//
	{
		BT_PROFILE("generateUniques_serial()");
		generateUniques(commandQueue, numFluidParticles, &gridPointers);
		
		error_code = clFinish(commandQueue);
		CHECK_CL_ERROR(error_code);
	}
}

void FluidSortingGrid_OpenCL_Program::generateValueIndexPairs(cl_command_queue commandQueue, int numFluidParticles, 
															  btScalar cellSize, void *fluidPositionsBufferAddress)
{
	cl_int error_code;

	///btScalar cellSize
	///__global btVector3 *fluidPositions
	///__global ValueIndexPair *out_pairs
	error_code = clSetKernelArg( kernel_generateValueIndexPairs, 0, sizeof(btScalar), &cellSize );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_generateValueIndexPairs, 1, sizeof(void*), fluidPositionsBufferAddress );
	CHECK_CL_ERROR(error_code);
	
	cl_mem valueIndexBuffer = m_valueIndexPairs->getBufferCL();
	error_code = clSetKernelArg( kernel_generateValueIndexPairs, 2, sizeof(void*), static_cast<void*>(&valueIndexBuffer) );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(commandQueue, kernel_generateValueIndexPairs, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
}
void FluidSortingGrid_OpenCL_Program::rearrangeParticleArrays(cl_command_queue commandQueue, int numFluidParticles, void *fluidBufferAddress)
{
	cl_int error_code;

	///__global ValueIndexPair *sortedPairs
	///__global btVector3 *rearrange
	///__global btVector3 *temporary
	cl_mem valueIndexBuffer = m_valueIndexPairs->getBufferCL();
	error_code = clSetKernelArg( kernel_rearrangeParticleArrays, 0, sizeof(void*), static_cast<void*>(&valueIndexBuffer) );
	CHECK_CL_ERROR(error_code);
	
	error_code = clSetKernelArg( kernel_rearrangeParticleArrays, 1, sizeof(void*), fluidBufferAddress );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_rearrangeParticleArrays, 2, sizeof(void*), buffer_temp.getAddress() );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(commandQueue, kernel_rearrangeParticleArrays, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
}

void FluidSortingGrid_OpenCL_Program::generateUniques(cl_command_queue commandQueue, int numFluidParticles, 
													  FluidSortingGrid_OpenCLPointers *gridPointers)
{
	cl_int error_code;

	///__global ValueIndexPair *sortedPairs
	///int numSortedPairs
	///__global SortGridValue *out_activeCells
	///__global FluidGridIterator *out_cellContents
	///__global int *out_numActiveCells
	cl_mem valueIndexBuffer = m_valueIndexPairs->getBufferCL();
	error_code = clSetKernelArg( kernel_generateUniques, 0, sizeof(void*),  static_cast<void*>(&valueIndexBuffer) );
	CHECK_CL_ERROR(error_code);
	
	error_code = clSetKernelArg( kernel_generateUniques, 1, sizeof(int), &numFluidParticles );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_generateUniques, 2, sizeof(void*), gridPointers->m_buffer_activeCells );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_generateUniques, 3, sizeof(void*), gridPointers->m_buffer_cellContents );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_generateUniques, 4, sizeof(void*), gridPointers->m_buffer_numActiveCells );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = 1;
	error_code = clEnqueueNDRangeKernel(commandQueue, kernel_generateUniques, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
}
