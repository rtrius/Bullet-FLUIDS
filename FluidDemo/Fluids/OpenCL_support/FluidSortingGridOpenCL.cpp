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
void FluidSortingGrid_OpenCL::writeToOpenCL(cl_context context, cl_command_queue commandQueue, FluidSortingGrid *sortingGrid)
{
	int numActiveCells = sortingGrid->getNumGridCells();
	
	if(m_maxActiveCells < numActiveCells)
	{
		deallocate();
		allocate(context, numActiveCells);
	}
	
	btAlignedObjectArray<SortGridValue> &activeCells = sortingGrid->internalGetActiveCells();
	btAlignedObjectArray<FluidGridIterator>&cellContents = sortingGrid->internalGetCellContents();
	
	m_buffer_numActiveCells.writeToBuffer( commandQueue, &numActiveCells, sizeof(int) );
	m_buffer_activeCells.writeToBuffer( commandQueue, &activeCells[0], sizeof(SortGridValue)*numActiveCells );
	m_buffer_cellContents.writeToBuffer( commandQueue, &cellContents[0], sizeof(FluidGridIterator)*numActiveCells );
	
	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}
void FluidSortingGrid_OpenCL::readFromOpenCL(cl_context context, cl_command_queue commandQueue, FluidSortingGrid *sortingGrid)
{
	btAssert(m_maxActiveCells != 0);
	
	btAlignedObjectArray<SortGridValue> &activeCells = sortingGrid->internalGetActiveCells();
	btAlignedObjectArray<FluidGridIterator>&cellContents = sortingGrid->internalGetCellContents();
	
	int numActiveCells;
	m_buffer_numActiveCells.readFromBuffer( commandQueue, &numActiveCells, sizeof(int) );
	activeCells.resize(numActiveCells);
	cellContents.resize(numActiveCells);
	
	m_buffer_activeCells.readFromBuffer( commandQueue, &activeCells[0], sizeof(SortGridValue)*numActiveCells );
	m_buffer_cellContents.readFromBuffer( commandQueue, &cellContents[0], sizeof(FluidGridIterator)*numActiveCells );	
	
	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}

FluidSortingGrid_OpenCLPointers FluidSortingGrid_OpenCL::getPointers()
{
	FluidSortingGrid_OpenCLPointers pointers;
	pointers.m_buffer_numActiveCells = m_buffer_numActiveCells.getAddress();
	pointers.m_buffer_activeCells = m_buffer_activeCells.getAddress();
	pointers.m_buffer_cellContents = m_buffer_cellContents.getAddress();
	
	return pointers;
}

int FluidSortingGrid_OpenCL::getNumActiveCells(cl_command_queue commandQueue)
{
	btAssert(m_maxActiveCells != 0);
	
	int numActiveCells;
	m_buffer_numActiveCells.readFromBuffer( commandQueue, &numActiveCells, sizeof(int) );
	
	cl_int errorCode = clFinish(commandQueue);
	CHECK_CL_ERROR(errorCode);
	
	return numActiveCells;
}

void FluidSortingGrid_OpenCL::resize(cl_context context, int maxGridCells)
{
	deallocate();
	allocate(context, maxGridCells);
}

void FluidSortingGrid_OpenCL::allocate(cl_context context, int maxGridCells)
{
	m_maxActiveCells = maxGridCells;
	
	m_buffer_numActiveCells.allocate( context, sizeof(int) );
	m_buffer_activeCells.allocate( context, sizeof(SortGridValue)*maxGridCells );
	m_buffer_cellContents.allocate( context, sizeof(FluidGridIterator)*maxGridCells );
}
void FluidSortingGrid_OpenCL::deallocate()
{
	m_maxActiveCells = 0;
	
	m_buffer_numActiveCells.deallocate();
	m_buffer_activeCells.deallocate();
	m_buffer_cellContents.deallocate();
}


// /////////////////////////////////////////////////////////////////////////////
//class FluidSortingGrid_OpenCL_Program
// /////////////////////////////////////////////////////////////////////////////
/*
FluidSortingGrid_OpenCL_Program::FluidSortingGrid_OpenCL_Program()
{
	sortingGrid_program = INVALID_PROGRAM;
	kernel_generateValueIndexPairs = INVALID_KERNEL;
	kernel_rearrangeParticleArrays = INVALID_KERNEL;
	kernel_rearrangeParticleArraysWriteBack = INVALID_KERNEL;
	kernel_generateUniques = INVALID_KERNEL;
}
//~FluidSortingGrid_OpenCL_Program() { buffer_temp.deallocate(); }

void FluidSortingGrid_OpenCL_Program::initialize(cl_context context, cl_device_id gpu_device)
{
	cl_int error_code;
	
	//Program
	const char CL_SORTING_GRID_PROGRAM_PATH[] = "INVALID_PATH INVALID_PATH ./Demos/FluidDemo/Fluids/OpenCL_support/fluidsSortingGrid.cl";
	sortingGrid_program = compileProgramOpenCL(context, gpu_device, CL_SORTING_GRID_PROGRAM_PATH);

	//Kernels
	kernel_generateValueIndexPairs = clCreateKernel(sortingGrid_program, "generateValueIndexPairs", &error_code);
	CHECK_CL_ERROR(error_code);
	kernel_rearrangeParticleArrays = clCreateKernel(sortingGrid_program, "rearrangeParticleArrays", &error_code);
	CHECK_CL_ERROR(error_code);
	kernel_rearrangeParticleArraysWriteBack = clCreateKernel(sortingGrid_program, "rearrangeParticleArraysWriteBack", &error_code);
	CHECK_CL_ERROR(error_code);
	kernel_generateUniques = clCreateKernel(sortingGrid_program, "generateUniques", &error_code);
	CHECK_CL_ERROR(error_code);
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
	buffer_pairs.deallocate();
	buffer_temp.deallocate();
	
	//
	sortingGrid_program = INVALID_PROGRAM;
	kernel_generateValueIndexPairs = INVALID_KERNEL;
	kernel_rearrangeParticleArrays = INVALID_KERNEL;
	kernel_rearrangeParticleArraysWriteBack = INVALID_KERNEL;
	kernel_generateUniques = INVALID_KERNEL;
}

void FluidSortingGrid_OpenCL_Program::insertParticlesIntoGrids(cl_context context, cl_command_queue commandQueue, int maxFluidParticles, 
														  	   btAlignedObjectArray<FluidSph*> *fluids, 
														 	   btAlignedObjectArray<Fluid_OpenCL*> *fluidData,
															   btAlignedObjectArray<FluidSortingGrid_OpenCL*> *gridData)
{
	unsigned int maxFluidParticlesSize = sizeof(btVector3)*maxFluidParticles;
	if( buffer_temp.getSize() < maxFluidParticlesSize )
	{
		buffer_temp.deallocate();
		buffer_temp.allocate(context, maxFluidParticlesSize);
	}
	
	unsigned int maxPairsSize = sizeof(ValueIndexPair)*maxFluidParticles;
	if( buffer_pairs.getSize() < maxPairsSize )
	{
		buffer_pairs.deallocate();
		buffer_pairs.allocate(context, maxPairsSize);
	}

	for(int i = 0; i < fluids->size(); ++i) insertParticlesSingleFluid( context, commandQueue, (*fluids)[i], (*fluidData)[i], (*gridData)[i] );
}

void FluidSortingGrid_OpenCL_Program::insertParticlesSingleFluid(cl_context context, cl_command_queue commandQueue,
																 FluidSph *fluid, Fluid_OpenCL *fluidData, FluidSortingGrid_OpenCL *gridData)
{
	btAssert( fluid->getGrid()->getGridType() == FluidGrid::FT_IndexRange );

	int numFluidParticles = fluid->numParticles();

	Fluid_OpenCLPointers fluidPointers = fluidData->getPointers();
	FluidSortingGrid_OpenCLPointers gridPointers = gridData->getPointers();
	
	//
	generateValueIndexPairs(commandQueue, numFluidParticles, fluidPointers.m_buffer_pos);
	
	//
	//RADIX_SORTER::RADIX_SORT();
	
	//
	rearrangeParticleArrays(commandQueue, numFluidParticles, fluidPointers.m_buffer_pos);
	rearrangeParticleArraysWriteBack(commandQueue, numFluidParticles, fluidPointers.m_buffer_pos);
	rearrangeParticleArrays(commandQueue, numFluidParticles, fluidPointers.m_buffer_vel);
	rearrangeParticleArraysWriteBack(commandQueue, numFluidParticles, fluidPointers.m_buffer_vel);
	rearrangeParticleArrays(commandQueue, numFluidParticles, fluidPointers.m_buffer_vel_eval);
	rearrangeParticleArraysWriteBack(commandQueue, numFluidParticles, fluidPointers.m_buffer_vel_eval);
	rearrangeParticleArrays(commandQueue, numFluidParticles, fluidPointers.m_buffer_externalAcceleration);
	rearrangeParticleArraysWriteBack(commandQueue, numFluidParticles, fluidPointers.m_buffer_externalAcceleration);
	
	//
	generateUniques(commandQueue, numFluidParticles, &gridPointers);
	
		//generateUniques() returns the number of active cells in FluidSortingGrid_OpenCL.m_buffer_numActiveCells
	int numActiveCells = gridData->getNumActiveCells(commandQueue);
	if( numActiveCells < gridData->getMaxActiveCells() )
	{
		gridData->resize(context, numActiveCells);
		gridPointers = gridData->getPointers();
		
		generateUniques(commandQueue, numFluidParticles, &gridPointers);
	}
	
	//Writeback
	FluidSortingGrid *sortingGrid = reinterpret_cast<FluidSortingGrid*>( fluid->internalGetGrid() );
	gridData->readFromOpenCL(context, commandQueue, sortingGrid);
	
	//cl_int error_code = clFinish(commandQueue);
	//CHECK_CL_ERROR(error_code);
}

void FluidSortingGrid_OpenCL_Program::generateValueIndexPairs(cl_command_queue commandQueue, int numFluidParticles, 
															  void *fluidPositionsBufferAddress)
{
	cl_int error_code;

	///__global btVector3 *fluidPositions
	///__global ValueIndexPair *out_pairs
	error_code = clSetKernelArg( kernel_generateValueIndexPairs, 0, sizeof(void*), fluidPositionsBufferAddress );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_generateValueIndexPairs, 1, sizeof(void*), buffer_pairs.getAddress() );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(commandQueue, kernel_generateValueIndexPairs, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}
void FluidSortingGrid_OpenCL_Program::rearrangeParticleArrays(cl_command_queue commandQueue, int numFluidParticles, void *fluidBufferAddress)
{
	cl_int error_code;

	///__global ValueIndexPair *sortedPairs
	///__global btVector3 *rearrange
	///__global btVector3 *temporary
	error_code = clSetKernelArg( kernel_rearrangeParticleArrays, 0, sizeof(void*), buffer_pairs.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_rearrangeParticleArrays, 1, sizeof(void*), fluidBufferAddress );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_rearrangeParticleArrays, 2, sizeof(void*), buffer_temp.getAddress() );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(commandQueue, kernel_rearrangeParticleArrays, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}
void FluidSortingGrid_OpenCL_Program::rearrangeParticleArraysWriteBack(cl_command_queue commandQueue, int numFluidParticles, 
																	   void *fluidBufferAddress)
{
	cl_int error_code;
	
	///__global btVector3 *rearrange
	///__global btVector3 *temporary
	error_code = clSetKernelArg( kernel_rearrangeParticleArraysWriteBack, 0, sizeof(void*), fluidBufferAddress );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_rearrangeParticleArraysWriteBack, 1, sizeof(void*), buffer_temp.getAddress() );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(commandQueue, kernel_rearrangeParticleArraysWriteBack, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(commandQueue);
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
	error_code = clSetKernelArg( kernel_generateUniques, 0, sizeof(void*), buffer_pairs.getAddress() );
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
	
	error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}
*/
