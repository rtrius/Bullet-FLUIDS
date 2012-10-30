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
#include "btFluidSortingGridOpenCL.h"

#include "btExperimentsOpenCL/btLauncherCL.h"

#include "../btFluidSortingGrid.h"

#include "btFluidSphOpenCL.h"

// /////////////////////////////////////////////////////////////////////////////
//class btFluidSortingGridOpenCL
// /////////////////////////////////////////////////////////////////////////////
void btFluidSortingGridOpenCL::writeToOpenCL(cl_command_queue queue, btFluidSortingGrid* sortingGrid)
{
	int numActiveCells = sortingGrid->internalGetActiveCells().size();

	m_numActiveCells.resize(1);
	m_numActiveCells.copyFromHostPointer(&numActiveCells, 1, 0, false);
	
	m_activeCells.copyFromHost( sortingGrid->internalGetActiveCells(), false );
	m_cellContents.copyFromHost( sortingGrid->internalGetCellContents(), false );
	
	clFinish(queue);
}
void btFluidSortingGridOpenCL::readFromOpenCL(cl_command_queue queue, btFluidSortingGrid* sortingGrid)
{
	m_activeCells.copyToHost( sortingGrid->internalGetActiveCells(), false );
	m_cellContents.copyToHost( sortingGrid->internalGetCellContents(), false );
	
	clFinish(queue);
}

int btFluidSortingGridOpenCL::getNumActiveCells() const
{
	int numActiveCells;
	m_numActiveCells.copyToHostPointer(&numActiveCells, 1, 0, true);
	
	return numActiveCells;
}

// /////////////////////////////////////////////////////////////////////////////
//class btFluidSortingGridOpenCLProgram
// /////////////////////////////////////////////////////////////////////////////
btFluidSortingGridOpenCLProgram::btFluidSortingGridOpenCLProgram(cl_context context, cl_command_queue queue, cl_device_id device)
: m_tempBuffer(context, queue), m_radixSorter(context, device, queue), m_valueIndexPairs(context, queue)
{
	const char CL_SORTING_GRID_PROGRAM_PATH[] = "./Demos/FluidDemo/Fluids/OpenCL_support/fluids.cl";
	sortingGrid_program = compileProgramOpenCL(context, device, CL_SORTING_GRID_PROGRAM_PATH);
	btAssert(sortingGrid_program);

	kernel_generateValueIndexPairs = clCreateKernel(sortingGrid_program, "generateValueIndexPairs", 0);
	btAssert(kernel_generateValueIndexPairs);
	kernel_rearrangeParticleArrays = clCreateKernel(sortingGrid_program, "rearrangeParticleArrays", 0);
	btAssert(kernel_rearrangeParticleArrays);
	kernel_generateUniques = clCreateKernel(sortingGrid_program, "generateUniques", 0);
	btAssert(kernel_generateUniques);
}
btFluidSortingGridOpenCLProgram::~btFluidSortingGridOpenCLProgram()
{
	clReleaseKernel(kernel_generateValueIndexPairs);
	clReleaseKernel(kernel_rearrangeParticleArrays);
	clReleaseKernel(kernel_generateUniques);
	clReleaseProgram(sortingGrid_program);
}


void rearrangeToMatchSortedValues2(const btAlignedObjectArray<btSortData>& sortedValues, btAlignedObjectArray<btVector3>& out_rearranged)
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
void btFluidSortingGridOpenCLProgram::insertParticlesIntoGrid(cl_context context, cl_command_queue commandQueue,
															  btFluidSph* fluid, btFluidSphOpenCL* fluidData, btFluidSortingGridOpenCL* gridData)
{
	BT_PROFILE("btFluidSortingGridOpenCLProgram::insertParticlesIntoGrid()");
	
	int numFluidParticles = fluid->numParticles();
	m_tempBuffer.resize(numFluidParticles);
	m_valueIndexPairs.resize(numFluidParticles);
	
	//Cannot check number of nonempty grid cells before generateUniques();
	//temporarily resize m_activeCells and m_cellContents
	//to handle the case where each particle occupies a different grid cell.
	gridData->m_numActiveCells.resize(1);
	gridData->m_activeCells.resize(numFluidParticles);
	gridData->m_cellContents.resize(numFluidParticles);
	
	//
	{
		BT_PROFILE("generateValueIndexPairs()");
		generateValueIndexPairs( commandQueue, numFluidParticles, fluid->getGrid().getCellSize(), fluidData->m_pos.getBufferCL() );
		
		clFinish(commandQueue);
	}
	
	//Note that btRadixSort32CL uses btSortData, while btFluidSortingGrid uses btValueIndexPair.
	//btSortData.m_key == btValueIndexPair.m_value (value to sort by)
	//btSortData.m_value == btValueIndexPair.m_index (fluid particle index)
	{
		BT_PROFILE("radix sort");
		m_radixSorter.execute(m_valueIndexPairs, 32);
	}
	
	//
	{
		BT_PROFILE("rearrange device");
		
		rearrangeParticleArrays( commandQueue, numFluidParticles, fluidData->m_pos.getBufferCL() );
		fluidData->m_pos.copyFromOpenCLArray(m_tempBuffer);
		
		rearrangeParticleArrays( commandQueue, numFluidParticles, fluidData->m_vel_eval.getBufferCL() );
		fluidData->m_vel_eval.copyFromOpenCLArray(m_tempBuffer);
		
		clFinish(commandQueue);
	}
	
	{
		BT_PROFILE("rearrange host");
		m_valueIndexPairs.copyToHost(m_valueIndexPairsHost, true);
		
		btFluidParticles& particles = fluid->internalGetParticles();
		rearrangeToMatchSortedValues2(m_valueIndexPairsHost, particles.m_pos);
		rearrangeToMatchSortedValues2(m_valueIndexPairsHost, particles.m_vel);
		rearrangeToMatchSortedValues2(m_valueIndexPairsHost, particles.m_vel_eval);
		rearrangeToMatchSortedValues2(m_valueIndexPairsHost, particles.m_externalAcceleration);
	}
	
	//
	{
		BT_PROFILE("generateUniques_serial()");
		generateUniques(commandQueue, numFluidParticles, gridData);
		
		clFinish(commandQueue);
	}
	
	int numActiveCells = gridData->getNumActiveCells();
	gridData->m_activeCells.resize(numActiveCells);
	gridData->m_cellContents.resize(numActiveCells);
}

void btFluidSortingGridOpenCLProgram::generateValueIndexPairs(cl_command_queue commandQueue, int numFluidParticles, 
															  btScalar cellSize, cl_mem fluidPositionsBuffer)
{
	btBufferInfoCL bufferInfo[] = 
	{
		btBufferInfoCL( fluidPositionsBuffer ),
		btBufferInfoCL( m_valueIndexPairs.getBufferCL() )
	};
	
	btLauncherCL launcher(commandQueue, kernel_generateValueIndexPairs);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
	launcher.setConst(cellSize);
	
	launcher.launchAutoSizedWorkGroups1D(numFluidParticles);
}
void btFluidSortingGridOpenCLProgram::rearrangeParticleArrays(cl_command_queue commandQueue, int numFluidParticles, cl_mem fluidBuffer)
{
	btBufferInfoCL bufferInfo[] = 
	{
		btBufferInfoCL( m_valueIndexPairs.getBufferCL() ),
		btBufferInfoCL( fluidBuffer ),
		btBufferInfoCL( m_tempBuffer.getBufferCL() )
	};
	
	btLauncherCL launcher(commandQueue, kernel_rearrangeParticleArrays);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
	
	launcher.launchAutoSizedWorkGroups1D(numFluidParticles);
}

void btFluidSortingGridOpenCLProgram::generateUniques(cl_command_queue commandQueue, int numFluidParticles, btFluidSortingGridOpenCL* gridData)
{
	btBufferInfoCL bufferInfo[] = 
	{
		btBufferInfoCL( m_valueIndexPairs.getBufferCL() ),
		btBufferInfoCL( gridData->m_activeCells.getBufferCL() ),
		btBufferInfoCL( gridData->m_cellContents.getBufferCL() ),
		btBufferInfoCL( gridData->m_numActiveCells.getBufferCL() )
	};
	
	btLauncherCL launcher(commandQueue, kernel_generateUniques);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
	launcher.setConst(numFluidParticles);
	
	launcher.launchAutoSizedWorkGroups1D(1);
}
