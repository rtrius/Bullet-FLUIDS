/* FluidSolverOpenCL.cpp
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
#include "FluidSolverOpenCL.h"

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

#include "btExperimentsOpenCL/btLauncherCL.h"

#include "../FluidParameters.h"

FluidSolverOpenCL::FluidSolverOpenCL(cl_context context, cl_command_queue queue, cl_device_id device)
: m_globalFluidParams(context, queue), m_sortingGridProgram(context, queue, device)
{
	m_context = context;
	m_commandQueue = queue;

	//
	const char CL_PROGRAM_PATH[] = "./Demos/FluidDemo/Fluids/OpenCL_support/fluids.cl";
	m_fluidsProgram = compileProgramOpenCL(context, device, CL_PROGRAM_PATH);
	btAssert(m_fluidsProgram);
	
	m_kernel_sphComputePressure = clCreateKernel(m_fluidsProgram, "sphComputePressure", 0);
	btAssert(m_kernel_sphComputePressure);
	m_kernel_sphComputeForce = clCreateKernel(m_fluidsProgram, "sphComputeForce", 0);
	btAssert(m_kernel_sphComputeForce);
}

FluidSolverOpenCL::~FluidSolverOpenCL()
{
	clReleaseKernel(m_kernel_sphComputePressure);
	clReleaseKernel(m_kernel_sphComputeForce);
	clReleaseProgram(m_fluidsProgram);
}

void FluidSolverOpenCL::stepSimulation(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids)
{	
	BT_PROFILE("FluidSolverOpenCL::stepSimulation()");
	
#ifdef BT_USE_DOUBLE_PRECISION
	printf("BT_USE_DOUBLE_PRECISION not supported on OpenCL.\n");
	btAssert(0);
	return;
#endif	

	btAlignedObjectArray<FluidSph*> validFluids;
	for(int i = 0; i < fluids->size(); ++i) 
	{
		if( (*fluids)[i]->numParticles() ) validFluids.push_back( (*fluids)[i] );
	}

	int numValidFluids = validFluids.size();
	
	for(int i = 0; i < numValidFluids; ++i) validFluids[i]->insertParticlesIntoGrid();
	
	//Write data from CPU to OpenCL
	m_globalFluidParams.resize(1);
	m_globalFluidParams.copyFromHostPointer(&FG, 1, 0, false);
	clFinish(m_commandQueue);

	//resize m_gridData to numValidFluids
	{
		while(m_gridData.size() > numValidFluids)
		{
			FluidSortingGridOpenCL *lastElement = m_gridData[m_gridData.size() - 1];
			lastElement->~FluidSortingGridOpenCL();
			btAlignedFree(lastElement);
			
			m_gridData.pop_back();
		}
		while(m_gridData.size() < numValidFluids)
		{
			void *ptr = btAlignedAlloc( sizeof(FluidSortingGridOpenCL), 16 );
			FluidSortingGridOpenCL *newElement = new(ptr) FluidSortingGridOpenCL(m_context, m_commandQueue);
			
			m_gridData.push_back(newElement);
		}
	}
	//resize m_fluidData to numValidFluids
	{
		while(m_fluidData.size() > numValidFluids)
		{
			FluidSphOpenCL *lastElement = m_fluidData[m_fluidData.size() - 1];
			lastElement->~FluidSphOpenCL();
			btAlignedFree(lastElement);
			
			m_fluidData.pop_back();
		}
		while(m_fluidData.size() < numValidFluids)
		{
			void *ptr = btAlignedAlloc( sizeof(FluidSphOpenCL), 16 );
			FluidSphOpenCL *newElement = new(ptr) FluidSphOpenCL(m_context, m_commandQueue);
			
			m_fluidData.push_back(newElement);
		}
	}
	
	{
		BT_PROFILE("stepSimulation() - writeToOpenCL");
		for(int i = 0; i < numValidFluids; ++i)
		{
			const FluidParametersLocal &FL = validFluids[i]->getLocalParameters();
			
			m_gridData[i]->writeToOpenCL( m_commandQueue, &validFluids[i]->internalGetGrid() );
			m_fluidData[i]->writeToOpenCL( m_commandQueue, FL, &validFluids[i]->internalGetFluidParticles() );
		}
	}

	//
	{
		BT_PROFILE("stepSimulation() - grid update, sph force");
	
		for(int i = 0; i < numValidFluids; ++i)
		{
			int numFluidParticles = validFluids[i]->numParticles();
			
			FluidSortingGridOpenCL *gridData = m_gridData[i];
			FluidSphOpenCL *fluidData = m_fluidData[i];
			
			//m_sortingGridProgram.insertParticlesIntoGrid(m_context, m_commandQueue, validFluids[i], fluidData, gridData);
			
			sphComputePressure( numFluidParticles, gridData, fluidData, validFluids[i]->getGrid().getCellSize() );
			sphComputeForce(numFluidParticles, gridData, fluidData);
		}
	}
	
	//Read data from OpenCL to CPU
	{
		BT_PROFILE("stepSimulation() - readFromOpenCL");
		for(int i = 0; i < numValidFluids; ++i)
		{
			m_gridData[i]->readFromOpenCL( m_commandQueue, &validFluids[i]->internalGetGrid() );
			m_fluidData[i]->readFromOpenCL( m_commandQueue, &validFluids[i]->internalGetFluidParticles() );
		}
	}
	
	//
	for(int i = 0; i < numValidFluids; ++i)
	{
		integrate( FG, validFluids[i]->getLocalParameters(), &validFluids[i]->internalGetFluidParticles() );
	}
}

void FluidSolverOpenCL::sphComputePressure(int numFluidParticles, FluidSortingGridOpenCL *gridData,
											FluidSphOpenCL *fluidData, btScalar cellSize) 
{
	BT_PROFILE("FluidSolverOpenCL::sphComputePressure()");
	
	btBufferInfoCL bufferInfo[] = 
	{ 
		btBufferInfoCL( m_globalFluidParams.getBufferCL() ), 
		btBufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
		btBufferInfoCL( fluidData->m_pos.getBufferCL() ),
		btBufferInfoCL( fluidData->m_pressure.getBufferCL() ),
		btBufferInfoCL( fluidData->m_invDensity.getBufferCL() ),
		btBufferInfoCL( fluidData->m_neighborTable.getBufferCL() ),
		btBufferInfoCL( gridData->m_numActiveCells.getBufferCL() ),
		btBufferInfoCL( gridData->m_activeCells.getBufferCL() ),
		btBufferInfoCL( gridData->m_cellContents.getBufferCL() )	};
	
	btLauncherCL launcher(m_commandQueue, m_kernel_sphComputePressure);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
	launcher.setConst(cellSize);
	
	launcher.launchAutoSizedWorkGroups1D(numFluidParticles);
	
	clFinish(m_commandQueue);
}
void FluidSolverOpenCL::sphComputeForce(int numFluidParticles, FluidSortingGridOpenCL *gridData, FluidSphOpenCL *fluidData) 
{
	BT_PROFILE("FluidSolverOpenCL::sphComputeForce()");
	
	btBufferInfoCL bufferInfo[] = 
	{ 
		btBufferInfoCL( m_globalFluidParams.getBufferCL() ), 
		btBufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
		btBufferInfoCL( fluidData->m_pos.getBufferCL() ),
		btBufferInfoCL( fluidData->m_vel_eval.getBufferCL() ),
		btBufferInfoCL( fluidData->m_sph_force.getBufferCL() ),
		btBufferInfoCL( fluidData->m_pressure.getBufferCL() ),
		btBufferInfoCL( fluidData->m_invDensity.getBufferCL() ),
		btBufferInfoCL( fluidData->m_neighborTable.getBufferCL() )
	};
	
	btLauncherCL launcher(m_commandQueue, m_kernel_sphComputeForce);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
	
	launcher.launchAutoSizedWorkGroups1D(numFluidParticles);
	
	clFinish(m_commandQueue);
}
