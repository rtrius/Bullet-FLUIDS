/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#include "btFluidSphSolverOpenCL.h"

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

#include "btExperimentsOpenCL/btLauncherCL.h"

#include "../btFluidSphParameters.h"

btFluidSphSolverOpenCL::btFluidSphSolverOpenCL(cl_context context, cl_command_queue queue, cl_device_id device)
: m_globalFluidParams(context, queue), m_sortingGridProgram(context, queue, device)
{
	m_context = context;
	m_commandQueue = queue;

	//
	const char CL_PROGRAM_PATH[] = "./Demos/FluidSphDemo/BulletFluids/Sph/OpenCL_support/fluids.cl";
	m_fluidsProgram = compileProgramOpenCL(context, device, CL_PROGRAM_PATH);
	btAssert(m_fluidsProgram);
	
	m_kernel_sphComputePressure = clCreateKernel(m_fluidsProgram, "sphComputePressure", 0);
	btAssert(m_kernel_sphComputePressure);
	m_kernel_sphComputeForce = clCreateKernel(m_fluidsProgram, "sphComputeForce", 0);
	btAssert(m_kernel_sphComputeForce);
}

btFluidSphSolverOpenCL::~btFluidSphSolverOpenCL()
{
	clReleaseKernel(m_kernel_sphComputePressure);
	clReleaseKernel(m_kernel_sphComputeForce);
	clReleaseProgram(m_fluidsProgram);
}

void btFluidSphSolverOpenCL::updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids)
{	
	BT_PROFILE("btFluidSphSolverOpenCL::updateGridAndCalculateSphForces()");
	
#ifdef BT_USE_DOUBLE_PRECISION
	btAssert(0 && "BT_USE_DOUBLE_PRECISION not supported on OpenCL.\n");
	return;
#endif	

	btAlignedObjectArray<btFluidSph*> validFluids;
	for(int i = 0; i < numFluids; ++i) 
	{
		if( fluids[i]->numParticles() ) validFluids.push_back( fluids[i] );
	}

	int numValidFluids = validFluids.size();
	
	//#define BT_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT is not supported on OpenCL
	//If there are less than 32768 particles, CPU performance is equal to or faster than OpenCL.
	const bool UPDATE_GRID_ON_GPU = false;
	
	if(!UPDATE_GRID_ON_GPU)
		for(int i = 0; i < numValidFluids; ++i) validFluids[i]->insertParticlesIntoGrid();
	
	//Write data from CPU to OpenCL
	m_globalFluidParams.resize(1);
	m_globalFluidParams.copyFromHostPointer(&FG, 1, 0, false);
	clFinish(m_commandQueue);

	//resize m_gridData to numValidFluids
	{
		while(m_gridData.size() > numValidFluids)
		{
			btFluidSortingGridOpenCL* lastElement = m_gridData[m_gridData.size() - 1];
			lastElement->~btFluidSortingGridOpenCL();
			btAlignedFree(lastElement);
			
			m_gridData.pop_back();
		}
		while(m_gridData.size() < numValidFluids)
		{
			void* ptr = btAlignedAlloc( sizeof(btFluidSortingGridOpenCL), 16 );
			btFluidSortingGridOpenCL* newElement = new(ptr) btFluidSortingGridOpenCL(m_context, m_commandQueue);
			
			m_gridData.push_back(newElement);
		}
	}
	//resize m_fluidData to numValidFluids
	{
		while(m_fluidData.size() > numValidFluids)
		{
			btFluidSphOpenCL* lastElement = m_fluidData[m_fluidData.size() - 1];
			lastElement->~btFluidSphOpenCL();
			btAlignedFree(lastElement);
			
			m_fluidData.pop_back();
		}
		while(m_fluidData.size() < numValidFluids)
		{
			void* ptr = btAlignedAlloc( sizeof(btFluidSphOpenCL), 16 );
			btFluidSphOpenCL* newElement = new(ptr) btFluidSphOpenCL(m_context, m_commandQueue);
			
			m_fluidData.push_back(newElement);
		}
	}
	
	{
		BT_PROFILE("writeToOpenCL");
		for(int i = 0; i < numValidFluids; ++i)
		{
			const btFluidSphParametersLocal& FL = validFluids[i]->getLocalParameters();
			
			m_gridData[i]->writeToOpenCL( m_commandQueue, validFluids[i]->internalGetGrid() );
			m_fluidData[i]->writeToOpenCL( m_commandQueue, FL, validFluids[i]->internalGetParticles() );
		}
	}

	//
	if(UPDATE_GRID_ON_GPU)
	{
		BT_PROFILE("update grid");
	
		for(int i = 0; i < numValidFluids; ++i)
			m_sortingGridProgram.insertParticlesIntoGrid(m_context, m_commandQueue, validFluids[i], m_fluidData[i], m_gridData[i]);
	}
	
	//
	{
		BT_PROFILE("calculate sph force");
	
		for(int i = 0; i < numValidFluids; ++i)
		{
			int numFluidParticles = validFluids[i]->numParticles();
			
			btFluidSortingGridOpenCL* gridData = m_gridData[i];
			btFluidSphOpenCL* fluidData = m_fluidData[i];
			
			sphComputePressure( numFluidParticles, gridData, fluidData, validFluids[i]->getGrid().getCellSize() );
			sphComputeForce( numFluidParticles, gridData, fluidData, validFluids[i]->getGrid().getCellSize() );
		}
	}
	
	//Read data from OpenCL to CPU
	{
		BT_PROFILE("readFromOpenCL, applySphForce()");
		for(int i = 0; i < numValidFluids; ++i)
		{
			if( m_tempSphForce.size() < validFluids[i]->numParticles() ) m_tempSphForce.resize( validFluids[i]->numParticles() );
		
			if(UPDATE_GRID_ON_GPU) m_gridData[i]->readFromOpenCL( m_commandQueue, validFluids[i]->internalGetGrid() );
			m_fluidData[i]->readFromOpenCL( m_commandQueue, m_tempSphForce );
			
			applySphForce(FG, validFluids[i], m_tempSphForce);
		}
	}
}

void btFluidSphSolverOpenCL::sphComputePressure(int numFluidParticles, btFluidSortingGridOpenCL* gridData, btFluidSphOpenCL* fluidData, btScalar cellSize) 
{
	BT_PROFILE("btFluidSphSolverOpenCL::sphComputePressure()");
	
	btBufferInfoCL bufferInfo[] = 
	{ 
		btBufferInfoCL( m_globalFluidParams.getBufferCL() ), 
		btBufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
		btBufferInfoCL( fluidData->m_pos.getBufferCL() ),
		btBufferInfoCL( fluidData->m_density.getBufferCL() ),
		btBufferInfoCL( gridData->m_numActiveCells.getBufferCL() ),
		btBufferInfoCL( gridData->m_activeCells.getBufferCL() ),
		btBufferInfoCL( gridData->m_cellContents.getBufferCL() )
	};
	
	btLauncherCL launcher(m_commandQueue, m_kernel_sphComputePressure);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
	launcher.setConst(cellSize);
	
	launcher.launchAutoSizedWorkGroups1D(numFluidParticles);
	
	clFinish(m_commandQueue);
}
void btFluidSphSolverOpenCL::sphComputeForce(int numFluidParticles, btFluidSortingGridOpenCL* gridData, btFluidSphOpenCL* fluidData, btScalar cellSize) 
{
	BT_PROFILE("btFluidSphSolverOpenCL::sphComputeForce()");
	
	btBufferInfoCL bufferInfo[] = 
	{ 
		btBufferInfoCL( m_globalFluidParams.getBufferCL() ), 
		btBufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
		btBufferInfoCL( fluidData->m_pos.getBufferCL() ),
		btBufferInfoCL( fluidData->m_vel_eval.getBufferCL() ),
		btBufferInfoCL( fluidData->m_sph_force.getBufferCL() ),
		btBufferInfoCL( fluidData->m_density.getBufferCL() ),
		btBufferInfoCL( gridData->m_numActiveCells.getBufferCL() ),
		btBufferInfoCL( gridData->m_activeCells.getBufferCL() ),
		btBufferInfoCL( gridData->m_cellContents.getBufferCL() )
	};
	
	btLauncherCL launcher(m_commandQueue, m_kernel_sphComputeForce);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
	launcher.setConst(cellSize);
	
	launcher.launchAutoSizedWorkGroups1D(numFluidParticles);
	
	clFinish(m_commandQueue);
}
