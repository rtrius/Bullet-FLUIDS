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

#include "LinearMath/btQuickprof.h"		//BT_PROFILE(name) macro

#include "btExperimentsOpenCL/btLauncherCL.h"
#include "btExperimentsOpenCL/btOpenCLUtils.h"

#include "BulletFluids/Sph/btFluidSphParameters.h"

#include "fluidSphCL.h"

btFluidSphSolverOpenCL::btFluidSphSolverOpenCL(cl_context context, cl_command_queue queue, cl_device_id device)
: m_globalFluidParams(context, queue), m_sortingGridProgram(context, queue, device)
{
	m_context = context;
	m_commandQueue = queue;

	//
	const char CL_PROGRAM_PATH[] = "./Demos/FluidSphDemo/BulletMultiThreaded/SphSolverOpenCL/fluidSph.cl";
	
	const char* kernelSource = fluidSphCL;	//fluidSphCL.h
	cl_int error;
	char* additionalMacros = 0;
	m_fluidsProgram = btOpenCLUtils::compileCLProgramFromString(context, device, kernelSource, &error, 
																	additionalMacros, CL_PROGRAM_PATH);
	btAssert(m_fluidsProgram);
	
	m_kernel_findNeighborCellsPerCell = clCreateKernel(m_fluidsProgram, "findNeighborCellsPerCell", 0);
	btAssert(m_kernel_findNeighborCellsPerCell);
	m_kernel_findGridCellIndexPerParticle = clCreateKernel(m_fluidsProgram, "findGridCellIndexPerParticle", 0);
	btAssert(m_kernel_findGridCellIndexPerParticle);
	m_kernel_sphComputePressure = clCreateKernel(m_fluidsProgram, "sphComputePressure", 0);
	btAssert(m_kernel_sphComputePressure);
	m_kernel_sphComputeForce = clCreateKernel(m_fluidsProgram, "sphComputeForce", 0);
	btAssert(m_kernel_sphComputeForce);
}

btFluidSphSolverOpenCL::~btFluidSphSolverOpenCL()
{
	clReleaseKernel(m_kernel_findNeighborCellsPerCell);
	clReleaseKernel(m_kernel_findGridCellIndexPerParticle);
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
		if( fluids[i]->numParticles() ) 
		{
			validFluids.push_back( fluids[i] );
		}
		else
		{
			//Update AABB
			btFluidSortingGrid& grid = fluids[i]->internalGetGrid();
			btVector3& pointAabbMin = grid.internalGetPointAabbMin();
			btVector3& pointAabbMax = grid.internalGetPointAabbMax();
			
			pointAabbMin.setValue(0,0,0);
			pointAabbMax.setValue(0,0,0);
		}
	}

	int numValidFluids = validFluids.size();
	
//BT_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT is not supported when using OpenCL grid update.
#ifdef BT_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT
	const bool UPDATE_GRID_ON_GPU = false;
#else
	const bool UPDATE_GRID_ON_GPU = true;
#endif
	
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
			
			if(!UPDATE_GRID_ON_GPU) m_gridData[i]->writeToOpenCL( m_commandQueue, validFluids[i]->internalGetGrid() );
			m_fluidData[i]->writeToOpenCL( m_commandQueue, FL, validFluids[i]->internalGetParticles() );
		}
	}
	
	//
	{
		BT_PROFILE("calculate sph force");
	
		for(int i = 0; i < numValidFluids; ++i)
		{
			btFluidSph* fluid = validFluids[i];
			btFluidSortingGridOpenCL* gridData = m_gridData[i];
			btFluidSphOpenCL* fluidData = m_fluidData[i];
			
			if(UPDATE_GRID_ON_GPU)
				m_sortingGridProgram.insertParticlesIntoGrid(m_context, m_commandQueue, fluid, m_fluidData[i], m_gridData[i]);
			
			int numActiveCells = gridData->getNumActiveCells();
			int numFluidParticles = fluid->numParticles();
			
			findNeighborCells( numActiveCells, numFluidParticles, gridData, fluidData);
			sphComputePressure( numFluidParticles, gridData, fluidData, fluid->getGrid().getCellSize() );
			sphComputeForce( numFluidParticles, gridData, fluidData, fluid->getGrid().getCellSize() );
			
			clFlush(m_commandQueue);
			
			//This branch is executed on CPU, while findNeighborCells() and sphComputePressure/Force() is simultaneously executed on GPU
			if(UPDATE_GRID_ON_GPU)
			{
				BT_PROFILE("simultaneous rearrange/AABB update");
			
				m_sortingGridProgram.rearrangeParticlesOnHost(fluid);
				
				{
					BT_PROFILE("Update AABB - CPU");
					btVector3 aabbMin(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);
					btVector3 aabbMax(-BT_LARGE_FLOAT, -BT_LARGE_FLOAT, -BT_LARGE_FLOAT);
					for(int i = 0; i < fluid->numParticles(); ++i) 
					{
						const btVector3& position = fluid->getPosition(i);
							
						aabbMin.setMin(position);
						aabbMax.setMax(position);
					}
					
					btFluidSortingGrid& grid = fluid->internalGetGrid();
					btVector3& pointAabbMin = grid.internalGetPointAabbMin();
					btVector3& pointAabbMax = grid.internalGetPointAabbMax();
					
					pointAabbMin = aabbMin;
					pointAabbMax = aabbMax;
				}
			}
			
			clFinish(m_commandQueue);
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


void btFluidSphSolverOpenCL::findNeighborCells(int numActiveGridCells, int numFluidParticles, 
												btFluidSortingGridOpenCL* gridData, btFluidSphOpenCL* fluidData)
{
	//BT_PROFILE("findNeighborCells");
	
	//Perform 9 binary searches per cell, to locate the 27 neighbor cells
	{
		gridData->m_foundCells.resize(numActiveGridCells);
		
		btBufferInfoCL bufferInfo[] = 
		{ 
			btBufferInfoCL( gridData->m_numActiveCells.getBufferCL() ),
			btBufferInfoCL( gridData->m_activeCells.getBufferCL() ),
			btBufferInfoCL( gridData->m_cellContents.getBufferCL() ),
			btBufferInfoCL( gridData->m_foundCells.getBufferCL() )
		};
		
		btLauncherCL launcher(m_commandQueue, m_kernel_findNeighborCellsPerCell);
		launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
		
		launcher.launch1D(numActiveGridCells);
	}
	
	//For each particle, locate the grid cell that they are contained in so that 
	//they can use the results from m_kernel_findNeighborCellsPerCell, executed above
	{
		btBufferInfoCL bufferInfo[] = 
		{
			btBufferInfoCL( gridData->m_numActiveCells.getBufferCL() ),
			btBufferInfoCL( gridData->m_cellContents.getBufferCL() ),
			btBufferInfoCL( fluidData->m_cellIndex.getBufferCL() )
		};
		
		btLauncherCL launcher(m_commandQueue, m_kernel_findGridCellIndexPerParticle);
		launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
		
		launcher.launch1D(numActiveGridCells);
	}
}
void btFluidSphSolverOpenCL::sphComputePressure(int numFluidParticles, btFluidSortingGridOpenCL* gridData, btFluidSphOpenCL* fluidData, btScalar cellSize) 
{
	//BT_PROFILE("sphComputePressure");
	
	btBufferInfoCL bufferInfo[] = 
	{ 
		btBufferInfoCL( m_globalFluidParams.getBufferCL() ), 
		btBufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
		btBufferInfoCL( fluidData->m_pos.getBufferCL() ),
		btBufferInfoCL( fluidData->m_density.getBufferCL() ),
		btBufferInfoCL( gridData->m_foundCells.getBufferCL() ),
		btBufferInfoCL( fluidData->m_cellIndex.getBufferCL() )
	};
	
	btLauncherCL launcher(m_commandQueue, m_kernel_sphComputePressure);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
	//clFinish(m_commandQueue);
}
void btFluidSphSolverOpenCL::sphComputeForce(int numFluidParticles, btFluidSortingGridOpenCL* gridData, btFluidSphOpenCL* fluidData, btScalar cellSize) 
{
	//BT_PROFILE("sphComputeForce");
	
	btBufferInfoCL bufferInfo[] = 
	{ 
		btBufferInfoCL( m_globalFluidParams.getBufferCL() ), 
		btBufferInfoCL( fluidData->m_localParameters.getBufferCL() ),
		btBufferInfoCL( fluidData->m_pos.getBufferCL() ),
		btBufferInfoCL( fluidData->m_vel_eval.getBufferCL() ),
		btBufferInfoCL( fluidData->m_sph_force.getBufferCL() ),
		btBufferInfoCL( fluidData->m_density.getBufferCL() ),
		btBufferInfoCL( gridData->m_foundCells.getBufferCL() ),
		btBufferInfoCL( fluidData->m_cellIndex.getBufferCL() )
	};
	
	btLauncherCL launcher(m_commandQueue, m_kernel_sphComputeForce);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
	launcher.setConst(numFluidParticles);
	
	launcher.launch1D(numFluidParticles);
	//clFinish(m_commandQueue);
}
