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

FluidSolverOpenCL::FluidSolverOpenCL()
{
	m_fluidsProgram = 0;
	m_kernel_sphComputePressure = 0;
	m_kernel_sphComputeForce = 0;
	
	//
	initialize();
}

void FluidSolverOpenCL::initialize()
{
	m_configCL.initialize();
	
	//
	cl_int error_code;
	
	//Load and build program
	if(m_configCL.m_context && m_configCL.m_commandQueue)
	{
		const char CL_PROGRAM_PATH[] = "./Demos/FluidDemo/Fluids/OpenCL_support/fluids.cl";
		m_fluidsProgram = compileProgramOpenCL(m_configCL.m_context, m_configCL.m_device, CL_PROGRAM_PATH);
		
		//Kernels
		m_kernel_sphComputePressure = clCreateKernel(m_fluidsProgram, "sphComputePressure", &error_code);
		CHECK_CL_ERROR(error_code);
		m_kernel_sphComputeForce = clCreateKernel(m_fluidsProgram, "sphComputeForce", &error_code);
		CHECK_CL_ERROR(error_code);
		
		//Buffers
		m_buffer_globalFluidParams.allocate( m_configCL.m_context, sizeof(FluidParametersGlobal) );
	}
	else printf("FluidSolverOpenCL::initialize_stage4_program_and_buffers() error: invalid m_configCL.m_context or command_queue. \n");
	
	//
	m_sortingGridProgram.initialize(m_configCL.m_context, m_configCL.m_device, m_configCL.m_commandQueue);
}

void FluidSolverOpenCL::deactivate()
{
	m_sortingGridProgram.deactivate();

	//
	cl_int error_code;
	
	//Buffers
	m_buffer_globalFluidParams.deallocate();
	
	//Kernels
	error_code = clReleaseKernel(m_kernel_sphComputePressure);
	CHECK_CL_ERROR(error_code);
	error_code = clReleaseKernel(m_kernel_sphComputeForce);
	CHECK_CL_ERROR(error_code);
		
	//Program
	error_code = clReleaseProgram(m_fluidsProgram);
	CHECK_CL_ERROR(error_code);
	
	//
	m_fluidsProgram = 0;
	m_kernel_sphComputePressure = 0;
	m_kernel_sphComputeForce = 0;
	
	//
	m_configCL.deactivate();
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
	
	//for(int i = 0; i < numValidFluids; ++i) validFluids[i]->insertParticlesIntoGrid();
	
	//Write data from CPU to OpenCL
	FluidParametersGlobal globalParameters = FG;
	m_buffer_globalFluidParams.writeToBuffer( m_configCL.m_commandQueue, &globalParameters, sizeof(FluidParametersGlobal) );
	cl_int error_code = clFinish(m_configCL.m_commandQueue);
	CHECK_CL_ERROR(error_code);

		//Resize m_gridData and m_fluidData to match numValidFluids:
		//Calling btAlignedObjectArray<T>::resize(n) with n > btAlignedObjectArray<T>::size()
		//may call the destructor, ~T(), for all existing objects in the array. As a result,
		//calling resize(numValidFluids) without resize(0) here may cause the OpenCL buffers
		//to be deallocated without setting them to 0, eventually causing
		//a segmentation fault.
	if( m_gridData.size() != numValidFluids )
	{
		m_gridData.resize(0);
		m_gridData.resize(numValidFluids);
	}	
	if( m_fluidData.size() != numValidFluids )
	{
		m_fluidData.resize(0);
		m_fluidData.resize(numValidFluids);
	}
	
	{
		BT_PROFILE("stepSimulation() - writeToOpenCL");
		for(int i = 0; i < numValidFluids; ++i)
		{
			const FluidParametersLocal &FL = validFluids[i]->getLocalParameters();
			
			m_gridData[i].writeToOpenCL( m_configCL.m_context, m_configCL.m_commandQueue, &validFluids[i]->internalGetGrid(), false );
			m_fluidData[i].writeToOpenCL( m_configCL.m_context, m_configCL.m_commandQueue, FL, &validFluids[i]->internalGetFluidParticles() );
		}
	}

	//
	{
		BT_PROFILE("stepSimulation() - grid update, sph force");
	
		for(int i = 0; i < numValidFluids; ++i)
		{
			int numFluidParticles = validFluids[i]->numParticles();
			
			//m_sortingGridProgram.insertParticlesIntoGrid() does not generate 
			//'cell processing groups', which are needed for symmetry optimizations.
			m_sortingGridProgram.insertParticlesIntoGrid(m_configCL.m_context, m_configCL.m_commandQueue, validFluids[i], &m_fluidData[i], &m_gridData[i]);
			
			FluidSortingGrid_OpenCL &gridData = m_gridData[i];
			Fluid_OpenCL &fluidData = m_fluidData[i];
			
			sphComputePressure( numFluidParticles, &gridData, &fluidData, validFluids[i]->getGrid().getCellSize() );
			sphComputeForce(numFluidParticles, &gridData, &fluidData);
		}
	}
	
	//Read data from OpenCL to CPU
	{
		BT_PROFILE("stepSimulation() - readFromOpenCL");
		for(int i = 0; i < numValidFluids; ++i)
		{
			m_gridData[i].readFromOpenCL( m_configCL.m_context, m_configCL.m_commandQueue, &validFluids[i]->internalGetGrid(), false );
			m_fluidData[i].readFromOpenCL( m_configCL.m_context, m_configCL.m_commandQueue, &validFluids[i]->internalGetFluidParticles() );
		}
	}
	
	//
	for(int i = 0; i < numValidFluids; ++i)
	{
		integrate( FG, validFluids[i]->getLocalParameters(), &validFluids[i]->internalGetFluidParticles() );
	}
}

void FluidSolverOpenCL::sphComputePressure(int numFluidParticles, FluidSortingGrid_OpenCL *gridData,
											Fluid_OpenCL *fluidData, btScalar cellSize) 
{
	BT_PROFILE("FluidSolverOpenCL::sphComputePressure()");
	
	btBufferInfoCL bufferInfo[] = 
	{ 
		btBufferInfoCL( m_buffer_globalFluidParams.getBuffer() ), 
		btBufferInfoCL( fluidData->m_buffer_localParameters.getBuffer() ),
		btBufferInfoCL( fluidData->m_buffer_pos.getBuffer() ),
		btBufferInfoCL( fluidData->m_buffer_pressure.getBuffer() ),
		btBufferInfoCL( fluidData->m_buffer_invDensity.getBuffer() ),
		btBufferInfoCL( fluidData->m_buffer_neighborTable.getBuffer() ),
		btBufferInfoCL( gridData->m_buffer_numActiveCells.getBuffer() ),
		btBufferInfoCL( gridData->m_buffer_activeCells.getBuffer() ),
		btBufferInfoCL( gridData->m_buffer_cellContents.getBuffer() )	};
	
	btLauncherCL launcher(m_configCL.m_commandQueue, m_kernel_sphComputePressure);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
	launcher.setConst(cellSize);
	
	launcher.launchAutoSizedWorkGroups1D(numFluidParticles);
	
	clFinish(m_configCL.m_commandQueue);
}
void FluidSolverOpenCL::sphComputeForce(int numFluidParticles, FluidSortingGrid_OpenCL *gridData, Fluid_OpenCL *fluidData) 
{
	BT_PROFILE("FluidSolverOpenCL::sphComputeForce()");
	
	btBufferInfoCL bufferInfo[] = 
	{ 
		btBufferInfoCL( m_buffer_globalFluidParams.getBuffer() ), 
		btBufferInfoCL( fluidData->m_buffer_localParameters.getBuffer() ),
		btBufferInfoCL( fluidData->m_buffer_pos.getBuffer() ),
		btBufferInfoCL( fluidData->m_buffer_vel_eval.getBuffer() ),
		btBufferInfoCL( fluidData->m_buffer_sph_force.getBuffer() ),
		btBufferInfoCL( fluidData->m_buffer_pressure.getBuffer() ),
		btBufferInfoCL( fluidData->m_buffer_invDensity.getBuffer() ),
		btBufferInfoCL( fluidData->m_buffer_neighborTable.getBuffer() )
	};
	
	btLauncherCL launcher(m_configCL.m_commandQueue, m_kernel_sphComputeForce);
	launcher.setBuffers( bufferInfo, sizeof(bufferInfo)/sizeof(btBufferInfoCL) );
	
	launcher.launchAutoSizedWorkGroups1D(numFluidParticles);
	
	clFinish(m_configCL.m_commandQueue);
}
