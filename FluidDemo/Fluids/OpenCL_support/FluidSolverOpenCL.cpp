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

#include "../FluidStaticGrid.h"
#include "../FluidParameters.h"

FluidSolverOpenCL::FluidSolverOpenCL()
{
	m_platformId = INVALID_PLATFORM_ID;
	m_context = INVALID_CONTEXT;
	
	m_device = INVALID_DEVICE_ID;
	m_commandQueue = INVALID_COMMAND_QUEUE;
	
	m_fluidsProgram = INVALID_PROGRAM;
	m_kernel_sphComputePressure = INVALID_KERNEL;
	m_kernel_sphComputeForce = INVALID_KERNEL;
	
	//
	initialize();
}

void FluidSolverOpenCL::initialize()
{
	initialize_stage1_platform();
	initialize_stage2_device();
	initialize_stage3_context_and_queue();
	initialize_stage4_program_and_buffer();
}

void FluidSolverOpenCL::deactivate()
{
	deactivate_stage1_program_and_buffer();
	deactivate_stage2_context_and_queue();
	
	//
	m_platformId = INVALID_PLATFORM_ID;
	m_context = INVALID_CONTEXT;
}

void FluidSolverOpenCL::stepSimulation(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids)
{	
	BT_PROFILE("FluidSolverOpenCL::stepSimulation()");
	
#ifdef BT_USE_DOUBLE_PRECISION
	printf("BT_USE_DOUBLE_PRECISION not supported on OpenCL.\n");
	btAssert(0);
	return;
#endif	

	//FluidSolverOpenCL requires use of 'class FluidSortingGrid'
	btAlignedObjectArray<FluidSph*> validFluids;
	for(int i = 0; i < fluids->size(); ++i) 
	{
		//GPU driver will crash if grid is incorrect type
		if( !(*fluids)[i]->numParticles() || (*fluids)[i]->getGrid()->getGridType() != FluidGrid::FT_IndexRange ) continue;
		
		validFluids.push_back( (*fluids)[i] );
	}

	int numValidFluids = validFluids.size();
	
	//	FluidSortingGrid::insertParticles() OpenCL program not implemented
	for(int i = 0; i < numValidFluids; ++i) validFluids[i]->insertParticlesIntoGrid();
	
	//Write data from CPU to OpenCL
	FluidParametersGlobal globalParameters = FG;
	m_buffer_globalFluidParams.writeToBuffer( m_commandQueue, &globalParameters, sizeof(FluidParametersGlobal) );
	cl_int error_code = clFinish(m_commandQueue);
	CHECK_CL_ERROR(error_code);
	
		//Resize m_gridData and m_fluidData to match numValidFluids:
		//Calling btAlignedObjectArray<T>::resize(n) with n > btAlignedObjectArray<T>::size()
		//may call the destructor, ~T(), for all existing objects in the array. As a result,
		//calling resize(numValidFluids) without resize(0) here may cause the OpenCL buffers
		//to be deallocated without setting them to INVALID_BUFFER, eventually causing
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
		
			m_gridData[i].writeToOpenCL( m_context, m_commandQueue, reinterpret_cast<FluidSortingGrid*>( validFluids[i]->internalGetGrid() ) );
			m_fluidData[i].writeToOpenCL( m_context, m_commandQueue, FL, &validFluids[i]->internalGetFluidParticles() );
		}
	}
	
	//
	{
		BT_PROFILE("stepSimulation() - grid update, sph force");
	
		for(int i = 0; i < numValidFluids; ++i)
		{
			int numFluidParticles = validFluids[i]->numParticles();
			FluidSortingGrid_OpenCLPointers gridPointers = m_gridData[i].getPointers();
			Fluid_OpenCLPointers fluidPointers = m_fluidData[i].getPointers();
			
			//grid_insertParticles(numFluidParticles, &gridPointers, &fluidPointers);
			sphComputePressure( numFluidParticles, &gridPointers, &fluidPointers, validFluids[i]->getGrid()->getCellSize() );
			sphComputeForce(numFluidParticles, &gridPointers, &fluidPointers);
		}
	}
	
	//Read data from OpenCL to CPU
	{
		BT_PROFILE("stepSimulation() - readFromOpenCL");
		for(int i = 0; i < numValidFluids; ++i)
		{
			//	FluidSortingGrid::insertParticles() OpenCL program not implemented
			//m_gridData[i].readFromOpenCL( m_context, m_commandQueue, reinterpret_cast<FluidSortingGrid*>( validFluids[i]->internalGetGrid() ) );
			m_fluidData[i].readFromOpenCL( m_context, m_commandQueue, &validFluids[i]->internalGetFluidParticles() );
		}
	}
	
	//
	for(int i = 0; i < numValidFluids; ++i)
	{
		integrate( FG, validFluids[i]->getLocalParameters(), &validFluids[i]->internalGetFluidParticles() );
	}
}

void FluidSolverOpenCL::initialize_stage1_platform()
{	
	cl_int error_code;
	const size_t MAX_STRING_LENGTH = 1024;
	char string[MAX_STRING_LENGTH];
	
	//Select platform
	cl_uint num_platforms;
	cl_platform_id platforms[MAX_PLATFORMS];
	error_code = clGetPlatformIDs(MAX_PLATFORMS, platforms, &num_platforms);
	CHECK_CL_ERROR(error_code);
	
	printf("Platforms available: \n");
	for(cl_uint i = 0; i < num_platforms; ++i)
	{
		error_code = clGetPlatformInfo(platforms[i], CL_PLATFORM_NAME, MAX_STRING_LENGTH, string, NULL);
		CHECK_CL_ERROR(error_code);
		printf("CL_PLATFORM_NAME: %s\n", string);
		
		error_code = clGetPlatformInfo(platforms[i], CL_PLATFORM_VENDOR, MAX_STRING_LENGTH, string, NULL);
		CHECK_CL_ERROR(error_code);
		printf("CL_PLATFORM_VENDOR: %s\n", string);
		
		error_code = clGetPlatformInfo(platforms[i], CL_PLATFORM_VERSION, MAX_STRING_LENGTH, string, NULL);
		CHECK_CL_ERROR(error_code);
		printf("CL_PLATFORM_VERSION: %s\n", string);
		
		error_code = clGetPlatformInfo(platforms[i], CL_PLATFORM_PROFILE, MAX_STRING_LENGTH, string, NULL);
		CHECK_CL_ERROR(error_code);
		printf("CL_PLATFORM_PROFILE: %s\n", string);
		
		error_code = clGetPlatformInfo(platforms[i], CL_PLATFORM_EXTENSIONS, MAX_STRING_LENGTH, string, NULL);
		CHECK_CL_ERROR(error_code);
		printf("CL_PLATFORM_EXTENSIONS: %s\n", string);
		
		//Select any platform with a GPU device
		cl_uint num_gpu_devices;
		cl_device_id devices[MAX_DEVICES];	//Replacing 'devices' with 'NULL' in clGetDeviceIDs() results in (error_code != CL_SUCCESS)
		error_code = clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_GPU, MAX_DEVICES, devices, &num_gpu_devices);
		CHECK_CL_ERROR(error_code);	
		
		if(m_platformId == INVALID_PLATFORM_ID && num_gpu_devices > 0) 
		{
			m_platformId = platforms[i];
			printf("-----Above platform selected.\n");
		}
	}
	printf("\n");
}
void FluidSolverOpenCL::initialize_stage2_device()
{
	cl_int error_code;
	const size_t MAX_STRING_LENGTH = 1024;
	char string[MAX_STRING_LENGTH];
	
	//Select device
	cl_uint num_devices;
	cl_device_id devices[MAX_DEVICES];
	if(m_platformId != INVALID_PLATFORM_ID)
	{
		//Get devices
		error_code = clGetDeviceIDs(m_platformId, CL_DEVICE_TYPE_GPU, MAX_DEVICES, devices, &num_devices);
		CHECK_CL_ERROR(error_code);	
		
		printf("num_devices: %d\n", num_devices);
		for(cl_uint i = 0; i < num_devices; ++i)
		{
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_NAME, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf("CL_DEVICE_NAME: %s\n", string);
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_PLATFORM, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf( "CL_DEVICE_PLATFORM: %d\n", *reinterpret_cast<const int*>(string) );

			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_VERSION, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf("CL_DEVICE_VERSION: %s\n", string);	
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_OPENCL_C_VERSION, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf("CL_DEVICE_OPENCL_C_VERSION: %s\n", string);	
			
			error_code = clGetDeviceInfo(devices[i], CL_DRIVER_VERSION, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf("CL_DRIVER_VERSION: %s\n", string);	
				
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_PROFILE, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf("CL_DEVICE_PROFILE: %s\n", string);
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_AVAILABLE, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf( "CL_DEVICE_AVAILABLE: %d\n", *reinterpret_cast<const cl_bool*>(string) );	
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_COMPILER_AVAILABLE, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf( "CL_DEVICE_COMPILER_AVAILABLE: %d\n", *reinterpret_cast<const cl_bool*>(string) );
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_ENDIAN_LITTLE, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf( "CL_DEVICE_ENDIAN_LITTLE: %d\n", *reinterpret_cast<const cl_bool*>(string) );
			if(*reinterpret_cast<const cl_bool*>(string) != CL_TRUE) printf(" Warning: device does not use little endian encoding.\n");
			
			error_code = clGetDeviceInfo(devices[i], CL_DEVICE_TYPE, MAX_STRING_LENGTH, string, NULL);
			CHECK_CL_ERROR(error_code);
			printf( "CL_DEVICE_TYPE: %s\n", get_cl_device_type(*reinterpret_cast<cl_device_type*>(string)) );
			
			//Select the first available GPU device
			if(m_device == INVALID_DEVICE_ID)
			{
				m_device = devices[i];
				printf("-----Above device selected.\n");
			}
			
			printf("\n");
		}
		printf("\n");
	}
	else
	{
		printf("FluidSolverOpenCL::initialize_stage2_device() error: invalid m_platformId. \n");
	}
}

void FluidSolverOpenCL::initialize_stage3_context_and_queue()
{
	cl_int error_code;
	
	if(m_device != INVALID_DEVICE_ID)
	{
		//Create context
		cl_context_properties context_properties[3] = { CL_CONTEXT_PLATFORM, cl_context_properties(m_platformId), 0 };
		m_context = clCreateContext(context_properties, 1, &m_device, NULL, NULL, &error_code);
		CHECK_CL_ERROR(error_code);	
		
		//Create command queue
		const cl_command_queue_properties COMMAND_QUEUE_PROPERTIES = 0;		//CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
		m_commandQueue = clCreateCommandQueue(m_context, m_device, COMMAND_QUEUE_PROPERTIES, &error_code);
		CHECK_CL_ERROR(error_code);
	}
	else
	{
		printf("FluidSolverOpenCL::initialize_stage3_context_and_queue() error: invalid m_device. \n");
	}
}

void FluidSolverOpenCL::initialize_stage4_program_and_buffer()
{
	cl_int error_code;
	
	//Load and build program
	if(m_context != INVALID_CONTEXT && m_commandQueue != INVALID_COMMAND_QUEUE)
	{
		const char CL_PROGRAM_PATH[] = "./Demos/FluidDemo/Fluids/OpenCL_support/fluids.cl";
		m_fluidsProgram = compileProgramOpenCL(m_context, m_device, CL_PROGRAM_PATH);
		
		//Kernels
		m_kernel_sphComputePressure = clCreateKernel(m_fluidsProgram, "sphComputePressure", &error_code);
		CHECK_CL_ERROR(error_code);
		m_kernel_sphComputeForce = clCreateKernel(m_fluidsProgram, "sphComputeForce", &error_code);
		CHECK_CL_ERROR(error_code);
		
		//Buffers
		m_buffer_globalFluidParams.allocate( m_context, sizeof(FluidParametersGlobal) );
	}
	else printf("FluidSolverOpenCL::initialize_stage4_program_and_buffers() error: invalid m_context or command_queue. \n");
}

void FluidSolverOpenCL::deactivate_stage1_program_and_buffer()
{
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
	m_fluidsProgram = INVALID_PROGRAM;
	m_kernel_sphComputePressure = INVALID_KERNEL;
	m_kernel_sphComputeForce = INVALID_KERNEL;
}
void FluidSolverOpenCL::deactivate_stage2_context_and_queue()
{
	cl_int error_code;
	
	//Command queues
	error_code = clReleaseCommandQueue(m_commandQueue);
	CHECK_CL_ERROR(error_code);
	
	//Context
	error_code = clReleaseContext(m_context);
	CHECK_CL_ERROR(error_code);
	
	//
	m_commandQueue = INVALID_COMMAND_QUEUE;
	m_context = INVALID_CONTEXT;
}

void FluidSolverOpenCL::sphComputePressure(int numFluidParticles, FluidSortingGrid_OpenCLPointers *gridPointers, 
														Fluid_OpenCLPointers *fluidPointers, btScalar cellSize) 
{
	BT_PROFILE("FluidSolverOpenCL::sphComputePressure()");
	
	cl_int error_code;
	
	///		__global FluidParametersGlobal *FG
	///		__global FluidParametersLocal *FL
	///		__global btVector3 *fluidPosition
	///		__global btScalar *fluidPressure
	///		__global btScalar *fluidInvDensity
	///		__global FluidNeighbors *fluidNeighbors
	///		__global int *numActiveCells
	///		__global SortGridValue *cellValues, 
	///		__global FluidGridIterator *cellContents
	///		btScalar cellSize)
	int arg = 0;
	error_code = clSetKernelArg( m_kernel_sphComputePressure, arg++, sizeof(void*), m_buffer_globalFluidParams.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputePressure, arg++, sizeof(void*), fluidPointers->m_buffer_localParameters );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputePressure, arg++, sizeof(void*), fluidPointers->m_buffer_pos );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputePressure, arg++, sizeof(void*), fluidPointers->m_buffer_pressure );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputePressure, arg++, sizeof(void*), fluidPointers->m_buffer_invDensity );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputePressure, arg++, sizeof(void*), fluidPointers->m_buffer_neighborTable );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputePressure, arg++, sizeof(void*), gridPointers->m_buffer_numActiveCells );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputePressure, arg++, sizeof(void*), gridPointers->m_buffer_activeCells );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputePressure, arg++, sizeof(void*), gridPointers->m_buffer_cellContents );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputePressure, arg++, sizeof(btScalar), static_cast<void*>(&cellSize) );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(m_commandQueue, m_kernel_sphComputePressure, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(m_commandQueue);
	CHECK_CL_ERROR(error_code);
}
void FluidSolverOpenCL::sphComputeForce(int numFluidParticles, FluidSortingGrid_OpenCLPointers *gridPointers, Fluid_OpenCLPointers *fluidPointers) 
{
	BT_PROFILE("FluidSolverOpenCL::sphComputeForce()");
	
	cl_int error_code;
	
	///		__global FluidParametersGlobal *FG
	///		__global FluidParametersLocal *FL
	///		__global btVector3 *fluidPosition
	///		__global btVector3 *fluidVelEval
	///		__global btVector3 *fluidSphForce
	///		__global btScalar *fluidPressure
	///		__global btScalar *fluidInvDensity
	///		__global FluidNeighbors *fluidNeighbors
	int arg = 0;
	error_code = clSetKernelArg( m_kernel_sphComputeForce, arg++, sizeof(void*), m_buffer_globalFluidParams.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputeForce, arg++, sizeof(void*), fluidPointers->m_buffer_localParameters );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputeForce, arg++, sizeof(void*), fluidPointers->m_buffer_pos );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputeForce, arg++, sizeof(void*), fluidPointers->m_buffer_vel_eval );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputeForce, arg++, sizeof(void*), fluidPointers->m_buffer_sph_force );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputeForce, arg++, sizeof(void*), fluidPointers->m_buffer_pressure );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputeForce, arg++, sizeof(void*), fluidPointers->m_buffer_invDensity );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( m_kernel_sphComputeForce, arg++, sizeof(void*), fluidPointers->m_buffer_neighborTable );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(m_commandQueue, m_kernel_sphComputeForce, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(m_commandQueue);
	CHECK_CL_ERROR(error_code);
}
