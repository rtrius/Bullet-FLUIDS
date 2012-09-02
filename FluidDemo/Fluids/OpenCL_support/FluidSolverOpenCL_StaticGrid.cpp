/* FluidSolverOpenCL_StaticGrid.cpp
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
#include "FluidSolverOpenCL_StaticGrid.h"

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

#include "../FluidStaticGrid.h"
#include "../FluidParameters.h"

FluidSolverOpenCL_StaticGrid::FluidSolverOpenCL_StaticGrid()
{
	platform_id = INVALID_PLATFORM_ID;
	context = INVALID_CONTEXT;
	
	gpu_device = INVALID_DEVICE_ID;
	gpu_command_queue = INVALID_COMMAND_QUEUE;
	
	fluids_program = INVALID_PROGRAM;
	kernel_grid_insertParticles = INVALID_KERNEL;
	kernel_sph_computePressure = INVALID_KERNEL;
	kernel_sph_computeForce = INVALID_KERNEL;
	kernel_integrate = INVALID_KERNEL;
	
	//
	initialize();
}

void FluidSolverOpenCL_StaticGrid::initialize()
{
	initialize_stage1_platform();
	initialize_stage2_device();
	initialize_stage3_context_and_queue();
	initialize_stage4_program_and_buffer();
}

void FluidSolverOpenCL_StaticGrid::deactivate()
{
	deactivate_stage1_program_and_buffer();
	deactivate_stage2_context_and_queue();
	
	//
	platform_id = INVALID_PLATFORM_ID;
	context = INVALID_CONTEXT;
}

void FluidSolverOpenCL_StaticGrid::stepSimulation(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids)
{	
	BT_PROFILE("FluidSolverOpenCL_StaticGrid::stepSimulation()");
	
#ifdef BT_USE_DOUBLE_PRECISION
		printf("BT_USE_DOUBLE_PRECISION not supported on OpenCL.\n");
		return;
#endif	

	//FluidSolverOpenCL_StaticGrid requires use of 'class FluidStaticGrid'
	btAlignedObjectArray<FluidSph*> validFluids;
	for(int i = 0; i < fluids->size(); ++i) 
	{
		//GPU driver may crash if grid is incorrect type
		if( !(*fluids)[i]->numParticles() || (*fluids)[i]->getGrid()->getGridType() != FluidGrid::FT_LinkedList ) continue;
		
		validFluids.push_back( (*fluids)[i] );
	}

	int numValidFluids = validFluids.size();
	
	//Grids must be cleared before calling grid_insertParticles()
	for(int i = 0; i < numValidFluids; ++i)validFluids[i]->internalGetGrid()->clear();
	
	//Write data to OpenCL
	FluidParametersGlobal globalParameters = FG;
	buffer_globalFluidParams.writeToBuffer( gpu_command_queue, &globalParameters, sizeof(FluidParametersGlobal) );
	cl_int error_code = clFinish(gpu_command_queue);
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
		for(int i = 0; i < numValidFluids; ++i) writeSingleFluidToOpenCL(validFluids[i], &m_fluidData[i], &m_gridData[i]);
	}
	
	//stepSimulation()
	{
		BT_PROFILE("stepSimulation() - grid, sph, integrate");
	
		for(int i = 0; i < numValidFluids; ++i)
		{
			int numFluidParticles = validFluids[i]->numParticles();
			FluidStaticGrid_OpenCLPointers gridPointers = m_gridData[i].getPointers();
			Fluid_OpenCLPointers fluidPointers = m_fluidData[i].getPointers();
			
			grid_insertParticles(numFluidParticles, &gridPointers, &fluidPointers);
			sph_computePressure(numFluidParticles, &gridPointers, &fluidPointers);
			sph_computeForce(numFluidParticles, &gridPointers, &fluidPointers);
			integrate(numFluidParticles, &fluidPointers);
		}
	}
	
	//Read data from OpenCL
	{
		BT_PROFILE("stepSimulation() - readFromOpenCL");
		for(int i = 0; i < numValidFluids; ++i) readSingleFluidFromOpenCL(validFluids[i], &m_fluidData[i], &m_gridData[i]);
	}
	
	//Clear externalAcceleration, instead of reading it from OpenCL( should all be (0,0,0) after integrate() )
	for(int i = 0; i < numValidFluids; ++i)
	{
		FluidParticles &particles = validFluids[i]->internalGetFluidParticles();
		for(int n = 0; n < particles.size(); ++n) particles.m_externalAcceleration[n].setValue(0,0,0);
	}
}

void FluidSolverOpenCL_StaticGrid::initialize_stage1_platform()
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
		
		if(platform_id == INVALID_PLATFORM_ID && num_gpu_devices > 0) 
		{
			platform_id = platforms[i];
			printf("-----Above platform selected.\n");
		}
	}
	printf("\n");
}
void FluidSolverOpenCL_StaticGrid::initialize_stage2_device()
{
	cl_int error_code;
	const size_t MAX_STRING_LENGTH = 1024;
	char string[MAX_STRING_LENGTH];
	
	//Select device
	cl_uint num_devices;
	cl_device_id devices[MAX_DEVICES];
	if(platform_id != INVALID_PLATFORM_ID)
	{
		//Get devices
		error_code = clGetDeviceIDs(platform_id, CL_DEVICE_TYPE_GPU, MAX_DEVICES, devices, &num_devices);
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
			if(gpu_device == INVALID_DEVICE_ID)
			{
				gpu_device = devices[i];
				printf("-----Above device selected.\n");
			}
			
			printf("\n");
		}
		printf("\n");
	}
	else
	{
		printf("FluidSolverOpenCL_StaticGrid::initialize_stage2_device() error: invalid platform_id. \n");
	}
}

void FluidSolverOpenCL_StaticGrid::initialize_stage3_context_and_queue()
{
	cl_int error_code;
	
	if(gpu_device != INVALID_DEVICE_ID)
	{
		//Create context
		cl_context_properties context_properties[3] = { CL_CONTEXT_PLATFORM, cl_context_properties(platform_id), 0 };
		context = clCreateContext(context_properties, 1, &gpu_device, NULL, NULL, &error_code);
		CHECK_CL_ERROR(error_code);	
		
		//Create command queue
		const cl_command_queue_properties COMMAND_QUEUE_PROPERTIES = 0;		//CL_QUEUE_OUT_OF_ORDER_EXEC_MODE_ENABLE;
		gpu_command_queue = clCreateCommandQueue(context, gpu_device, COMMAND_QUEUE_PROPERTIES, &error_code);
		CHECK_CL_ERROR(error_code);
	}
	else
	{
		printf("FluidSolverOpenCL_StaticGrid::initialize_stage3_context_and_queue() error: invalid gpu_device. \n");
	}
}

void FluidSolverOpenCL_StaticGrid::initialize_stage4_program_and_buffer()
{
	cl_int error_code;
	
	//Load and build program
	if(context != INVALID_CONTEXT && gpu_command_queue != INVALID_COMMAND_QUEUE)
	{
		const char CL_PROGRAM_PATH[] = "./Demos/FluidDemo/Fluids/OpenCL_support/fluids_staticGrid.cl";
		fluids_program = compileProgramOpenCL(context, gpu_device, CL_PROGRAM_PATH);
		
		//Kernels
		kernel_grid_insertParticles = clCreateKernel(fluids_program, "grid_insertParticles", &error_code);
		CHECK_CL_ERROR(error_code);
		kernel_sph_computePressure = clCreateKernel(fluids_program, "sph_computePressure", &error_code);
		CHECK_CL_ERROR(error_code);
		kernel_sph_computeForce = clCreateKernel(fluids_program, "sph_computeForce", &error_code);
		CHECK_CL_ERROR(error_code);
		kernel_integrate = clCreateKernel(fluids_program, "integrate", &error_code);
		CHECK_CL_ERROR(error_code);
		
		//Buffers
		buffer_globalFluidParams.allocate( context, sizeof(FluidParametersGlobal) );
	}
	else printf("FluidSolverOpenCL_StaticGrid::initialize_stage4_program_and_buffers() error: invalid context or command_queue. \n");
}

void FluidSolverOpenCL_StaticGrid::deactivate_stage1_program_and_buffer()
{
	cl_int error_code;
	
	//Buffers
	buffer_globalFluidParams.deallocate();
	
	//Kernels
	error_code = clReleaseKernel(kernel_grid_insertParticles);
	CHECK_CL_ERROR(error_code);
	error_code = clReleaseKernel(kernel_sph_computePressure);
	CHECK_CL_ERROR(error_code);
	error_code = clReleaseKernel(kernel_sph_computeForce);
	CHECK_CL_ERROR(error_code);
	error_code = clReleaseKernel(kernel_integrate);
	CHECK_CL_ERROR(error_code);
		
	//Program
	error_code = clReleaseProgram(fluids_program);
	CHECK_CL_ERROR(error_code);
	
	//
	fluids_program = INVALID_PROGRAM;
	kernel_grid_insertParticles = INVALID_KERNEL;
	kernel_sph_computePressure = INVALID_KERNEL;
	kernel_sph_computeForce = INVALID_KERNEL;
	kernel_integrate = INVALID_KERNEL;
}
void FluidSolverOpenCL_StaticGrid::deactivate_stage2_context_and_queue()
{
	cl_int error_code;
	
	//Command queues
	error_code = clReleaseCommandQueue(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
	
	//Context
	error_code = clReleaseContext(context);
	CHECK_CL_ERROR(error_code);
	
	//
	gpu_command_queue = INVALID_COMMAND_QUEUE;
	context = INVALID_CONTEXT;
}

void FluidSolverOpenCL_StaticGrid::writeSingleFluidToOpenCL(FluidSph* fluid, Fluid_OpenCL *fluidData, FluidStaticGrid_OpenCL *gridData)
{
	const bool WRITE_ALL_DATA = false;
	const FluidParametersLocal &FL = fluid->getLocalParameters();
		
	gridData->writeToOpenCL( context, gpu_command_queue, static_cast<FluidStaticGrid*>(fluid->internalGetGrid()) );
	fluidData->writeToOpenCL( context, gpu_command_queue, FL, &fluid->internalGetFluidParticles(), WRITE_ALL_DATA );
}
void FluidSolverOpenCL_StaticGrid::readSingleFluidFromOpenCL(FluidSph* fluid, Fluid_OpenCL *fluidData, FluidStaticGrid_OpenCL *gridData)
{	
	const bool READ_ALL_DATA = false;
	
	gridData->readFromOpenCL( context, gpu_command_queue, static_cast<FluidStaticGrid*>(fluid->internalGetGrid()) );
	fluidData->readFromOpenCL( context, gpu_command_queue, &fluid->internalGetFluidParticles(), READ_ALL_DATA );
}

void FluidSolverOpenCL_StaticGrid::grid_insertParticles(int numFluidParticles, FluidStaticGrid_OpenCLPointers *gridPointers, 
											 Fluid_OpenCLPointers *fluidPointers) 
{
	BT_PROFILE("FluidSolverOpenCL_StaticGrid::grid_insertParticles()");

	cl_int error_code;
	///		__global btVector3 *fluidPositions
	///		__global int *fluidNextIndicies
	///		__global FluidStaticGridParameters *gridParams
	///		__global volatile int *gridCells
	///		__global volatile int *gridCellsNumFluids
	error_code = clSetKernelArg( kernel_grid_insertParticles, 0, sizeof(void*), fluidPointers->m_buffer_pos );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_grid_insertParticles, 1, sizeof(void*), fluidPointers->m_buffer_nextFluidIndex );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_grid_insertParticles, 2, sizeof(void*), gridPointers->m_buffer_gridParams );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_grid_insertParticles, 3, sizeof(void*), gridPointers->m_buffer_gridCells );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_grid_insertParticles, 4, sizeof(void*), gridPointers->m_buffer_gridCellsNumFluids );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(gpu_command_queue, kernel_grid_insertParticles, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}
void FluidSolverOpenCL_StaticGrid::sph_computePressure(int numFluidParticles, FluidStaticGrid_OpenCLPointers *gridPointers, 
											Fluid_OpenCLPointers *fluidPointers) 
{
	BT_PROFILE("FluidSolverOpenCL_StaticGrid::sph_computePressure()");
	
	cl_int error_code;
	
	///		__global FluidParametersGlobal *FG
	///		__global FluidParametersLocal *FL
	///		__global btVector3 *fluidPosition
	///		__global btScalar *fluidPressure
	///		__global btScalar *fluidInvDensity
	///		__global int *fluidNextIndicies
	///		__global FluidNeighbors *fluidNeighbors
	///		__global FluidStaticGridParameters *gridParams
	///		__global int *gridCells
	error_code = clSetKernelArg( kernel_sph_computePressure, 0, sizeof(void*), buffer_globalFluidParams.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 1, sizeof(void*), fluidPointers->m_buffer_localParameters );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 2, sizeof(void*), fluidPointers->m_buffer_pos );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 3, sizeof(void*), fluidPointers->m_buffer_pressure );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 4, sizeof(void*), fluidPointers->m_buffer_invDensity );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 5, sizeof(void*), fluidPointers->m_buffer_nextFluidIndex );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 6, sizeof(void*), fluidPointers->m_buffer_neighborTable );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 7, sizeof(void*), gridPointers->m_buffer_gridParams );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 8, sizeof(void*), gridPointers->m_buffer_gridCells );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(gpu_command_queue, kernel_sph_computePressure, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}
void FluidSolverOpenCL_StaticGrid::sph_computeForce(int numFluidParticles, FluidStaticGrid_OpenCLPointers *gridPointers, Fluid_OpenCLPointers *fluidPointers) 
{
	BT_PROFILE("FluidSolverOpenCL_StaticGrid::sph_computeForce()");
	
	cl_int error_code;
	
	///		__global FluidParametersGlobal *FG
	///		__global FluidParametersLocal *FL
	///		__global btVector3 *fluidPosition
	///		__global btVector3 *fluidVelEval
	///		__global btVector3 *fluidSphForce
	///		__global btScalar *fluidPressure
	///		__global btScalar *fluidInvDensity
	///		__global FluidNeighbors *fluidNeighbors
	error_code = clSetKernelArg( kernel_sph_computeForce, 0, sizeof(void*), buffer_globalFluidParams.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 1, sizeof(void*), fluidPointers->m_buffer_localParameters );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 2, sizeof(void*), fluidPointers->m_buffer_pos );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 3, sizeof(void*), fluidPointers->m_buffer_vel_eval );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 4, sizeof(void*), fluidPointers->m_buffer_sph_force );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 5, sizeof(void*), fluidPointers->m_buffer_pressure );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 6, sizeof(void*), fluidPointers->m_buffer_invDensity );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 7, sizeof(void*), fluidPointers->m_buffer_neighborTable );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(gpu_command_queue, kernel_sph_computeForce, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}
void FluidSolverOpenCL_StaticGrid::integrate(int numFluidParticles, Fluid_OpenCLPointers *fluidPointers) 
{
	BT_PROFILE("FluidSolverOpenCL_StaticGrid::integrate()");
	
	cl_int error_code;
	
	///		__global FluidParametersGlobal *FG
	///		__global FluidParametersLocal *FL
	///		__global btVector3 *fluidPosition
	///		__global btVector3 *fluidVel
	///		__global btVector3 *fluidVelEval
	///		__global btVector3 *fluidSphForce
	///		__global btVector3 *fluidExternalAcceleration
	error_code = clSetKernelArg( kernel_integrate, 0, sizeof(void*), buffer_globalFluidParams.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_integrate, 1, sizeof(void*), fluidPointers->m_buffer_localParameters );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_integrate, 2, sizeof(void*), fluidPointers->m_buffer_pos );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_integrate, 3, sizeof(void*), fluidPointers->m_buffer_vel );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_integrate, 4, sizeof(void*), fluidPointers->m_buffer_vel_eval );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_integrate, 5, sizeof(void*), fluidPointers->m_buffer_sph_force );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_integrate, 6, sizeof(void*), fluidPointers->m_buffer_externalAcceleration );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(gpu_command_queue, kernel_integrate, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}

