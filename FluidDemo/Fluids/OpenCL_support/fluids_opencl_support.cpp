/** fluids_opencl_support.cpp
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

#include "fluids_opencl_support.h"

#include "../grid.h"
#include "../fluid.h"

FluidSystem_OpenCL::FluidSystem_OpenCL()
{
	platform_id = INVALID_PLATFORM_ID;
	context = INVALID_CONTEXT;
	
	gpu_device = INVALID_DEVICE_ID;
	gpu_command_queue = INVALID_COMMAND_QUEUE;
	
	fluids_program = INVALID_PROGRAM;
	kernel_grid_insertParticles = INVALID_KERNEL;
	kernel_sph_computePressure = INVALID_KERNEL;
	kernel_sph_computeForce = INVALID_KERNEL;
	kernel_advance = INVALID_KERNEL;
}

void FluidSystem_OpenCL::initialize()
{
	initialize_stage1_platform();
	initialize_stage2_device();
	initialize_stage3_context_and_queue();
	initialize_stage4_program_and_buffers();
}

void FluidSystem_OpenCL::deactivate()
{
	deactivate_stage1_program_and_buffers();
	deactivate_stage2_context_and_queue();
	
	//
	platform_id = INVALID_PLATFORM_ID;
	context = INVALID_CONTEXT;
}

void FluidSystem_OpenCL::stepSimulation(FluidParameters_float *fluidParams, GridParameters *gridParams, 
										Fluids *fluids, int *gridCells, int *gridCellsNumFluids,
										bool transferAllData)
{	
	int numFluidParticles = fluids->size();
	int numGridCells = gridParams->m_numCells;

	//GPU driver may crash if MAX_FLUID_PARTICLES or MAX_GRID_CELLS is exceeded
	if(numFluidParticles > MAX_FLUID_PARTICLES)
	{
		printf("FluidSystem_OpenCL error: numFluidsParticles exceeds MAX_FLUID_PARTICLES(%d > %d).\n", numFluidParticles, MAX_FLUID_PARTICLES);
		return;
	}
	
	if(numGridCells > MAX_GRID_CELLS)
	{
		printf("FluidSystem_OpenCL error: numGridCells exceeds MAX_GRID_CELLS(%d > %d).\n", numGridCells, MAX_GRID_CELLS);
		return;
	}

	//Clear grid before calling grid_insertParticles()
	const int INVALID_PARTICLE_INDEX = -1;
	
	for(int i = 0; i < numGridCells; ++i) gridCells[i] = INVALID_PARTICLE_INDEX;
	for(int i = 0; i < numGridCells; ++i) gridCellsNumFluids[i] = 0;
	
	//
	
	writeToOpencl(fluidParams, gridParams, fluids, gridCells, gridCellsNumFluids, numFluidParticles, numGridCells);
	
	// must write(clear) buffer_gridCells; should also write buffer_gridCellsNumFluids, buffer_externalAcceleration
	//if(transferAllData)writeToOpencl(fluidParams, gridParams, fluids, gridCells, gridCellsNumFluids, numFluidParticles, numGridCells);
	//else writeToOpenclPostFirstFrame(fluidParams, gridParams, fluids, gridCells, gridCellsNumFluids, numFluidParticles, numGridCells);
	
	//
	grid_insertParticles(numFluidParticles);
	sph_computePressure(numFluidParticles);
	sph_computeForce(numFluidParticles);
	advance(numFluidParticles);
	
	//
	readFromOpencl(fluidParams, gridParams, fluids, gridCells, gridCellsNumFluids, numFluidParticles, numGridCells);
}


void FluidSystem_OpenCL::initialize_stage1_platform()
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
	for(int i = 0; i < num_platforms; ++i)
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
void FluidSystem_OpenCL::initialize_stage2_device()
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
		for(int i = 0; i < num_devices; ++i)
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
		printf("FluidSystem_OpenCL::initialize_stage2_device() error: invalid platform_id. \n");
	}
}

void FluidSystem_OpenCL::initialize_stage3_context_and_queue()
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
		printf("FluidSystem_OpenCL::initialize_stage3_context_and_queue() error: invalid gpu_device. \n");
	}
}

void FluidSystem_OpenCL::initialize_stage4_program_and_buffers()
{
	cl_int error_code;
	
	//Load and build program
	if(context != INVALID_CONTEXT && gpu_command_queue != INVALID_COMMAND_QUEUE)
	{
		std::string program_text = load_text_file(CL_PROGRAM_PATH);
		const char *pProgram_data = program_text.c_str();
		size_t program_length = program_text.length();
		
		//Program
		fluids_program = clCreateProgramWithSource(context, 1, const_cast<const char**>(&pProgram_data), NULL, &error_code);
		CHECK_CL_ERROR(error_code);
	
		error_code = clBuildProgram(fluids_program, 1, &gpu_device, NULL, NULL, NULL);
		CHECK_CL_ERROR(error_code);
		
		//if(error_code != CL_SUCCESS)
		{
			const size_t MAX_STRING_LENGTH = 65536;
			char *pString = new char[MAX_STRING_LENGTH];
		
			error_code = clGetDeviceInfo(gpu_device, CL_DEVICE_NAME, MAX_STRING_LENGTH, pString, NULL);
			CHECK_CL_ERROR(error_code);
			printf("for CL_DEVICE_NAME: %s\n", pString);
			
			error_code = clGetProgramBuildInfo( fluids_program, gpu_device, CL_PROGRAM_BUILD_LOG, 
												MAX_STRING_LENGTH, pString, NULL );
			CHECK_CL_ERROR(error_code);
			
			printf("----------CL Program Build Log - start----------------------\n");
			printf("%s\n", pString);
			printf("----------CL Program Build Log - end------------------------\n");
			printf("\n");
			
			delete[] pString;
		}
		
		if(error_code == CL_SUCCESS) printf("%s compiled successfully.\n", CL_PROGRAM_PATH);
	
		//Kernels
		kernel_grid_insertParticles = clCreateKernel(fluids_program, "grid_insertParticles", &error_code);
		CHECK_CL_ERROR(error_code);
		kernel_sph_computePressure = clCreateKernel(fluids_program, "sph_computePressure", &error_code);
		CHECK_CL_ERROR(error_code);
		kernel_sph_computeForce = clCreateKernel(fluids_program, "sph_computeForce", &error_code);
		CHECK_CL_ERROR(error_code);
		kernel_advance = clCreateKernel(fluids_program, "advance", &error_code);
		CHECK_CL_ERROR(error_code);
		
		//Buffers
		buffer_gridParams.allocate( context, sizeof(GridParameters) );
		buffer_fluidParams.allocate( context, sizeof(FluidParameters_float) );
		
		buffer_gridCells.allocate( context, sizeof(int) * MAX_GRID_CELLS );
		buffer_gridCellsNumFluids.allocate( context, sizeof(int) * MAX_GRID_CELLS );
		
		buffer_pos.allocate( context, sizeof(btVector3) * MAX_FLUID_PARTICLES );
		buffer_vel.allocate( context, sizeof(btVector3) * MAX_FLUID_PARTICLES );
		buffer_vel_eval.allocate( context, sizeof(btVector3) * MAX_FLUID_PARTICLES );
		buffer_sph_force.allocate( context, sizeof(btVector3) * MAX_FLUID_PARTICLES );
		buffer_externalAcceleration.allocate( context, sizeof(btVector3) * MAX_FLUID_PARTICLES );
		buffer_prev_pos.allocate( context, sizeof(btVector3) * MAX_FLUID_PARTICLES );
		buffer_pressure.allocate( context, sizeof(float) * MAX_FLUID_PARTICLES );
		buffer_density.allocate( context, sizeof(float) * MAX_FLUID_PARTICLES );
		buffer_nextFluidIndex.allocate( context, sizeof(int) * MAX_FLUID_PARTICLES );
		buffer_neighborTables.allocate( context, sizeof(Neighbors) * MAX_FLUID_PARTICLES );
	}
	else
	{
		printf("FluidSystem_OpenCL::initialize_stage4_program_and_buffers() error: invalid context or command_queue. \n");
	}
}

void FluidSystem_OpenCL::deactivate_stage1_program_and_buffers()
{
	cl_int error_code;
	
	//Buffers
	buffer_gridParams.deallocate();
	buffer_fluidParams.deallocate();
	
	buffer_gridCells.deallocate();
	buffer_gridCellsNumFluids.deallocate();
	
	buffer_pos.deallocate();
	buffer_vel.deallocate();
	buffer_vel_eval.deallocate();
	buffer_sph_force.deallocate();
	buffer_externalAcceleration.deallocate();
	buffer_prev_pos.deallocate();
	buffer_pressure.deallocate();
	buffer_density.deallocate();
	buffer_nextFluidIndex.deallocate();
	buffer_neighborTables.deallocate();
	
	//Kernels
	error_code = clReleaseKernel(kernel_grid_insertParticles);
	CHECK_CL_ERROR(error_code);
	error_code = clReleaseKernel(kernel_sph_computePressure);
	CHECK_CL_ERROR(error_code);
	error_code = clReleaseKernel(kernel_sph_computeForce);
	CHECK_CL_ERROR(error_code);
	error_code = clReleaseKernel(kernel_advance);
	CHECK_CL_ERROR(error_code);
		
	//Program
	error_code = clReleaseProgram(fluids_program);
	CHECK_CL_ERROR(error_code);
	
	//
	fluids_program = INVALID_PROGRAM;
	kernel_grid_insertParticles = INVALID_KERNEL;
	kernel_sph_computePressure = INVALID_KERNEL;
	kernel_sph_computeForce = INVALID_KERNEL;
	kernel_advance = INVALID_KERNEL;
}
void FluidSystem_OpenCL::deactivate_stage2_context_and_queue()
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

void FluidSystem_OpenCL::writeToOpencl( FluidParameters_float *fluidParams, GridParameters *gridParams, 
										Fluids *fluids, int *gridCells, int *gridCellsNumFluids,
										int numFluidParticles, int numGridCells )
{
	BT_PROFILE("FluidSystem_OpenCL::writeToOpencl()");

	cl_int error_code;
	
	buffer_gridParams.writeToBuffer( gpu_command_queue, gridParams, sizeof(GridParameters) );
	buffer_fluidParams.writeToBuffer( gpu_command_queue, fluidParams, sizeof(FluidParameters_float) );
	
	buffer_gridCells.writeToBuffer( gpu_command_queue, gridCells, sizeof(int)*numGridCells );
	buffer_gridCellsNumFluids.writeToBuffer( gpu_command_queue, gridCellsNumFluids, sizeof(int)*numGridCells );
	
	buffer_pos.writeToBuffer( gpu_command_queue, &fluids->m_pos[0], sizeof(btVector3)*numFluidParticles );
	buffer_vel.writeToBuffer( gpu_command_queue, &fluids->m_vel[0], sizeof(btVector3)*numFluidParticles );
	buffer_vel_eval.writeToBuffer( gpu_command_queue, &fluids->m_vel_eval[0], sizeof(btVector3)*numFluidParticles );
	buffer_sph_force.writeToBuffer( gpu_command_queue, &fluids->m_sph_force[0], sizeof(btVector3)*numFluidParticles );
	buffer_externalAcceleration.writeToBuffer( gpu_command_queue, &fluids->m_externalAcceleration[0], sizeof(btVector3)*numFluidParticles );
	buffer_prev_pos.writeToBuffer( gpu_command_queue, &fluids->m_prev_pos[0], sizeof(btVector3)*numFluidParticles );
	buffer_pressure.writeToBuffer( gpu_command_queue, &fluids->m_pressure[0], sizeof(float)*numFluidParticles );
	buffer_density.writeToBuffer( gpu_command_queue, &fluids->m_density[0], sizeof(float)*numFluidParticles );
	buffer_nextFluidIndex.writeToBuffer( gpu_command_queue, &fluids->m_nextFluidIndex[0], sizeof(int)*numFluidParticles );
	buffer_neighborTables.writeToBuffer( gpu_command_queue, &fluids->m_neighborTable[0], sizeof(Neighbors)*numFluidParticles );
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}
void FluidSystem_OpenCL::readFromOpencl( FluidParameters_float *fluidParams, GridParameters *gridParams,
										 Fluids *fluids, int *gridCells, int *gridCellsNumFluids,
										 int numFluidParticles, int numGridCells )
{
	BT_PROFILE("FluidSystem_OpenCL::readFromOpencl()");
	
	cl_int error_code;
	
	buffer_gridParams.readFromBuffer( gpu_command_queue, gridParams, sizeof(GridParameters) );
	buffer_fluidParams.readFromBuffer( gpu_command_queue, fluidParams, sizeof(FluidParameters_float) );
	
	buffer_gridCells.readFromBuffer( gpu_command_queue, gridCells, sizeof(int)*numGridCells );
	buffer_gridCellsNumFluids.readFromBuffer( gpu_command_queue, gridCellsNumFluids, sizeof(int)*numGridCells );
	
	buffer_pos.readFromBuffer( gpu_command_queue, &fluids->m_pos[0], sizeof(btVector3)*numFluidParticles );
	buffer_vel.readFromBuffer( gpu_command_queue, &fluids->m_vel[0], sizeof(btVector3)*numFluidParticles );
	buffer_vel_eval.readFromBuffer( gpu_command_queue, &fluids->m_vel_eval[0], sizeof(btVector3)*numFluidParticles );
	buffer_sph_force.readFromBuffer( gpu_command_queue, &fluids->m_sph_force[0], sizeof(btVector3)*numFluidParticles );
	buffer_externalAcceleration.readFromBuffer( gpu_command_queue, &fluids->m_externalAcceleration[0], sizeof(btVector3)*numFluidParticles );
	buffer_prev_pos.readFromBuffer( gpu_command_queue, &fluids->m_prev_pos[0], sizeof(btVector3)*numFluidParticles );
	buffer_pressure.readFromBuffer( gpu_command_queue, &fluids->m_pressure[0], sizeof(float)*numFluidParticles );
	buffer_density.readFromBuffer( gpu_command_queue, &fluids->m_density[0], sizeof(float)*numFluidParticles );
	buffer_nextFluidIndex.readFromBuffer( gpu_command_queue, &fluids->m_nextFluidIndex[0], sizeof(int)*numFluidParticles );
	buffer_neighborTables.readFromBuffer( gpu_command_queue, &fluids->m_neighborTable[0], sizeof(Neighbors)*numFluidParticles );
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}

void FluidSystem_OpenCL::grid_insertParticles(int numFluidParticles) 
{
	BT_PROFILE("FluidSystem_OpenCL::grid_insertParticles()");

	cl_int error_code;
	///		__global btVector3 *fluidPositions
	///		__global int *fluidNextIndicies
	///		__global GridParameters *gridParams
	///		__global volatile int *gridCells
	///		__global volatile int *gridCellsNumFluids
	error_code = clSetKernelArg( kernel_grid_insertParticles, 0, sizeof(void*), buffer_pos.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_grid_insertParticles, 1, sizeof(void*), buffer_nextFluidIndex.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_grid_insertParticles, 2, sizeof(void*), buffer_gridParams.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_grid_insertParticles, 3, sizeof(void*), buffer_gridCells.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_grid_insertParticles, 4, sizeof(void*), buffer_gridCellsNumFluids.getAddress() );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(gpu_command_queue, kernel_grid_insertParticles, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}
void FluidSystem_OpenCL::sph_computePressure(int numFluidParticles) 
{
	BT_PROFILE("FluidSystem_OpenCL::sph_computePressure()");
	
	cl_int error_code;
	
	///		__global FluidParameters_float *fluidParams
	///		__global btVector3 *fluidPosition
	///		__global float *fluidPressure
	///		__global float *fluidDensity
	///		__global int *fluidNextIndicies
	///		__global Neighbors *fluidNeighbors
	///		__global GridParameters *gridParams
	///		__global int *gridCells
	error_code = clSetKernelArg( kernel_sph_computePressure, 0, sizeof(void*), buffer_fluidParams.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 1, sizeof(void*), buffer_pos.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 2, sizeof(void*), buffer_pressure.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 3, sizeof(void*), buffer_density.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 4, sizeof(void*), buffer_nextFluidIndex.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 5, sizeof(void*), buffer_neighborTables.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 6, sizeof(void*), buffer_gridParams.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 7, sizeof(void*), buffer_gridCells.getAddress() );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(gpu_command_queue, kernel_sph_computePressure, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}
void FluidSystem_OpenCL::sph_computeForce(int numFluidParticles) 
{
	BT_PROFILE("FluidSystem_OpenCL::sph_computeForce()");
	
	cl_int error_code;
	
	///		__global FluidParameters_float *fluidParams
	///		__global btVector3 *fluidPosition
	///		__global btVector3 *fluidVelEval
	///		__global btVector3 *fluidSphForce
	///		__global float *fluidPressure
	///		__global float *fluidDensity
	///		__global Neighbors *fluidNeighbors
	error_code = clSetKernelArg( kernel_sph_computeForce, 0, sizeof(void*), buffer_fluidParams.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 1, sizeof(void*), buffer_pos.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 2, sizeof(void*), buffer_vel_eval.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 3, sizeof(void*), buffer_sph_force.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 4, sizeof(void*), buffer_pressure.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 5, sizeof(void*), buffer_density.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 6, sizeof(void*), buffer_neighborTables.getAddress() );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(gpu_command_queue, kernel_sph_computeForce, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}
void FluidSystem_OpenCL::advance(int numFluidParticles) 
{
	BT_PROFILE("FluidSystem_OpenCL::advance()");
	
	cl_int error_code;
	
	///		__global FluidParameters_float *fluidParams,
	///		__global btVector3 *fluidPosition
	///		__global btVector3 *fluidVel
	///		__global btVector3 *fluidVelEval
	///		__global btVector3 *fluidSphForce,
	///		__global btVector3 *fluidExternalAcceleration
	///		__global btVector3 *fluidPrevPosition
	error_code = clSetKernelArg( kernel_advance, 0, sizeof(void*), buffer_fluidParams.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_advance, 1, sizeof(void*), buffer_pos.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_advance, 2, sizeof(void*), buffer_vel.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_advance, 3, sizeof(void*), buffer_vel_eval.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_advance, 4, sizeof(void*), buffer_sph_force.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_advance, 5, sizeof(void*), buffer_externalAcceleration.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_advance, 6, sizeof(void*), buffer_prev_pos.getAddress() );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(gpu_command_queue, kernel_advance, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}
