/** FluidSolverOpenCL.cpp
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

#include "../grid.h"
#include "../FluidParameters.h"
#include "../FluidParticles.h"
#include "../FluidSph.h"

FluidSolverOpenCL::FluidSolverOpenCL()
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
	platform_id = INVALID_PLATFORM_ID;
	context = INVALID_CONTEXT;
}

void FluidSolverOpenCL::stepSimulation(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids)
{	
	BT_PROFILE("FluidSolverOpenCL::stepSimulation()");
	
#ifdef BT_USE_DOUBLE_PRECISION
		printf("BT_USE_DOUBLE_PRECISION not supported on OpenCL.\n");
		return;
#endif	

	//FluidSolverOpenCL requires use of 'class Grid'
	btAlignedObjectArray<FluidSph*> validFluids;
	for(int i = 0; i < fluids->size(); ++i) 
	{
		//GPU driver may crash if grid is incorrect type
		if( !(*fluids)[i]->numParticles() || (*fluids)[i]->getGrid()->getGridType() != FT_LinkedList ) continue;
			
		validFluids.push_back( (*fluids)[i] );
	}

	int numValidFluids = validFluids.size();
	
	//Clear grids before calling grid_insertParticles()
	for(int i = 0; i < numValidFluids; ++i)validFluids[i]->internalGetGrid()->clear();
	
	//Write data to OpenCL
	FluidParametersGlobal globalParameters = FG;
	buffer_globalFluidParams.writeToBuffer( gpu_command_queue, &globalParameters, sizeof(FluidParametersGlobal) );
	cl_int error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
	
	m_gridData.resize(numValidFluids);
	m_fluidData.resize(numValidFluids);

	const bool TRANSFER_ALL_DATA = false;
	for(int i = 0; i < numValidFluids; ++i)
	{
		FluidParametersLocal FL = validFluids[i]->getLocalParameters();
	
		m_gridData[i].writeToOpenCl( context, gpu_command_queue, static_cast<Grid*>(validFluids[i]->internalGetGrid()) );
		m_fluidData[i].writeToOpenCl( context, gpu_command_queue, &FL, &validFluids[i]->internalGetFluidParticles(), TRANSFER_ALL_DATA );
	}
	
	//stepSimulation()
	for(int i = 0; i < numValidFluids; ++i)
	{
		int numFluidParticles = validFluids[i]->numParticles();
		Grid_OpenCLPointers gridPointers = m_gridData[i].getPointers();
		Fluid_OpenCLPointers fluidPointers = m_fluidData[i].getPointers();
		
		grid_insertParticles(numFluidParticles, &gridPointers, &fluidPointers);
		sph_computePressure(numFluidParticles, &gridPointers, &fluidPointers);
		sph_computeForce(numFluidParticles, &gridPointers, &fluidPointers);
		integrate(numFluidParticles, &fluidPointers);
	}
	
	//Read data from OpenCL
	for(int i = 0; i < numValidFluids; ++i)
	{
		m_gridData[i].readFromOpenCl( context, gpu_command_queue, static_cast<Grid*>(validFluids[i]->internalGetGrid()) );
		m_fluidData[i].readFromOpenCl( context, gpu_command_queue, &validFluids[i]->internalGetFluidParticles(), TRANSFER_ALL_DATA );
	}
	
	//Clear externalAcceleration, instead of reading it from OpenCL( should all be (0,0,0) after integrate() )
	for(int i = 0; i < numValidFluids; ++i)
	{
		FluidParticles &particles = validFluids[i]->internalGetFluidParticles();
		for(int n = 0; n < particles.size(); ++n) particles.m_externalAcceleration[n].setValue(0,0,0);
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
void FluidSolverOpenCL::initialize_stage2_device()
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
		printf("FluidSolverOpenCL::initialize_stage2_device() error: invalid platform_id. \n");
	}
}

void FluidSolverOpenCL::initialize_stage3_context_and_queue()
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
		printf("FluidSolverOpenCL::initialize_stage3_context_and_queue() error: invalid gpu_device. \n");
	}
}

void FluidSolverOpenCL::initialize_stage4_program_and_buffer()
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
		kernel_integrate = clCreateKernel(fluids_program, "integrate", &error_code);
		CHECK_CL_ERROR(error_code);
		
		//Buffers
		buffer_globalFluidParams.allocate( context, sizeof(FluidParametersGlobal) );
		
	}
	else
	{
		printf("FluidSolverOpenCL::initialize_stage4_program_and_buffers() error: invalid context or command_queue. \n");
	}
}

void FluidSolverOpenCL::deactivate_stage1_program_and_buffer()
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
void FluidSolverOpenCL::deactivate_stage2_context_and_queue()
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


void FluidSolverOpenCL::grid_insertParticles(int numFluidParticles, Grid_OpenCLPointers *gridPointers, Fluid_OpenCLPointers *fluidPointers) 
{
	BT_PROFILE("FluidSolverOpenCL::grid_insertParticles()");

	cl_int error_code;
	///		__global btVector3 *fluidPositions
	///		__global int *fluidNextIndicies
	///		__global GridParameters *gridParams
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
void FluidSolverOpenCL::sph_computePressure(int numFluidParticles, Grid_OpenCLPointers *gridPointers, Fluid_OpenCLPointers *fluidPointers) 
{
	BT_PROFILE("FluidSolverOpenCL::sph_computePressure()");
	
	cl_int error_code;
	
	///		__global FluidParametersGlobal *FG
	///		__global FluidParametersLocal *FL
	///		__global btVector3 *fluidPosition
	///		__global btScalar *fluidPressure
	///		__global btScalar *fluidDensity
	///		__global int *fluidNextIndicies
	///		__global Neighbors *fluidNeighbors
	///		__global GridParameters *gridParams
	///		__global int *gridCells
	error_code = clSetKernelArg( kernel_sph_computePressure, 0, sizeof(void*), buffer_globalFluidParams.getAddress() );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 1, sizeof(void*), fluidPointers->m_buffer_localParameters );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 2, sizeof(void*), fluidPointers->m_buffer_pos );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 3, sizeof(void*), fluidPointers->m_buffer_pressure );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computePressure, 4, sizeof(void*), fluidPointers->m_buffer_density );
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
void FluidSolverOpenCL::sph_computeForce(int numFluidParticles, Grid_OpenCLPointers *gridPointers, Fluid_OpenCLPointers *fluidPointers) 
{
	BT_PROFILE("FluidSolverOpenCL::sph_computeForce()");
	
	cl_int error_code;
	
	///		__global FluidParametersGlobal *FG
	///		__global FluidParametersLocal *FL
	///		__global btVector3 *fluidPosition
	///		__global btVector3 *fluidVelEval
	///		__global btVector3 *fluidSphForce
	///		__global btScalar *fluidPressure
	///		__global btScalar *fluidDensity
	///		__global Neighbors *fluidNeighbors
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
	error_code = clSetKernelArg( kernel_sph_computeForce, 6, sizeof(void*), fluidPointers->m_buffer_density );
	CHECK_CL_ERROR(error_code);
	error_code = clSetKernelArg( kernel_sph_computeForce, 7, sizeof(void*), fluidPointers->m_buffer_neighborTable );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(gpu_command_queue, kernel_sph_computeForce, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}
void FluidSolverOpenCL::integrate(int numFluidParticles, Fluid_OpenCLPointers *fluidPointers) 
{
	BT_PROFILE("FluidSolverOpenCL::integrate()");
	
	cl_int error_code;
	
	///		__global FluidParametersGlobal *FG
	///		__global FluidParametersLocal *FL
	///		__global btVector3 *fluidPosition
	///		__global btVector3 *fluidVel
	///		__global btVector3 *fluidVelEval
	///		__global btVector3 *fluidSphForce
	///		__global btVector3 *fluidExternalAcceleration
	///		__global btVector3 *fluidPrevPosition
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
	error_code = clSetKernelArg( kernel_integrate, 7, sizeof(void*), fluidPointers->m_buffer_prev_pos );
	CHECK_CL_ERROR(error_code);
	
	size_t global_work_size = numFluidParticles;
	error_code = clEnqueueNDRangeKernel(gpu_command_queue, kernel_integrate, 1, NULL, &global_work_size, NULL, 0, NULL, NULL);
	CHECK_CL_ERROR(error_code);
	
	error_code = clFinish(gpu_command_queue);
	CHECK_CL_ERROR(error_code);
}


////////////////////////////////////////////////////////////////////////////////
/// class Grid_OpenCL
////////////////////////////////////////////////////////////////////////////////
void Grid_OpenCL::writeToOpenCl(cl_context context, cl_command_queue commandQueue, Grid *grid)
{
	GridParameters GP = grid->getParameters();
	
	int currentNumCells = GP.m_numCells;
	if(m_numGridCells != currentNumCells)
	{
		deallocate();
		allocate(context, currentNumCells);
	}
	
	m_buffer_gridParams.writeToBuffer( commandQueue, &GP, sizeof(GridParameters) );
	
	m_buffer_gridCells.writeToBuffer( commandQueue, grid->getCellsPointer(), sizeof(int)*currentNumCells );
	m_buffer_gridCellsNumFluids.writeToBuffer( commandQueue, grid->getCellsNumFluidsPointer(), sizeof(int)*currentNumCells );

	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}
void Grid_OpenCL::readFromOpenCl(cl_context context, cl_command_queue commandQueue, Grid *grid)
{
	int currentNumCells = grid->getParameters().m_numCells;
	
	m_buffer_gridCells.readFromBuffer( commandQueue, grid->getCellsPointer(), sizeof(int)*currentNumCells );
	m_buffer_gridCellsNumFluids.readFromBuffer( commandQueue, grid->getCellsNumFluidsPointer(), sizeof(int)*currentNumCells );

	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}

Grid_OpenCLPointers Grid_OpenCL::getPointers()
{
	Grid_OpenCLPointers result;
	
	result.m_buffer_gridParams = m_buffer_gridParams.getAddress();
	
	result.m_buffer_gridCells = m_buffer_gridCells.getAddress();
	result.m_buffer_gridCellsNumFluids = m_buffer_gridCellsNumFluids.getAddress();

	return result;
}

void Grid_OpenCL::allocate(cl_context context, int numGridCells)
{
	m_buffer_gridParams.allocate( context, sizeof(GridParameters) );
	
	m_numGridCells = numGridCells;
	m_buffer_gridCells.allocate( context, sizeof(int) * numGridCells );
	m_buffer_gridCellsNumFluids.allocate( context, sizeof(int) * numGridCells );
}
void Grid_OpenCL::deallocate()
{
	m_buffer_gridParams.deallocate();
	
	m_numGridCells = 0;
	m_buffer_gridCells.deallocate();
	m_buffer_gridCellsNumFluids.deallocate();
}


////////////////////////////////////////////////////////////////////////////////
/// class Fluid_OpenCL
////////////////////////////////////////////////////////////////////////////////
void Fluid_OpenCL::writeToOpenCl(cl_context context, cl_command_queue commandQueue, 
								 FluidParametersLocal *localParameters, FluidParticles *particles, bool transferAllData)
{
	if(m_maxParticles != particles->m_maxParticles)
	{
		deallocate();
		allocate(context, particles->m_maxParticles);
	}
	
	m_buffer_localParameters.writeToBuffer( commandQueue, localParameters, sizeof(FluidParametersLocal) );
	
	int numParticles = particles->size();
	
	m_buffer_pos.writeToBuffer( commandQueue, &particles->m_pos[0], sizeof(btVector3)*numParticles );
	m_buffer_vel.writeToBuffer( commandQueue, &particles->m_vel[0], sizeof(btVector3)*numParticles );
	m_buffer_vel_eval.writeToBuffer( commandQueue, &particles->m_vel_eval[0], sizeof(btVector3)*numParticles );
	
	m_buffer_externalAcceleration.writeToBuffer( commandQueue, &particles->m_externalAcceleration[0], sizeof(btVector3)*numParticles );
	
	if(transferAllData)
	{
		m_buffer_sph_force.writeToBuffer( commandQueue, &particles->m_sph_force[0], sizeof(btVector3)*numParticles );
		m_buffer_prev_pos.writeToBuffer( commandQueue, &particles->m_prev_pos[0], sizeof(btVector3)*numParticles );
		m_buffer_pressure.writeToBuffer( commandQueue, &particles->m_pressure[0], sizeof(btScalar)*numParticles );
		m_buffer_density.writeToBuffer( commandQueue, &particles->m_density[0], sizeof(btScalar)*numParticles );
		m_buffer_nextFluidIndex.writeToBuffer( commandQueue, &particles->m_nextFluidIndex[0], sizeof(int)*numParticles );
		m_buffer_neighborTable.writeToBuffer( commandQueue, &particles->m_neighborTable[0], sizeof(Neighbors)*numParticles );
	}
	
	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}

void Fluid_OpenCL::readFromOpenCl(cl_context context, cl_command_queue commandQueue,
								  FluidParticles *particles, bool transferAllData)
{
	int numParticles = particles->size();

	m_buffer_pos.readFromBuffer( commandQueue, &particles->m_pos[0], sizeof(btVector3)*numParticles );
	m_buffer_vel.readFromBuffer( commandQueue, &particles->m_vel[0], sizeof(btVector3)*numParticles );
	m_buffer_vel_eval.readFromBuffer( commandQueue, &particles->m_vel_eval[0], sizeof(btVector3)*numParticles );
	
	m_buffer_prev_pos.readFromBuffer( commandQueue, &particles->m_prev_pos[0], sizeof(btVector3)*numParticles );
	m_buffer_nextFluidIndex.readFromBuffer( commandQueue, &particles->m_nextFluidIndex[0], sizeof(int)*numParticles );
	
	if(transferAllData)
	{
		m_buffer_sph_force.readFromBuffer( commandQueue, &particles->m_sph_force[0], sizeof(btVector3)*numParticles );
		m_buffer_externalAcceleration.readFromBuffer( commandQueue, &particles->m_externalAcceleration[0], sizeof(btVector3)*numParticles );
		m_buffer_pressure.readFromBuffer( commandQueue, &particles->m_pressure[0], sizeof(btScalar)*numParticles );
		m_buffer_density.readFromBuffer( commandQueue, &particles->m_density[0], sizeof(btScalar)*numParticles );
		m_buffer_neighborTable.readFromBuffer( commandQueue, &particles->m_neighborTable[0], sizeof(Neighbors)*numParticles );
	}
	
	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}

Fluid_OpenCLPointers Fluid_OpenCL::getPointers()
{
	Fluid_OpenCLPointers result;
	
	result.m_buffer_localParameters = m_buffer_localParameters.getAddress();
	
	result.m_buffer_pos = m_buffer_pos.getAddress();
	result.m_buffer_vel = m_buffer_vel.getAddress();
	result.m_buffer_vel_eval = m_buffer_vel_eval.getAddress();
	result.m_buffer_sph_force = m_buffer_sph_force.getAddress();
	result.m_buffer_externalAcceleration = m_buffer_externalAcceleration.getAddress();
	result.m_buffer_prev_pos = m_buffer_prev_pos.getAddress();
	result.m_buffer_pressure = m_buffer_pressure.getAddress();
	result.m_buffer_density = m_buffer_density.getAddress();
	result.m_buffer_nextFluidIndex = m_buffer_nextFluidIndex.getAddress();
	
	result.m_buffer_neighborTable = m_buffer_neighborTable.getAddress();

	return result;
}

void Fluid_OpenCL::allocate(cl_context context, int maxParticles)
{
	m_maxParticles = maxParticles;
	
	m_buffer_localParameters.allocate( context, sizeof(FluidParametersLocal) );
	
	m_buffer_pos.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_vel.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_vel_eval.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_sph_force.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_externalAcceleration.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_prev_pos.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_pressure.allocate( context, sizeof(btScalar) * maxParticles );
	m_buffer_density.allocate( context, sizeof(btScalar) * maxParticles );
	m_buffer_nextFluidIndex.allocate( context, sizeof(int) * maxParticles );
	
	m_buffer_neighborTable.allocate( context, sizeof(Neighbors) * maxParticles );
}
void Fluid_OpenCL::deallocate()
{	
	m_maxParticles = 0;
	
	m_buffer_localParameters.deallocate();
	
	m_buffer_pos.deallocate();
	m_buffer_vel.deallocate();
	m_buffer_vel_eval.deallocate();
	m_buffer_sph_force.deallocate();
	m_buffer_externalAcceleration.deallocate();
	m_buffer_prev_pos.deallocate();
	m_buffer_pressure.deallocate();
	m_buffer_density.deallocate();
	m_buffer_nextFluidIndex.deallocate();
	
	m_buffer_neighborTable.deallocate();
}
