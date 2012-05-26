/** fluids_opencl_support.h
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

#ifndef FLUIDS_OPENCL_H_INCLUDED
#define FLUIDS_OPENCL_H_INCLUDED

#include <string>
#include <cstdio>

#include <CL/cl.h>

#include "opencl_support.h"

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro


struct FluidParameters_float;
struct GridParameters;
struct Fluid;
class Neighbors;

//Loaded into FluidSystem_OpenCL.fluids_program
const char CL_PROGRAM_PATH[] = "./Demos/FluidDemo/Fluids/OpenCL_support/fluids.cl";
		
class FluidSystem_OpenCL
{
	static const int MAX_FLUID_PARTICLES = 32768;	//Determines size of buffer_fluids and buffer_neighborTables
	static const int MAX_GRID_CELLS = 65536;		//Determines size of buffer_gridCells and buffer_gridCellsNumFluids

	static const cl_uint MAX_PLATFORMS = 16;		//Arbitrary value
	static const cl_uint MAX_DEVICES = 16;			//Arbitrary value
	
	cl_platform_id platform_id;
	cl_context context;

	cl_device_id gpu_device;
	cl_command_queue gpu_command_queue;
	
	cl_program fluids_program;
	cl_kernel kernel_grid_insertParticles;
	cl_kernel kernel_sph_computePressure;
	cl_kernel kernel_sph_computeForce;
	cl_kernel kernel_advance;

	cl_mem buffer_gridParams;				//GridParameters
	cl_mem buffer_fluidParams;				//FluidParams_float
	
	cl_mem buffer_gridCells;				//int[]
	cl_mem buffer_gridCellsNumFluids;		//int[]
	cl_mem buffer_fluids;					//Fluid[]
	cl_mem buffer_neighborTables;			//Neighbors[]
	
public:	
	FluidSystem_OpenCL();
	
	void initialize();
	void deactivate();
	
	void stepSimulation(FluidParameters_float *fluidParams, GridParameters *gridParams,
						int numFluidParticles, Fluid *fluids, Neighbors *neighbors,
						int numGridCells, int *gridCells, int *gridCellsNumFluids);
						
private:
	void initialize_stage1_platform();
	void initialize_stage2_device();
	void initialize_stage3_context_and_queue();
	void initialize_stage4_program_and_buffers();

	void deactivate_stage1_program_and_buffers();
	void deactivate_stage2_context_and_queue();
	
	void writeToOpencl(	FluidParameters_float *fluidParams, GridParameters *gridParams,
						int numFluidParticles, Fluid *fluids, Neighbors *neighbors,
						int numGridCells, int *gridCells, int *gridCellsNumFluids	);
	void readFromOpencl(FluidParameters_float *fluidParams, GridParameters *gridParams,
						int numFluidParticles, Fluid *fluids, Neighbors *neighbors,
						int numGridCells, int *gridCells, int *gridCellsNumFluids	);
	
	void grid_insertParticles(int numFluidParticles);
	void sph_computePressure(int numFluidParticles);
	void sph_computeForce(int numFluidParticles);
	void advance(int numFluidParticles);
};

//#define COMPILE_MARCHING_CUBES_OPENCL
#ifdef COMPILE_MARCHING_CUBES_OPENCL

//	draft; incomplete and untested
#include "../vector3df.h"
class MarchingCubes_OpenCL
{
	static const int MAX_MARCHING_CUBES_CELLS = 32768;

	cl_kernel kernel_loadMarchingCubeVertex;
	cl_mem buffer_marching_cubes_cells;		//float[]
	
public:
	MarchingCubes_OpenCL()
	{
		kernel_loadMarchingCubeVertex = INVALID_KERNEL;
		buffer_marching_cubes_cells = INVALID_BUFFER;
	}

	void initialize_marching_cubes(cl_program fluids_program, cl_context context)
	{
		cl_int error_code;
		
		kernel_loadMarchingCubeVertex = clCreateKernel(fluids_program, "loadMarchingCubeVertex", &error_code);
		CHECK_CL_ERROR(error_code);
			
		buffer_marching_cubes_cells = clCreateBuffer(context, CL_MEM_READ_WRITE, sizeof(float)*MAX_MARCHING_CUBES_CELLS, NULL, &error_code);
		CHECK_CL_ERROR(error_code);

	}
	void deactivate_marching_cubes()
	{
		cl_int error_code;
		
		error_code = clReleaseMemObject(buffer_marching_cubes_cells);
		CHECK_CL_ERROR(error_code);
		error_code = clReleaseKernel(kernel_loadMarchingCubeVertex);
		CHECK_CL_ERROR(error_code);
		
		kernel_loadMarchingCubeVertex = INVALID_KERNEL;
		buffer_marching_cubes_cells = INVALID_BUFFER;
	}
	
	void load_marching_cubes_cells(cl_command_queue gpu_command_queue,
								   btVector3 min, btVector3 cell_size,
								   Fluid *fluid,  GridParameters *gridParams, int *gridCells, float gridCellSize,
								   float *marching_cubes_cells, int marching_cubes_cells_per_edge)
	{
		cl_int error_code;
		
		const int NUM_CELLS = marching_cubes_cells_per_edge * marching_cubes_cells_per_edge * marching_cubes_cells_per_edge;
		
		if(NUM_CELLS > MAX_MARCHING_CUBES_CELLS)
		{
			printf("FluidSystem_OpenCL::load_marching_cubes_cells() error: NUM_CELLS > MAX_MARCHING_CUBES_CELLS(%d > %d)\n",
				   NUM_CELLS, MAX_MARCHING_CUBES_CELLS);
				   
			return;
		}
		
		//
		error_code = clEnqueueWriteBuffer(gpu_command_queue, buffer_marching_cubes_cells, CL_TRUE, 0, 
										  sizeof(float)*NUM_CELLS, marching_cubes_cells, 0, NULL, NULL);
		CHECK_CL_ERROR(error_code);
		
		//
		{
			BT_PROFILE("FluidSystem_OpenCL::load_marching_cubes_cells()");
			
			///		const btVector3 min
			///		const btVector3 cell_size
			///		 __global Fluid *fluids
			///		 __global GridParameters *gridParams
			///		 __global int *gridCells
			///		float gridCellSize
			///		__global float *march_cells
			///		int march_cells_per_edge
			error_code = clSetKernelArg( kernel_loadMarchingCubeVertex, 0, sizeof(btVector3), reinterpret_cast<void*>(&min) );
			CHECK_CL_ERROR(error_code);
			error_code = clSetKernelArg( kernel_loadMarchingCubeVertex, 1, sizeof(btVector3), reinterpret_cast<void*>(&cell_size) );
			CHECK_CL_ERROR(error_code);
			error_code = clSetKernelArg( kernel_loadMarchingCubeVertex, 2, sizeof(void*), reinterpret_cast<void*>(&fluid) );
			CHECK_CL_ERROR(error_code);
			error_code = clSetKernelArg( kernel_loadMarchingCubeVertex, 3, sizeof(void*), reinterpret_cast<void*>(&gridParams) );
			CHECK_CL_ERROR(error_code);
			error_code = clSetKernelArg( kernel_loadMarchingCubeVertex, 4, sizeof(int*), reinterpret_cast<void*>(&gridCells) );
			CHECK_CL_ERROR(error_code);
			error_code = clSetKernelArg( kernel_loadMarchingCubeVertex, 5, sizeof(float), reinterpret_cast<void*>(&gridCellSize) );
			CHECK_CL_ERROR(error_code);
			error_code = clSetKernelArg( kernel_loadMarchingCubeVertex, 6, sizeof(float*), reinterpret_cast<void*>(&marching_cubes_cells) );
			CHECK_CL_ERROR(error_code);
			error_code = clSetKernelArg( kernel_loadMarchingCubeVertex, 7, sizeof(int), reinterpret_cast<void*>(&marching_cubes_cells_per_edge) );
			CHECK_CL_ERROR(error_code);

			
			size_t global_work_size[3] = { marching_cubes_cells_per_edge, marching_cubes_cells_per_edge, marching_cubes_cells_per_edge};
			error_code = clEnqueueNDRangeKernel(gpu_command_queue, kernel_loadMarchingCubeVertex, 3, NULL, global_work_size, NULL, 0, NULL, NULL);
			CHECK_CL_ERROR(error_code);
			
			error_code = clFinish(gpu_command_queue);
			CHECK_CL_ERROR(error_code);
		}
		
		//
		error_code = clEnqueueReadBuffer(gpu_command_queue, buffer_marching_cubes_cells, CL_TRUE, 0, 
										 sizeof(float)*NUM_CELLS, marching_cubes_cells, 0, NULL, NULL);
		CHECK_CL_ERROR(error_code);
	}
};
#endif


#endif


