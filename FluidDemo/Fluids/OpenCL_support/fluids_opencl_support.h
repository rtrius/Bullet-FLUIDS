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

#include <cstdio>

#include <CL/cl.h>

#include "opencl_support.h"

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro


struct FluidParameters;
struct GridParameters;
struct Fluids;
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

	OpenCLBuffer buffer_gridParams;				//GridParameters
	OpenCLBuffer buffer_fluidParams;			//FluidParameters
	
	OpenCLBuffer buffer_gridCells;				//int[]
	OpenCLBuffer buffer_gridCellsNumFluids;		//int[]
	
	OpenCLBuffer buffer_pos;						//btVector3[]
	OpenCLBuffer buffer_vel;						//btVector3[]
	OpenCLBuffer buffer_vel_eval;					//btVector3[]
	OpenCLBuffer buffer_sph_force;					//btVector3[]
	OpenCLBuffer buffer_externalAcceleration;		//btVector3[]
	OpenCLBuffer buffer_prev_pos;					//btVector3[]
	OpenCLBuffer buffer_pressure;					//btScalar[]
	OpenCLBuffer buffer_density;					//btScalar[]
	OpenCLBuffer buffer_nextFluidIndex;				//int[]
	OpenCLBuffer buffer_neighborTables;				//Neighbors[]
	
public:	
	FluidSystem_OpenCL();
	
	void initialize();
	void deactivate();
	
	void stepSimulation(FluidParameters *fluidParams, GridParameters *gridParams, 
						Fluids *fluids, int *gridCells, int *gridCellsNumFluids,
						bool transferAllData);
						
private:
	void initialize_stage1_platform();
	void initialize_stage2_device();
	void initialize_stage3_context_and_queue();
	void initialize_stage4_program_and_buffers();

	void deactivate_stage1_program_and_buffers();
	void deactivate_stage2_context_and_queue();
	
	void writeToOpencl(	FluidParameters *fluidParams, GridParameters *gridParams,
						Fluids *fluids, int *gridCells, int *gridCellsNumFluids,
						int numFluidParticles, int numGridCells );
	void readFromOpencl(FluidParameters *fluidParams, GridParameters *gridParams,
						Fluids *fluids, int *gridCells, int *gridCellsNumFluids,
						int numFluidParticles, int numGridCells );
	
	void grid_insertParticles(int numFluidParticles);
	void sph_computePressure(int numFluidParticles);
	void sph_computeForce(int numFluidParticles);
	void advance(int numFluidParticles);
};


#endif


