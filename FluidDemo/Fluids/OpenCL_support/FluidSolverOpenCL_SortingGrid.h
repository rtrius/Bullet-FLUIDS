/* FluidSolverOpenCL_SortingGrid.h
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

#ifndef FLUID_SOLVER_OPENCL_SORTING_GRID_H
#define FLUID_SOLVER_OPENCL_SORTING_GRID_H

#include <CL/cl.h>

#include "LinearMath/btAlignedObjectArray.h"

#include "../FluidSolver.h"

#include "opencl_support.h"
#include "FluidOpenCL.h"
#include "FluidSortingGridOpenCL.h"

struct FluidParametersGlobal;
class FluidSph;

//Loaded into FluidSolverOpenCL_SortingGrid.fluids_program
//const char CL_PROGRAM_PATH[] = "./Demos/FluidDemo/Fluids/OpenCL_support/fluids_sortingGrid.cl";

///@brief Implementation of FluidSolverGridNeighbor for GPU.
///@remarks
///Requires use of FluidGrid::FT_IndexRange FluidSph created using 
///FluidGrid::FT_LinkedList are excluded from the calculations and not updated.
///@par
///Does not implement fluid-fluid interactions.
///@par
///Large world support not implemented; SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED must not be defined.
class FluidSolverOpenCL_SortingGrid : public FluidSolver
{
	static const cl_uint MAX_PLATFORMS = 16;		//Arbitrary value
	static const cl_uint MAX_DEVICES = 16;			//Arbitrary value
	
	cl_platform_id platform_id;
	cl_context context;

	cl_device_id gpu_device;
	cl_command_queue gpu_command_queue;
	
	cl_program fluids_program;
	cl_kernel kernel_sph_computePressure;
	cl_kernel kernel_sph_computeForce;
	cl_kernel kernel_integrate;

	OpenCLBuffer buffer_globalFluidParams;		//FluidParametersGlobal
	
	btAlignedObjectArray<Fluid_OpenCL> m_fluidData;
	btAlignedObjectArray<FluidSortingGrid_OpenCL> m_gridData;
	
public:	
	FluidSolverOpenCL_SortingGrid();
	~FluidSolverOpenCL_SortingGrid() { deactivate(); }
	
	virtual void stepSimulation(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids);
	
private:
	void initialize();
	void deactivate();

	void initialize_stage1_platform();
	void initialize_stage2_device();
	void initialize_stage3_context_and_queue();
	void initialize_stage4_program_and_buffer();

	void deactivate_stage1_program_and_buffer();
	void deactivate_stage2_context_and_queue();
	
	void writeSingleFluidToOpenCL(FluidSph* fluid, Fluid_OpenCL *fluidData, FluidSortingGrid_OpenCL *gridData);
	void readSingleFluidFromOpenCL(FluidSph* fluid, Fluid_OpenCL *fluidData, FluidSortingGrid_OpenCL *gridData);
	
	//void grid_insertParticles(int numFluidParticles, FluidSortingGrid_OpenCLPointers *gridPointers, Fluid_OpenCLPointers *fluidPointers);
	void sph_computePressure(int numFluidParticles, FluidSortingGrid_OpenCLPointers *gridPointers, 
							 Fluid_OpenCLPointers *fluidPointers, btScalar cellSize);
	void sph_computeForce(int numFluidParticles, FluidSortingGrid_OpenCLPointers *gridPointers, Fluid_OpenCLPointers *fluidPointers);
	void integrate(int numFluidParticles, Fluid_OpenCLPointers *fluidPointers);
};

#endif


