/* FluidSolverOpenCL.h
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

#ifndef FLUID_SOLVER_OPENCL_H
#define FLUID_SOLVER_OPENCL_H

#include <CL/cl.h>

#include "LinearMath/btAlignedObjectArray.h"

#include "opencl_support.h"
#include "../FluidSolver.h"

struct FluidParametersGlobal;
struct FluidParametersLocal;
struct FluidParticles;
class FluidSph;
class FluidStaticGrid;


struct FluidStaticGrid_OpenCLPointers
{
	void *m_buffer_gridParams;
	
	void *m_buffer_gridCells;
	void *m_buffer_gridCellsNumFluids;
};

///@brief Manages OpenCL buffers corresponding to a FluidStaticGrid.
class FluidStaticGrid_OpenCL
{
	OpenCLBuffer m_buffer_gridParams;			//FluidStaticGridParameters
	
	int m_numGridCells;
	OpenCLBuffer m_buffer_gridCells;			//int[]
	OpenCLBuffer m_buffer_gridCellsNumFluids;	//int[]	
	
public:
	FluidStaticGrid_OpenCL() : m_numGridCells(0) {}
	~FluidStaticGrid_OpenCL() { deallocate(); }
	
	void writeToOpenCL(cl_context context, cl_command_queue commandQueue, FluidStaticGrid *grid);
	void readFromOpenCL(cl_context context, cl_command_queue commandQueue, FluidStaticGrid *grid);
	
	FluidStaticGrid_OpenCLPointers getPointers();
	
private:
	void allocate(cl_context context, int numGridCells);
	void deallocate();
};

struct Fluid_OpenCLPointers
{
	void *m_buffer_localParameters;
	
	void *m_buffer_pos;
	void *m_buffer_vel;
	void *m_buffer_vel_eval;
	void *m_buffer_sph_force;
	void *m_buffer_externalAcceleration;
	void *m_buffer_pressure;
	void *m_buffer_density;
	void *m_buffer_nextFluidIndex;
	
	void *m_buffer_neighborTable;
};

///@brief Manages OpenCL buffers corresponding to FluidParticles and a FluidParametersLocal.
class Fluid_OpenCL
{
	OpenCLBuffer m_buffer_localParameters;			//FluidParametersLocal
	
	//
	int m_maxParticles;
	
	OpenCLBuffer m_buffer_pos;						//btVector3[]
	OpenCLBuffer m_buffer_vel;						//btVector3[]
	OpenCLBuffer m_buffer_vel_eval;					//btVector3[]
	OpenCLBuffer m_buffer_sph_force;				//btVector3[]
	OpenCLBuffer m_buffer_externalAcceleration;		//btVector3[]
	OpenCLBuffer m_buffer_pressure;					//float[]
	OpenCLBuffer m_buffer_density;					//float[]
	OpenCLBuffer m_buffer_nextFluidIndex;			//int[]
	
	OpenCLBuffer m_buffer_neighborTable;			//NeighborTable[]

public:
	Fluid_OpenCL() : m_maxParticles(0) {}
	~Fluid_OpenCL() { deallocate(); }
	
	void writeToOpenCL(cl_context context, cl_command_queue commandQueue, 
					   const FluidParametersLocal &FL, FluidParticles *particles, bool transferAllData);
	void readFromOpenCL(cl_context context, cl_command_queue commandQueue,
						FluidParticles *particles, bool transferAllData);
	
	Fluid_OpenCLPointers getPointers();
	
private:
	void allocate(cl_context context, int maxParticles);
	void deallocate();
};


//Loaded into FluidSolverOpenCL.fluids_program
const char CL_PROGRAM_PATH[] = "./Demos/FluidDemo/Fluids/OpenCL_support/fluids.cl";

///@brief Implementation of FluidSolverGridNeighbor for GPU.
///@remarks
///Requires use of FluidGrid::FT_LinkedList; FluidSph created using 
///FluidGrid::FT_IndexRange are excluded from the calculations and not updated.
///@par
///Does not implement fluid-fluid interactions.
class FluidSolverOpenCL : public FluidSolver
{
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
	cl_kernel kernel_integrate;

	OpenCLBuffer buffer_globalFluidParams;		//FluidParametersGlobal
	
	btAlignedObjectArray<Fluid_OpenCL> m_fluidData;
	btAlignedObjectArray<FluidStaticGrid_OpenCL> m_gridData;
	
public:	
	FluidSolverOpenCL();
	~FluidSolverOpenCL() { deactivate(); }
	
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
	
	void grid_insertParticles(int numFluidParticles, FluidStaticGrid_OpenCLPointers *gridPointers, Fluid_OpenCLPointers *fluidPointers);
	void sph_computePressure(int numFluidParticles, FluidStaticGrid_OpenCLPointers *gridPointers, Fluid_OpenCLPointers *fluidPointers);
	void sph_computeForce(int numFluidParticles, FluidStaticGrid_OpenCLPointers *gridPointers, Fluid_OpenCLPointers *fluidPointers);
	void integrate(int numFluidParticles, Fluid_OpenCLPointers *fluidPointers);
};

#endif


