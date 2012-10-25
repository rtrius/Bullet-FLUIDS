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

#include "../FluidSolver.h"

#include "opencl_support.h"
#include "FluidOpenCL.h"
#include "FluidSortingGridOpenCL.h"

struct FluidParametersGlobal;
class FluidSph;

///@brief Solver that uses the GPU to accelerate SPH force calculation.
///@remarks
///Does not implement fluid-fluid interactions.
class FluidSolverOpenCL : public FluidSolver
{
	OpenCLConfig m_configCL;
	
	cl_program m_fluidsProgram;
	cl_kernel m_kernel_sphComputePressure;
	cl_kernel m_kernel_sphComputeForce;

	OpenCLBuffer m_buffer_globalFluidParams;		//FluidParametersGlobal
	
	btAlignedObjectArray<Fluid_OpenCL> m_fluidData;
	btAlignedObjectArray<FluidSortingGrid_OpenCL> m_gridData;
	
	FluidSortingGrid_OpenCL_Program m_sortingGridProgram;
	
public:	
	FluidSolverOpenCL();
	~FluidSolverOpenCL() { deactivate(); }
	
	virtual void stepSimulation(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids);
	
private:
	void initialize();
	void deactivate();
	
	void sphComputePressure(int numFluidParticles, FluidSortingGrid_OpenCL *gridData, 
							 Fluid_OpenCL *fluidData, btScalar cellSize);
	void sphComputeForce(int numFluidParticles, FluidSortingGrid_OpenCL *gridData, Fluid_OpenCL *fluidData);
};

#endif


