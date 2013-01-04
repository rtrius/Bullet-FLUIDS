/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_FLUID_SPH_SOLVER_OPENCL_H
#define BT_FLUID_SPH_SOLVER_OPENCL_H

#include "LinearMath/btAlignedObjectArray.h"

#include "../../BulletFluids/Sph/btFluidSphSolver.h"

#include "btFluidSphOpenCL.h"
#include "btFluidSortingGridOpenCL.h"

struct btFluidSphParametersGlobal;
class btFluidSph;

///@brief Solver that uses the GPU to accelerate SPH force calculation.
///@remarks
///Does not implement fluid-fluid interactions.
class btFluidSphSolverOpenCL : public btFluidSphSolver
{
	cl_context m_context;
	cl_command_queue m_commandQueue;
	
	cl_program m_fluidsProgram;
	cl_kernel m_kernel_sphComputePressure;
	cl_kernel m_kernel_sphComputeForce;

	btOpenCLArray<btFluidSphParametersGlobal> m_globalFluidParams;
	
	btAlignedObjectArray<btFluidSphOpenCL*> m_fluidData;
	btAlignedObjectArray<btFluidSortingGridOpenCL*> m_gridData;
	
	btFluidSortingGridOpenCLProgram m_sortingGridProgram;
	
	btAlignedObjectArray<btVector3> m_tempSphForce;
	
public:	
	btFluidSphSolverOpenCL(cl_context context, cl_command_queue queue, cl_device_id device);
	virtual ~btFluidSphSolverOpenCL();
	
	virtual void updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids);
	
private:
	void sphComputePressure(int numFluidParticles, btFluidSortingGridOpenCL* gridData, btFluidSphOpenCL* fluidData, btScalar cellSize);
	void sphComputeForce(int numFluidParticles, btFluidSortingGridOpenCL* gridData, btFluidSphOpenCL* fluidData, btScalar cellSize);
};

#endif


