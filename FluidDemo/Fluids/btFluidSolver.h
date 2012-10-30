/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. 
   If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
//Portions of this file based on FLUIDS v.2 - SPH Fluid Simulator for CPU and GPU
//Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com
#ifndef BT_FLUID_SOLVER_H
#define BT_FLUID_SOLVER_H

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

#include "btFluidSph.h"

///@brief Interface for particle motion computation. 
///@remarks
///Determines how the positions and velocities of fluid particles change from 
///one simulation step to the next.
class btFluidSolver
{
public:
	virtual void stepSimulation(const btFluidParametersGlobal& FG, btAlignedObjectArray<btFluidSph*>* fluids) = 0;
	
protected:
	virtual void integrate(const btFluidParametersGlobal& FG, const btFluidParametersLocal& FL, btFluidParticles* fluids);
};

///@brief Standard CPU fluid solver; solves the incompressible Navier-Stokes equations using SPH(Smoothed Particle Hydrodynamics).
///@remarks
///Pressure is calculated using a btFluidSortingGrid, and force using btFluidNeighbors 
///table generated during the pressure calculation. Symmetry is exploited by checking
///only 14 of 27 surrounding grid cells, halving the number of calculations.
///@par
///Experimental multithreading support is implemented for this solver.
///In testing, performance decreases when over 3 threads are used.
///@par
///Does not implement fluid-fluid interactions.
class btFluidSolverSph : public btFluidSolver
{
public:
	virtual void stepSimulation(const btFluidParametersGlobal& FG, btAlignedObjectArray<btFluidSph*>* fluids)
	{
		BT_PROFILE("btFluidSolverSph::stepSimulation()");
	
		//
		for(int i = 0; i < fluids->size(); ++i) (*fluids)[i]->insertParticlesIntoGrid();
		
		//
		for(int i = 0; i < fluids->size(); ++i) sphComputePressure( FG, (*fluids)[i] );
			
		for(int i = 0; i < fluids->size(); ++i) sphComputeForce( FG, (*fluids)[i] );
			
		for(int i = 0; i < fluids->size(); ++i)
			integrate( FG, (*fluids)[i]->getLocalParameters(), &(*fluids)[i]->internalGetParticles() );
	}
	
protected:
	virtual void sphComputePressure(const btFluidParametersGlobal& FG, btFluidSph* fluid);
	virtual void sphComputeForce(const btFluidParametersGlobal& FG, btFluidSph* fluid);
};

#endif


