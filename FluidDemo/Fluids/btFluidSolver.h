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
	virtual void stepSimulation(const btFluidParametersGlobal& FG, btFluidSph** fluids, int numFluids) = 0;
	
protected:
	virtual void integrate(const btFluidParametersGlobal& FG, const btFluidParametersLocal& FL, btFluidParticles& particles);
	
	static void applySphForce(const btFluidParametersGlobal& FG, btFluidSph* fluid, const btAlignedObjectArray<btVector3>& sphForce)
	{
		BT_PROFILE("applySphForce()");
		
		const btFluidParametersLocal& FL = fluid->getLocalParameters();
		btScalar speedLimitSquared = FG.m_speedLimit*FG.m_speedLimit;
		for(int n = 0; n < fluid->numParticles(); ++n) 
		{
			btVector3 acceleration = sphForce[n];
					
			//Limit speed
			btScalar speedSquared = acceleration.length2();
			if(speedSquared > speedLimitSquared) acceleration *= FG.m_speedLimit / btSqrt(speedSquared);
					
			fluid->applyForce(n, acceleration * FL.m_particleMass);
		}
	}
};

///@brief Standard CPU fluid solver; solves the Navier-Stokes equations using SPH(Smoothed Particle Hydrodynamics).
///@remarks
///Pressure is calculated using a btFluidSortingGrid, and force using btFluidNeighbors 
///table generated during the pressure calculation. Symmetry is exploited by checking
///only 14 of 27 surrounding grid cells, halving the number of calculations.
///@par
///Experimental multithreading support is implemented for this solver.
///In testing, performance decreases when over 3 threads are used.
///@par
///In order to maximize performance, this solver only implements basic SPH.
///Fluid-fluid interaction and surface tension are not implemented.
class btFluidSolverSph : public btFluidSolver
{
public:
	///Contains parallel arrays that 'extend' btFluidParticles with SPH specific data
	struct SphParticles
	{
		btAlignedObjectArray<btVector3> m_sphForce;		///<Sum of pressure and viscosity forces; simulation scale.
		btAlignedObjectArray<btScalar> m_pressure;		///<Value of the pressure scalar field at the particle's position.
		btAlignedObjectArray<btScalar> m_invDensity;	///<Inverted value of the density scalar field at the particle's position.
		
		btAlignedObjectArray<btFluidNeighbors> m_neighborTable;
		
		int size() const { return m_sphForce.size(); }
		void resize(int newSize)
		{
			m_sphForce.resize(newSize);
			m_pressure.resize(newSize);
			m_invDensity.resize(newSize);
			
			m_neighborTable.resize(newSize);
		}
	};

protected:
	btAlignedObjectArray<btFluidSolverSph::SphParticles> m_sphData;

public:
	virtual void stepSimulation(const btFluidParametersGlobal& FG, btFluidSph** fluids, int numFluids)
	{
		BT_PROFILE("btFluidSolverSph::stepSimulation()");
		
		//SPH data is discarded/recalculated every frame, so only 1
		//set of arrays are needed if there is no fluid-fluid interaction.
		if( m_sphData.size() != 1 )m_sphData.resize(1);
		
		for(int i = 0; i < numFluids; ++i) 
		{
			btFluidSph* fluid = fluids[i];
			btFluidSolverSph::SphParticles& sphData = m_sphData[0];
			if( fluid->numParticles() > sphData.size() ) sphData.resize( fluid->numParticles() );
			
			fluid->insertParticlesIntoGrid();
			
			sphComputePressure(FG, fluid, sphData);
			
			sphComputeForce(FG, fluid, sphData);
			
			applySphForce(FG, fluid, sphData.m_sphForce);
			
			integrate( FG, fluids[i]->getLocalParameters(), fluid->internalGetParticles() );
		}
	}
	
protected:
	virtual void sphComputePressure(const btFluidParametersGlobal& FG, btFluidSph* fluid, btFluidSolverSph::SphParticles& sphData);
	virtual void sphComputeForce(const btFluidParametersGlobal& FG, btFluidSph* fluid, btFluidSolverSph::SphParticles& sphData);
};

#endif


