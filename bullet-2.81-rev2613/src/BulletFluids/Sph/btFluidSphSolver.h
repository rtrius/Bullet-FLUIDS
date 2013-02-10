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
//Portions of this file based on FLUIDS v.2 - SPH Fluid Simulator for CPU and GPU
//Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com
#ifndef BT_FLUID_SPH_SOLVER_H
#define BT_FLUID_SPH_SOLVER_H

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

#include "btFluidSph.h"

///@brief Interface for particle motion computation. 
///@remarks
///Determines how the positions and velocities of fluid particles change from 
///one simulation step to the next.
class btFluidSphSolver
{
public:
	virtual void updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids) = 0;
	
	static void applyForcesSingleFluid(const btFluidSphParametersGlobal& FG, btFluidSph* fluid);
	static void integratePositionsSingleFluid(const btFluidSphParametersGlobal& FG, btFluidParticles& particles);
	
	
protected:
	static void applySphForce(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, const btAlignedObjectArray<btVector3>& sphForce)
	{
		BT_PROFILE("applySphForce()");
		
		const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
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

///@brief Stores a table of particle indicies and their distances from a single fluid particle.
///@remarks
///Only particles within the SPH interaction radius are included.
///@par
///The table is generated during the pressure calculation step in order to avoid recalculating
///neighboring particles and their distances during force computation.
class btFluidSphNeighbors
{
public:
	static const int MAX_NEIGHBORS_HALVED = 80;
	
private:
	int m_count;
	int m_particleIndicies[MAX_NEIGHBORS_HALVED];
	btScalar m_distances[MAX_NEIGHBORS_HALVED];
	
public:
	inline int numNeighbors() const { return m_count; }
	inline int getNeighborIndex(int index) const { return m_particleIndicies[index]; }
	inline btScalar getDistance(int index) const { return m_distances[index]; }
	inline bool isFilled() const { return (m_count >= MAX_NEIGHBORS_HALVED); }
	
	inline void clear() { m_count = 0; }
	inline void addNeighbor(int neighborIndex, btScalar distance)
	{
		m_particleIndicies[m_count] = neighborIndex;
		m_distances[m_count] = distance;
		++m_count;
	}
	
	inline void updateDistance(int index, btScalar distance) { m_distances[index] = distance; }
};

///@brief Standard CPU fluid solver; solves the Navier-Stokes equations using SPH(Smoothed Particle Hydrodynamics).
///@remarks
///Pressure is calculated using a btFluidSortingGrid, and force using btFluidSphNeighbors 
///table generated during the pressure calculation. Symmetry is exploited by checking
///only 14 of 27 surrounding grid cells, halving the number of calculations.
///@par
///In order to maximize performance, this solver only implements basic SPH.
///Fluid-fluid interaction and surface tension are not implemented.
///@par
///A short introduction to SPH fluids may be found in: \n
///"Particle-Based Fluid Simulation for Interactive Applications". \n
///M. Muller, D. Charypar, M. Gross. Proceedings of 2003 ACM SIGGRAPH Symposium on Computer Animations, p.154-159, 2003. \n
///\n
///For a more extensive overview, see: \n
///"Lagrangian Fluid Dynamics Using Smoothed Particle Hydrodynamics". \n
///M. Kelager. Master's thesis, University of Copenhagen, Department of Computer Science. January 2006. \n
class btFluidSphSolverDefault : public btFluidSphSolver
{
public:
	///Contains parallel arrays that 'extend' btFluidParticles with SPH specific data
	struct SphParticles
	{
		btAlignedObjectArray<btVector3> m_sphForce;		///<Sum of pressure and viscosity forces; simulation scale.
		btAlignedObjectArray<btScalar> m_pressure;		///<Value of the pressure scalar field at the particle's position.
		btAlignedObjectArray<btScalar> m_invDensity;	///<Inverted value of the density scalar field at the particle's position.
		
		btAlignedObjectArray<btFluidSphNeighbors> m_neighborTable;
		
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
	btAlignedObjectArray<btFluidSphSolverDefault::SphParticles> m_sphData;

public:
	virtual void updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids)
	{
		BT_PROFILE("btFluidSphSolverDefault::updateGridAndCalculateSphForces()");
		
		//SPH data is discarded/recalculated every frame, so only 1
		//set of arrays are needed if there is no fluid-fluid interaction.
		if( m_sphData.size() != 1 )m_sphData.resize(1);
		
		for(int i = 0; i < numFluids; ++i) 
		{
			btFluidSph* fluid = fluids[i];
			btFluidSphSolverDefault::SphParticles& sphData = m_sphData[0];
			if( fluid->numParticles() > sphData.size() ) sphData.resize( fluid->numParticles() );
			
			fluid->insertParticlesIntoGrid();
			
			sphComputePressure(FG, fluid, sphData);
			
			sphComputeForce(FG, fluid, sphData);
			
			applySphForce(FG, fluid, sphData.m_sphForce);
		}
	}
	
protected:
	virtual void sphComputePressure(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, btFluidSphSolverDefault::SphParticles& sphData);
	virtual void sphComputeForce(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, btFluidSphSolverDefault::SphParticles& sphData);
	
	//Necessary to use separate classes for multithreaded and single threaded solver
	virtual void computeSumsInMultithreadingGroup(const btFluidSphParametersGlobal& FG, const btAlignedObjectArray<int>& multithreadingGroup,
												const btFluidSortingGrid& grid, btFluidParticles& particles, 
												btFluidSphSolverDefault::SphParticles& sphData)
	{	
		for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
			calculateSumsInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, sphData);
	}
	virtual void computeForcesInMultithreadingGroup(const btFluidSphParametersGlobal& FG, const btScalar vterm, 
													const btAlignedObjectArray<int>& multithreadingGroup, const btFluidSortingGrid& grid, 
													btFluidParticles& particles, btFluidSphSolverDefault::SphParticles& sphData)
	{
		for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
			calculateForcesInCellSymmetric(FG, vterm, multithreadingGroup[cell], grid, particles, sphData);
	}
	
public:
	static void calculateSumsInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, const btFluidSortingGrid& grid, 
											btFluidParticles& particles, btFluidSphSolverDefault::SphParticles& sphData);
	static void calculateForcesInCellSymmetric(const btFluidSphParametersGlobal& FG, const btScalar vterm,
											int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
											btFluidSphSolverDefault::SphParticles& sphData);
};

#endif


