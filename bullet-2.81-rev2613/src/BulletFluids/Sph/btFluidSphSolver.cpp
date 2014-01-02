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
#include "btFluidSphSolver.h"

#include "btFluidSortingGrid.h"

void btFluidSphSolver::applyForcesSingleFluid(const btFluidSphParametersGlobal& FG, btFluidSph* fluid)
{
	BT_PROFILE("btFluidSphSolver::applyForcesSingleFluid()");

	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	btFluidParticles& particles = fluid->internalGetParticles();
	
	const btScalar invParticleMass = btScalar(1.0) / FL.m_particleMass;
	
	for(int i = 0; i < particles.size(); ++i)
	{
		btVector3& vel = particles.m_vel[i];
		btVector3& vel_eval = particles.m_vel_eval[i];
	
		btVector3 acceleration = FL.m_gravity + (particles.m_accumulatedForce[i] * invParticleMass);

		//Leapfrog integration
		btVector3 vnext = vel + acceleration * FG.m_timeStep;	//v(t+1/2) = v(t-1/2) + a(t) dt	
		vel_eval = (vel + vnext) * btScalar(0.5);				//v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute (sph)forces later
		vel = vnext;
	}
	
	for(int i = 0; i < particles.size(); ++i) particles.m_accumulatedForce[i].setValue(0, 0, 0);
}
void btFluidSphSolver::integratePositionsSingleFluid(const btFluidSphParametersGlobal& FG, btFluidParticles& particles)
{
	BT_PROFILE("btFluidSphSolver::integratePositionsSingleFluid()");
	
	//Velocity is at simulation scale; divide by simulation scale to convert to world scale
	btScalar timeStepDivSimScale = FG.m_timeStep / FG.m_simulationScale;
	
	//Leapfrog integration
	//p(t+1) = p(t) + v(t+1/2)*dt
	for(int i = 0; i < particles.size(); ++i) particles.m_pos[i] += particles.m_vel[i] * timeStepDivSimScale;
}

void btFluidSphSolverDefault::sphComputePressure(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, btFluidSphSolverDefault::SphParticles& sphData)
{
	BT_PROFILE("btFluidSphSolverDefault::sphComputePressure()");
	
	const int numParticles = fluid->numParticles();
	
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	const btFluidSortingGrid& grid = fluid->getGrid();
	btFluidParticles& particles = fluid->internalGetParticles();
	
	{
		BT_PROFILE("sphComputePressure() - reset sums, clear table");
		
		const btScalar poly6ZeroDistance = FG.m_sphRadiusSquared * FG.m_sphRadiusSquared * FG.m_sphRadiusSquared;
		const btScalar initialSphSum = poly6ZeroDistance * FL.m_initialSum;
		for(int i = 0; i < numParticles; ++i) sphData.m_invDensity[i] = initialSphSum;
		for(int i = 0; i < numParticles; ++i) sphData.m_neighborTable[i].clear();
	}
	
	{
		BT_PROFILE("sphComputePressure() - compute sums");
		
		for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
		{
			const btAlignedObjectArray<int>& currentGroup = grid.internalGetMultithreadingGroup(group);
			if( !currentGroup.size() ) continue;
			
			computeSumsInMultithreadingGroup(FG, currentGroup, grid, particles, sphData);
		}
	}
	
	{
		BT_PROFILE("sphComputePressure() - compute pressure/density");
		
		for(int i = 0; i < numParticles; ++i)
		{
			btScalar density = sphData.m_invDensity[i] * FL.m_sphParticleMass * FG.m_poly6KernCoeff;
			sphData.m_pressure[i] = (density - FL.m_restDensity) * FL.m_stiffness;
			sphData.m_invDensity[i] = btScalar(1.0) / density;
		}
	}
}

void btFluidSphSolverDefault::sphComputeForce(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, btFluidSphSolverDefault::SphParticles& sphData)
{
	BT_PROFILE("btFluidSphSolverDefault::sphComputeForce()");
	
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	const btFluidSortingGrid& grid = fluid->getGrid();
	btFluidParticles& particles = fluid->internalGetParticles();
	
	btScalar vterm = FG.m_viscosityKernLapCoeff * FL.m_viscosity;
	
	for(int i = 0; i < particles.size(); ++i)sphData.m_sphForce[i].setValue(0, 0, 0);
	
	for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
	{
		const btAlignedObjectArray<int>& currentGroup = grid.internalGetMultithreadingGroup(group);
		if( !currentGroup.size() ) continue;
		
		computeForcesInMultithreadingGroup(FG, vterm, currentGroup, grid, particles, sphData);
	}
	
	for(int i = 0; i < particles.size(); ++i)sphData.m_sphForce[i] *= FL.m_sphParticleMass;
}

void btFluidSphSolverDefault::calculateSumsInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
															const btFluidSortingGrid& grid, btFluidParticles& particles,
															btFluidSphSolverDefault::SphParticles& sphData)
{
	
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	if(currentCell.m_firstIndex <= currentCell.m_lastIndex)	//if cell is not empty
	{
		btFluidSortingGrid::FoundCells foundCells;
		grid.findCellsSymmetric(particles.m_pos[currentCell.m_firstIndex], foundCells);
		
		for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
		{
			//Remove particle, with index i, from grid cell to prevent self-particle interactions
			++foundCells.m_iterators[0].m_firstIndex;	//Local cell; currentCell == foundCells.m_iterators[0]
			
			for(int cell = 0; cell < btFluidSortingGrid::NUM_FOUND_CELLS_SYMMETRIC; cell++) 
			{
				btFluidGridIterator& FI = foundCells.m_iterators[cell];
				
				for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
				{
					//Simulation-scale distance
					btVector3 difference = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		
					btScalar distanceSquared = difference.length2();
					
					if(FG.m_sphRadiusSquared > distanceSquared)
					{
						btScalar c = FG.m_sphRadiusSquared - distanceSquared;
						btScalar poly6KernPartialResult = c * c * c;
						sphData.m_invDensity[i] += poly6KernPartialResult;
						sphData.m_invDensity[n] += poly6KernPartialResult;
						
						btScalar distance = btSqrt(distanceSquared);
						if( !sphData.m_neighborTable[i].isFilled() ) sphData.m_neighborTable[i].addNeighbor(n, distance);
						else if( !sphData.m_neighborTable[n].isFilled() ) sphData.m_neighborTable[n].addNeighbor(i, distance);
						else 
						{
							cell = btFluidSortingGrid::NUM_FOUND_CELLS_SYMMETRIC;	//Break out of outer loop 
							break;
						}
					}
				}
			}
		}
	}
}

void computeForceNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, const btScalar vterm, int particleIndex, 
										btFluidParticles& particles, btFluidSphSolverDefault::SphParticles& sphData)
{
	int i = particleIndex;
	
	for(int j = 0; j < sphData.m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = sphData.m_neighborTable[i].getNeighborIndex(j);
		
		btVector3 difference = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
		btScalar distance = sphData.m_neighborTable[i].getDistance(j);
		
		btScalar c = FG.m_sphSmoothRadius - distance;
		btScalar pterm = btScalar(-0.5) * c * FG.m_spikyKernGradCoeff * (sphData.m_pressure[i] + sphData.m_pressure[n]);
		pterm /= (distance < SIMD_EPSILON) ? SIMD_EPSILON : distance;
		
		btScalar dterm = c * sphData.m_invDensity[i] * sphData.m_invDensity[n];

		btVector3 force(  (pterm * difference.x() + vterm * (particles.m_vel_eval[n].x() - particles.m_vel_eval[i].x())) * dterm,
						  (pterm * difference.y() + vterm * (particles.m_vel_eval[n].y() - particles.m_vel_eval[i].y())) * dterm,
						  (pterm * difference.z() + vterm * (particles.m_vel_eval[n].z() - particles.m_vel_eval[i].z())) * dterm );
		
		sphData.m_sphForce[i] += force;
		sphData.m_sphForce[n] += -force;
	}
}
void btFluidSphSolverDefault::calculateForcesInCellSymmetric(const btFluidSphParametersGlobal& FG, const btScalar vterm,
															int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
															btFluidSphSolverDefault::SphParticles& sphData)
{
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
	{
		computeForceNeighborTableSymmetric(FG, vterm, i, particles, sphData);
	}
}
