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
#include "btFluidSphSolverPCISPH.h"

#include "btFluidSortingGrid.h"

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

void determineNeighborsInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
											const btFluidSortingGrid& grid, btFluidParticles& particles,
											btFluidSphSolverPCISPH::PciSphParticles& pciSphData)
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
					
					if(distanceSquared < FG.m_sphRadiusSquared)
					{
						btScalar distance = btSqrt(distanceSquared);
						if( !pciSphData.m_neighborTable[i].isFilled() ) pciSphData.m_neighborTable[i].addNeighbor(n, distance);
						else if( !pciSphData.m_neighborTable[n].isFilled() ) pciSphData.m_neighborTable[n].addNeighbor(i, distance);
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

void computeViscosityForceNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, int particleIndex, 
													btFluidParticles& particles, btFluidSphSolverPCISPH::PciSphParticles& pciSphData)
{
	int i = particleIndex;
	
	for(int j = 0; j < pciSphData.m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = pciSphData.m_neighborTable[i].getNeighborIndex(j);
		btScalar distance = pciSphData.m_neighborTable[i].getDistance(j);
		
		if(distance >= FG.m_sphSmoothRadius) continue;
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		
		//Assume that particles are at rest density
		btScalar viscosityScalar = closeness;	//Partial result of Laplacian of viscosity kernel
		//btScalar viscosityScalar = closeness / pciSphData.m_density[n];
		
		btVector3 viscosityForce = (particles.m_vel_eval[n] - particles.m_vel_eval[i]) * viscosityScalar;
		
		pciSphData.m_viscosityForce[i] += viscosityForce;
		pciSphData.m_viscosityForce[n] += -viscosityForce;
	}
}
void computeViscosityForceInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
										const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverPCISPH::PciSphParticles& pciSphData)
{
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computeViscosityForceNeighborTableSymmetric(FG, particleIndex, particles, pciSphData);
	}
}

void computeSumsNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, int particleIndex, 
										btFluidParticles& particles, btFluidSphSolverPCISPH::PciSphParticles& pciSphData,
										const btAlignedObjectArray<btVector3>& particlePositions)
{
	int i = particleIndex;
	
	for(int j = 0; j < pciSphData.m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = pciSphData.m_neighborTable[i].getNeighborIndex(j);
		
		btVector3 difference = (particlePositions[i] - particlePositions[n]) * FG.m_simulationScale;
		
		btScalar distanceSquared = difference.length2();
		pciSphData.m_neighborTable[i].updateDistance( j, btSqrt(distanceSquared) );
		
		if(distanceSquared >= FG.m_sphRadiusSquared) continue;
		btScalar c = FG.m_sphRadiusSquared - distanceSquared;
		
		btScalar poly6KernPartialResult = c * c * c;
		
		pciSphData.m_density[i] += poly6KernPartialResult;
		pciSphData.m_density[n] += poly6KernPartialResult;
	}
}
void computeSumsInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
								const btFluidSortingGrid& grid, btFluidParticles& particles,
								btFluidSphSolverPCISPH::PciSphParticles& pciSphData,
								const btAlignedObjectArray<btVector3>& particlePositions)
{
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computeSumsNeighborTableSymmetric(FG, particleIndex, particles, pciSphData, particlePositions);
	}
}


void computePressureForceNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, int particleIndex, 
													btFluidParticles& particles, btFluidSphSolverPCISPH::PciSphParticles& pciSphData)
{
	int i = particleIndex;
	
	for(int j = 0; j < pciSphData.m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = pciSphData.m_neighborTable[i].getNeighborIndex(j);
		btScalar distance = pciSphData.m_neighborTable[i].getDistance(j);
		if(distance >= FG.m_sphSmoothRadius) continue;
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		
		btScalar pressureScalar = (closeness*closeness) * (pciSphData.m_pressure[i] + pciSphData.m_pressure[n]) / pciSphData.m_density[n];
		
		//btVector3 simScaleNormal = (particles.m_pos[i] - particles.m_pos[n]) * (FG.m_simulationScale / distance);
		btVector3 simScaleNormal = (pciSphData.m_predictedPosition[i] - pciSphData.m_predictedPosition[n]) * (FG.m_simulationScale / distance);
		
		btVector3 pressureForce = simScaleNormal * pressureScalar;
		
		pciSphData.m_pressureForce[i] += pressureForce;
		pciSphData.m_pressureForce[n] += -pressureForce;
	}
}
void computePressureForceInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
										const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverPCISPH::PciSphParticles& pciSphData)
{
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computePressureForceNeighborTableSymmetric(FG, particleIndex, particles, pciSphData);
	}
}

void btFluidSphSolverPCISPH::updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids)
{
	BT_PROFILE("btFluidSphSolverPCISPH::updateGridAndCalculateSphForces()");
	
	//SPH data is discarded/recalculated every frame, so only 1
	//set of arrays are needed if there is no fluid-fluid interaction.
	if( m_pciSphData.size() != 1 ) m_pciSphData.resize(1);
	
	for(int i = 0; i < numFluids; ++i)
	{
		btFluidSph* fluid = fluids[i];
		const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
		const btFluidSortingGrid& grid = fluid->getGrid();
		btFluidParticles& particles = fluid->internalGetParticles();
		btFluidSphSolverPCISPH::PciSphParticles& pciSphData = m_pciSphData[0];
		
		if( fluid->numParticles() > pciSphData.size() ) pciSphData.resize( fluid->numParticles() );
		
		fluid->insertParticlesIntoGrid();
		
		//Generate neighbor tables
		for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_neighborTable[n].clear();
		for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
		{
			const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
		
			for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
				determineNeighborsInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, pciSphData);
		}
		
		//Calculate viscosity force
		//Other forces(gravity, collision) should also be included here
		{
			for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_viscosityForce[n].setValue(0,0,0);
			
			for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
			{
				const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
			
				for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
					computeViscosityForceInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, pciSphData);
			}
			
			//Assume that particles are at rest density
			const btScalar viscosityForceConstants = FG.m_viscosityKernLapCoeff * FL.m_viscosity * FL.m_particleMass / FL.m_restDensity;
			for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_viscosityForce[n] *= viscosityForceConstants;
		}
		
		//Initialize pressure, pressure force
		for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_pressureForce[n].setValue(0,0,0);
		for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_pressure[n] = btScalar(0.0);
		
		const int MIN_ITERATIONS = 4;
		for(int j = 0; j < MIN_ITERATIONS; ++j)
		{
			//Predict velocity
			for(int n = 0; n < fluid->numParticles(); ++n) 
			{
				btVector3 acceleration = (pciSphData.m_viscosityForce[n] + pciSphData.m_pressureForce[n]) / FL.m_restDensity;
			
				pciSphData.m_predictedVelocity[n] = particles.m_vel[n] + acceleration * FG.m_timeStep;	
			}
				
			//Predict position
				//Velocity is at simulation scale; divide by simulation scale to convert to world scale
			btScalar timeStepDivSimScale = FG.m_timeStep / FG.m_simulationScale;
			for(int n = 0; n < fluid->numParticles(); ++n) 
				pciSphData.m_predictedPosition[n] = particles.m_pos[n] + pciSphData.m_predictedVelocity[n] * timeStepDivSimScale;
			
			//Predict density
			{
				for(int n = 0; n < fluid->numParticles(); ++n) 
					pciSphData.m_density[n] = FG.m_initialSum;
				
				for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
				{
					const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
					
					for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
						computeSumsInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, pciSphData, pciSphData.m_predictedPosition);
				}
				
				for(int n = 0; n < fluid->numParticles(); ++n) 
					pciSphData.m_density[n] *= FL.m_sphParticleMass * FG.m_poly6KernCoeff;
			}
			
			//Predict density variation
			for(int n = 0; n < fluid->numParticles(); ++n) 
				pciSphData.m_densityError[n] = pciSphData.m_density[n] - FL.m_restDensity;
			
			const bool PRINT_DENSITY_ERROR = true;
			if( PRINT_DENSITY_ERROR && fluid->numParticles() )
			{
				btScalar minDensityError(BT_LARGE_FLOAT);
				btScalar maxDensityError(-BT_LARGE_FLOAT);
				btScalar totalDensityError(0.0);
				for(int n = 0; n < fluid->numParticles(); ++n)
				{
					minDensityError = btMin(minDensityError, pciSphData.m_densityError[n]);
					maxDensityError = btMax(maxDensityError, pciSphData.m_densityError[n]);
					totalDensityError += pciSphData.m_densityError[n];
				}
				
				printf("iteration: %d \n", j);
				printf("min, max density error: %f, %f \n", minDensityError, maxDensityError);
				printf("total density error: %f, %f \n", totalDensityError);
				printf( "avg density error: %f, %f \n", totalDensityError  / static_cast<btScalar>(fluid->numParticles()) );
				printf("\n");
			}
			
			//Update pressure
			//btScalar beta_PCISPH = (btScalar(2.0)*FG.m_timeStep*FG.m_timeStep*FL.m_particleMass*FL.m_particleMass) / (FL.m_restDensity*FL.m_restDensity);
			//btScalar delta_PCISPH = btScalar(1.0) / beta;
			
			//btScalar beta = FG.m_timeStep*FG.m_timeStep / FL.m_restDensity;	//	should include simulation scale?
			//btScalar precomputed = btScalar(1.0) / beta;
			//btScalar delta = btScalar(-1.0) / (beta * -precomputed);
			
			btScalar delta = btScalar(0.5);
			for(int n = 0; n < fluid->numParticles(); ++n) 
				pciSphData.m_pressure[n] += pciSphData.m_densityError[n] * delta;
			
			//Compute pressure force
			for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_pressureForce[n].setValue(0,0,0);
			for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
			{
				const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
			
				for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
					computePressureForceInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, pciSphData);
			}
			
			const btScalar pressureForceConstants = btScalar(-0.5) * FG.m_spikyKernGradCoeff * FL.m_particleMass;
			for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_pressureForce[n] *= pressureForceConstants;
		}
		
		//Apply SPH force to particles
		for(int n = 0; n < fluid->numParticles(); ++n) 
		{
			btVector3 sphAcceleration = (pciSphData.m_viscosityForce[n] + pciSphData.m_pressureForce[n]) / pciSphData.m_density[n];
			
			fluid->applyForce(n, sphAcceleration * FL.m_particleMass);
		}
	}
}
