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

#include "BulletFluids/Sph/btFluidSortingGrid.h"

#include "LinearMath/btQuickprof.h"		//BT_PROFILE(name) macro

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

#define ASSUME_REST_DENSITY_FOR_VISCOSITY_FORCE
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
		
#ifdef ASSUME_REST_DENSITY_FOR_VISCOSITY_FORCE
		//Assume that particles are at rest density
		btScalar viscosityScalar = closeness;	//Partial result of Laplacian of viscosity kernel
#else
		btScalar viscosityScalar = closeness / pciSphData.m_density[n];
#endif
	
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
													btFluidParticles& particles, btFluidSphSolverPCISPH::PciSphParticles& pciSphData,
													const btAlignedObjectArray<btVector3>& particlePositions)
{
	int i = particleIndex;
	
	for(int j = 0; j < pciSphData.m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = pciSphData.m_neighborTable[i].getNeighborIndex(j);
		btScalar distance = pciSphData.m_neighborTable[i].getDistance(j);
		if(distance >= FG.m_sphSmoothRadius) continue;
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;
	
		btScalar pressureScalar = (closeness*closeness) * (pciSphData.m_pressure[i] + pciSphData.m_pressure[n]) / pciSphData.m_density[n];
		
		//Alternate pressure force
		//btScalar pressureScalar = (closeness*closeness);
		//pressureScalar *= pciSphData.m_pressure[i] / (pciSphData.m_density[i]*pciSphData.m_density[i]); 
		//pressureScalar *= pciSphData.m_pressure[n] / (pciSphData.m_density[n]*pciSphData.m_density[n]);
		
		btVector3 simScaleNormal = (particlePositions[i] - particlePositions[n]) * (FG.m_simulationScale / distance);
		btVector3 pressureForce = simScaleNormal * pressureScalar;
		
		pciSphData.m_pressureForce[i] += pressureForce;
		pciSphData.m_pressureForce[n] += -pressureForce;
	}
}
void computePressureForceInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
										const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverPCISPH::PciSphParticles& pciSphData,
										const btAlignedObjectArray<btVector3>& particlePositions)
{
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computePressureForceNeighborTableSymmetric(FG, particleIndex, particles, pciSphData, particlePositions);
	}
}


inline void resolveAabbCollision(const btFluidSphParametersLocal& FL, const btVector3& vel_eval,
								 btVector3* acceleration, const btVector3& normal, btScalar distance)
{
	if( distance < btScalar(0.0) )	//Negative distance indicates penetration
	{
		btScalar penetrationDepth = -distance;
	
		btScalar accelerationMagnitude = FL.m_boundaryStiff * penetrationDepth - FL.m_boundaryDamp * normal.dot(vel_eval);
		
		*acceleration += normal * accelerationMagnitude;
	}
}
btVector3 determineAabbAcceleration(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, btScalar simScaleParticleRadius,
									btFluidParticles& particles, int particleIndex)
{
	int i = particleIndex;
	
	const btScalar radius = simScaleParticleRadius;
	const btScalar simScale = FG.m_simulationScale;
	
	const btVector3& min = FL.m_aabbBoundaryMin;
	const btVector3& max = FL.m_aabbBoundaryMax;
	
	const btVector3& pos = particles.m_pos[i];
	const btVector3& vel_eval = particles.m_vel_eval[i];
	
	btVector3 acceleration(0,0,0);
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3( 1.0, 0.0, 0.0), ( pos.x() - min.x() )*simScale - radius );
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3(-1.0, 0.0, 0.0), ( max.x() - pos.x() )*simScale - radius );
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3(0.0,  1.0, 0.0), ( pos.y() - min.y() )*simScale - radius );
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3(0.0, -1.0, 0.0), ( max.y() - pos.y() )*simScale - radius );
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3(0.0, 0.0,  1.0), ( pos.z() - min.z() )*simScale - radius );
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3(0.0, 0.0, -1.0), ( max.z() - pos.z() )*simScale - radius );
	
	return acceleration;
}

//Calculates a recommended(but not ideal) stiffness for the fluid.
//This should be called when the fluid is at rest, and stacked high enough(but not too high)
//so that particles in the center are fully surrounded. A particle in the center should
//have other particles around it up to at least FG.m_sphSmoothRadius.
btScalar computeStiffness(const btFluidSphParametersGlobal& FG, const btFluidSph* fluid)
{
	if( !fluid->numParticles() ) return btScalar(1.0);

	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	const btFluidSortingGrid& grid = fluid->getGrid();
	const btFluidParticles& particles = fluid->getParticles();
	
	const btScalar poly6KernGradCoeff = btScalar(-945.0) / ( btScalar(32.0) * SIMD_PI * btPow(FG.m_sphSmoothRadius, 9) );
	
	//int maxNeighbors = 0;
	btScalar kernelResultMax;
	btScalar kernelResultMin;
	for(int i = 0; i < particles.size(); ++i)
	{
		int numNeighbors = 0;
		btVector3 poly6KernelGradientSum(0,0,0);
		btVector3 spikyKernelGradientSum(0,0,0);
		btScalar kernelDotProductSum(0.0);
		
		btFluidSortingGrid::FoundCells foundCells;
		grid.findCells(particles.m_pos[i], foundCells);
		for(int cell = 0; cell < btFluidSortingGrid::NUM_FOUND_CELLS; ++cell) 
		{
			const btFluidGridIterator& FI = foundCells.m_iterators[cell];
			for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
			{
				btVector3 difference = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
				btScalar distanceSquared = difference.length2();
				
				if(distanceSquared < FG.m_sphRadiusSquared) 
				{
					++numNeighbors;
				
					btScalar distance = btSqrt(distanceSquared);
					if(distance < SIMD_EPSILON) distance = SIMD_EPSILON;
					
					btScalar closeness = FG.m_sphSmoothRadius - distance;
					btScalar squaredCloseness = FG.m_sphRadiusSquared - distanceSquared;
					
					btVector3 normal = difference / distance;
					
					btVector3 poly6KernelGradient = difference * (poly6KernGradCoeff * squaredCloseness * squaredCloseness);
					btVector3 spikyKernelGradient = normal * (FG.m_spikyKernGradCoeff * closeness * closeness);
					
					poly6KernelGradientSum += poly6KernelGradient;
					spikyKernelGradientSum += spikyKernelGradient;
					kernelDotProductSum += poly6KernelGradient.dot(spikyKernelGradient);
				}
			}
		}
		
		btScalar particleKernelResult = -spikyKernelGradientSum.dot(poly6KernelGradientSum) - kernelDotProductSum;
		
		//kernelResult should be calculated using a particle with a filled neighborhood
		if(i != 0)
		{
			kernelResultMin = (particleKernelResult < kernelResultMin) ? particleKernelResult : kernelResultMin;
			kernelResultMax = (particleKernelResult > kernelResultMax) ? particleKernelResult : kernelResultMax;
		}
		else
		{
			kernelResultMin = particleKernelResult;
			kernelResultMax = particleKernelResult;
		}
	}
	
	
	btScalar beta = (btScalar(2.0)*FG.m_timeStep*FG.m_timeStep*FL.m_particleMass*FL.m_particleMass) / (FL.m_restDensity*FL.m_restDensity);
	
	btScalar deltaMin = btScalar(-1.0) / (beta * kernelResultMin);
	btScalar deltaMax = btScalar(-1.0) / (beta * kernelResultMax);
	
	printf("beta: %f \n", beta);
	printf("kernelResultMin: %f \n", kernelResultMin);
	printf("kernelResultMax: %f \n", kernelResultMax);
	printf("beta * kernelResultMin: %f \n", beta * kernelResultMin);
	printf("beta * kernelResultMax: %f \n", beta * kernelResultMax);
	printf("deltaMin: %f \n", deltaMin);
	printf("deltaMax: %f \n", deltaMax);
	//printf("maxNeighbors: %d \n", maxNeighbors);
	
	return deltaMin;
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
		if( !fluid->numParticles() ) continue;
		
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
		

#ifndef ASSUME_REST_DENSITY_FOR_VISCOSITY_FORCE
		{
			for(int n = 0; n < fluid->numParticles(); ++n) 
				pciSphData.m_density[n] = FG.m_initialSum;
			
			for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
			{
				const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
				
				for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
					computeSumsInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, pciSphData, particles.m_pos);
			}
			
			for(int n = 0; n < fluid->numParticles(); ++n) 
				pciSphData.m_density[n] *= FL.m_sphParticleMass * FG.m_poly6KernCoeff;
		}
#endif			
		
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
			

#ifdef ASSUME_REST_DENSITY_FOR_VISCOSITY_FORCE
			//Assume that particles are at rest density
			const btScalar viscosityForceConstants = FG.m_viscosityKernLapCoeff * FL.m_viscosity * FL.m_particleMass / FL.m_restDensity;
#else
			const btScalar viscosityForceConstants = FG.m_viscosityKernLapCoeff * FL.m_viscosity * FL.m_particleMass;
#endif
			for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_viscosityForce[n] *= viscosityForceConstants;
		}
		
		//Initialize pressure, pressure force
		for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_pressureForce[n].setValue(0,0,0);
		for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_pressure[n] = btScalar(0.0);
		
		
		btScalar currentDensityError(BT_LARGE_FLOAT);	//Arbitrary initial value, higher than MAX_DENSITY_ERROR
		const btScalar MAX_DENSITY_ERROR(0.0);			//Percentage; set to 0 to keep iterating even if the fluid is below the density threshold
		const int MIN_ITERATIONS = 4;
		const int MAX_ITERATIONS = 4;
		for(int j = 0; (j < MIN_ITERATIONS || currentDensityError > MAX_DENSITY_ERROR) && j < MAX_ITERATIONS; ++j)
		{
			//Predict velocity
			for(int n = 0; n < fluid->numParticles(); ++n) 
			{
#ifdef ASSUME_REST_DENSITY_FOR_VISCOSITY_FORCE
				btScalar density = (j == 0) ? FL.m_restDensity : pciSphData.m_density[n];
				btVector3 acceleration = (pciSphData.m_viscosityForce[n] + pciSphData.m_pressureForce[n]) / density;
#else
				//Alternate pressure force
				//btVector3 viscosityAccel = (pciSphData.m_viscosityForce[n]) / pciSphData.m_density[n];
				//btVector3 pressureAccel = pciSphData.m_pressureForce[n];
				//btVector3 acceleration = viscosityAccel + pressureAccel;

				btVector3 acceleration = (pciSphData.m_viscosityForce[n] + pciSphData.m_pressureForce[n]) / pciSphData.m_density[n];
#endif

				const bool CONSIDER_EXTERNAL_FORCES = false;
				if(CONSIDER_EXTERNAL_FORCES)
				{
					acceleration += FL.m_gravity;
				
					const btScalar simScaleParticleRadius = FL.m_particleRadius * FG.m_simulationScale;
					acceleration += determineAabbAcceleration(FG, FL, simScaleParticleRadius, particles, i);
				}
				
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
			btScalar maxDensityError(0.0);
			for(int n = 0; n < fluid->numParticles(); ++n)
			{
				btScalar densityError = pciSphData.m_density[n] - FL.m_restDensity;
			
				maxDensityError = btMax(maxDensityError, densityError);
				pciSphData.m_densityError[n] = densityError;
			}
			currentDensityError = (maxDensityError / FL.m_restDensity) * btScalar(100.0);
			
			const bool PRINT_DENSITY_ERROR = true;
			if( PRINT_DENSITY_ERROR && fluid->numParticles() )
			{
				btScalar totalPositiveDensityError(0.0);
				for(int n = 0; n < fluid->numParticles(); ++n) totalPositiveDensityError += btMax( btScalar(0.0), pciSphData.m_densityError[n]);
				
				printf("iteration: %d \n", j);
				printf("currentDensityError: %f%% (maxDensityError: %f) \n", currentDensityError, maxDensityError);
				printf("total positive density error: %f \n", totalPositiveDensityError);
				printf("average positive density error: %f \n", totalPositiveDensityError  / static_cast<btScalar>(fluid->numParticles()));
				
				//computeStiffness(FG, fluid);
				
				printf("\n");
			}
			
			
			//Update pressure
			btScalar stiffness = btScalar(0.7);
			//btScalar stiffness = btScalar(0.035281);
			for(int n = 0; n < fluid->numParticles(); ++n) 
				pciSphData.m_pressure[n] += pciSphData.m_densityError[n] * stiffness;
			
			//Compute pressure force
			for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_pressureForce[n].setValue(0,0,0);
			for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
			{
				const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
			
				for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
					computePressureForceInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, pciSphData, pciSphData.m_predictedPosition);
			}
			
			const btScalar pressureForceConstants = btScalar(-0.5) * FG.m_spikyKernGradCoeff * FL.m_particleMass;
			
			//Alternate pressure force
			//const btScalar pressureForceConstants = -FG.m_spikyKernGradCoeff * FL.m_particleMass;
			
			for(int n = 0; n < fluid->numParticles(); ++n) pciSphData.m_pressureForce[n] *= pressureForceConstants;
		}
		
		//Apply SPH force to particles
		for(int n = 0; n < fluid->numParticles(); ++n) 
		{
			btVector3 sphAcceleration = (pciSphData.m_viscosityForce[n] + pciSphData.m_pressureForce[n]) / pciSphData.m_density[n];
		
			//Alternate pressure force
			//btVector3 viscosityAccel = (pciSphData.m_viscosityForce[n]) / pciSphData.m_density[n];
			//btVector3 pressureAccel = pciSphData.m_pressureForce[n];
			//btVector3 sphAcceleration = viscosityAccel + pressureAccel;
			
			fluid->applyForce(n, sphAcceleration * FL.m_particleMass);
		}
	}
}
