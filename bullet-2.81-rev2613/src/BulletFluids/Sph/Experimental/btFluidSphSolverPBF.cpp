/*
Bullet-FLUIDS 
Copyright (c) 2012-2013 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#include "btFluidSphSolverPBF.h"

#include "BulletFluids/Sph/btFluidSortingGrid.h"

#include "LinearMath/btQuickprof.h"		//BT_PROFILE(name) macro


void computeSphSumsNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, int particleIndex, 
													btFluidParticles& particles, btFluidSphSolverPBF::PbfParticles& pbfData)
{
	int i = particleIndex;
	
	for(int j = 0; j < pbfData.m_neighborTable[i].numNeighbors(); j++) 
	{
		int n = pbfData.m_neighborTable[i].getNeighborIndex(j);
		
		//During the first iteration, the neighbor table contains distance using particles.m_pos instead of m_predictedPosition
		//btScalar distance = pbfData.m_neighborTable[i].getDistance(j);
		
		btVector3 difference = (pbfData.m_predictedPosition[i] - pbfData.m_predictedPosition[n]) * FG.m_simulationScale;		
		btScalar distanceSquared = difference.length2();
		btScalar distance = btSqrt(distanceSquared);
		distance = (distance < SIMD_EPSILON) ? SIMD_EPSILON : distance;
		pbfData.m_neighborTable[i].updateDistance(j, distance);
			
		if(FG.m_sphRadiusSquared > distanceSquared)
		{
			btScalar squaredCloseness = FG.m_sphRadiusSquared - distanceSquared;
			
			btScalar poly6KernelPartialResult = squaredCloseness * squaredCloseness * squaredCloseness;
			pbfData.m_density[i] += poly6KernelPartialResult;
			pbfData.m_density[n] += poly6KernelPartialResult;
		}
	}
}
void computeSphSumsInCellSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL,
										int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverPBF::PbfParticles& pbfData)
{
	BT_PROFILE("computeSphSumsInCellSymmetric()");
	
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computeSphSumsNeighborTableSymmetric(FG, FL, particleIndex, particles, pbfData);
	}
}

void calculateScalingFactorDenomNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, int particleIndex, 
														btFluidParticles& particles, btFluidSphSolverPBF::PbfParticles& pbfData)
{
	int i = particleIndex;
	
	for(int j = 0; j < pbfData.m_neighborTable[i].numNeighbors(); j++) 
	{
		int n = pbfData.m_neighborTable[i].getNeighborIndex(j);
		btScalar distance = pbfData.m_neighborTable[i].getDistance(j);
		if(distance >= FG.m_sphSmoothRadius) continue;
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		
		btScalar spikyKernelGradientScalar = closeness * closeness;
		btVector3 simScaleNormal = (pbfData.m_predictedPosition[i] - pbfData.m_predictedPosition[n]) * (FG.m_simulationScale / distance);
		btVector3 spikyKernelGradient_in = simScaleNormal * spikyKernelGradientScalar;
		
		btScalar spikyKernGradLenSq_in = spikyKernelGradient_in.dot(spikyKernelGradient_in);
		//btScalar spikyKernGradLenSq_ni = (-spikyKernelGradient_in).dot(-spikyKernelGradient_in);	//spikyKernGradLenSq_in == spikyKernGradLenSq_ni
		
		pbfData.m_scalingFactorDenominator[i] += spikyKernGradLenSq_in;
		pbfData.m_scalingFactorDenominator[n] += spikyKernGradLenSq_in;
	}
}
void calculateScalingFactorDenomInCellSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL,
										int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverPBF::PbfParticles& pbfData)
{
	BT_PROFILE("calculateScalingFactorDenomInCellSymmetric()");
	
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		calculateScalingFactorDenomNeighborTableSymmetric(FG, FL, particleIndex, particles, pbfData);
	}
}

void calculateDeltaPositionNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, int particleIndex, 
														btFluidParticles& particles, btFluidSphSolverPBF::PbfParticles& pbfData)
{
	int i = particleIndex;
	
	for(int j = 0; j < pbfData.m_neighborTable[i].numNeighbors(); j++) 
	{
		int n = pbfData.m_neighborTable[i].getNeighborIndex(j);
		btScalar distance = pbfData.m_neighborTable[i].getDistance(j);
		if(distance >= FG.m_sphSmoothRadius) continue;
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		
		btScalar spikyKernelGradientScalar = closeness * closeness;
		btVector3 simScaleNormal = (pbfData.m_predictedPosition[i] - pbfData.m_predictedPosition[n]) * (FG.m_simulationScale / distance);
		btVector3 spikyKernelGradient_in = simScaleNormal * spikyKernelGradientScalar;
		
		btScalar scalingFactorSum = (pbfData.m_scalingFactor[i] + pbfData.m_scalingFactor[n]);
		
		const bool ARTIFICIAL_PRESSURE = true;
		if(ARTIFICIAL_PRESSURE)
		{
			btScalar squaredCloseness = FG.m_sphRadiusSquared - distance * distance;
			btScalar poly6KernResult = FG.m_poly6KernCoeff * squaredCloseness * squaredCloseness * squaredCloseness;
		
			const btScalar k(0.001 * 1e-6);
			btScalar deltaQLength = (btScalar(0.0) * FG.m_sphSmoothRadius);		//Range: [0.0, 0.3] * FG.m_sphSmoothRadius
			
			btScalar squaredClosenessQ = FG.m_sphRadiusSquared - deltaQLength * deltaQLength;
			btScalar poly6KernResultQ = FG.m_poly6KernCoeff * squaredClosenessQ * squaredClosenessQ * squaredClosenessQ;
		
			btScalar s_corr_fraction = (poly6KernResult / poly6KernResultQ);
			btScalar s_corr = k * (s_corr_fraction * s_corr_fraction) * (s_corr_fraction * s_corr_fraction);
			
			scalingFactorSum += s_corr;
		}
		
		btVector3 deltaPosition_in = spikyKernelGradient_in * scalingFactorSum;
		
		pbfData.m_deltaPosition[i] += deltaPosition_in;
		pbfData.m_deltaPosition[n] += -deltaPosition_in;
	}
}
void calculateDeltaPositionInCellSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL,
										int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverPBF::PbfParticles& pbfData)
{
	BT_PROFILE("calculateDeltaPositionInCellSymmetric()");
	
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		calculateDeltaPositionNeighborTableSymmetric(FG, FL, particleIndex, particles, pbfData);
	}
}

void calculateXsphViscosityNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, int particleIndex, 
														btFluidParticles& particles, btFluidSphSolverPBF::PbfParticles& pbfData)
{
	int i = particleIndex;
	
	for(int j = 0; j < pbfData.m_neighborTable[i].numNeighbors(); j++) 
	{
		int n = pbfData.m_neighborTable[i].getNeighborIndex(j);
		btScalar distance = pbfData.m_neighborTable[i].getDistance(j);
		if(distance >= FG.m_sphSmoothRadius) continue;
		
		btScalar squaredCloseness = FG.m_sphRadiusSquared - distance * distance;
		btScalar poly6KernPartialResult = squaredCloseness * squaredCloseness * squaredCloseness;
		
		//This is the opposite of IISPH relativeVelocity_in, which is defined as velocity[i] - velocity[n]
		btVector3 relativeVelocity_in = pbfData.m_nextVelocity[n] - pbfData.m_nextVelocity[i];	
		
		btVector3 xsphViscosity_in = relativeVelocity_in * poly6KernPartialResult;
		
		pbfData.m_xsphViscosity[i] += xsphViscosity_in;
		pbfData.m_xsphViscosity[n] += -xsphViscosity_in;
	}
}
void calculateXsphViscosityInCellSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL,
										int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverPBF::PbfParticles& pbfData)
{
	BT_PROFILE("calculateXsphViscosityNeighborTableSymmetric()");
	
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		calculateXsphViscosityNeighborTableSymmetric(FG, FL, particleIndex, particles, pbfData);
	}
}

void btFluidSphSolverPBF::updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids)
{
	BT_PROFILE("btFluidSphSolverPBF::updateGridAndCalculateSphForces()");
	
	//SPH data is discarded/recalculated every frame, so only 1
	//set of arrays are needed if there is no fluid-fluid interaction.
	if( m_pbfData.size() != 1 ) m_pbfData.resize(1);
	
	for(int fluidIndex = 0; fluidIndex < numFluids; ++fluidIndex)
	{
		btFluidSph* fluid = fluids[fluidIndex];
		int numParticles = fluid->numParticles();
		if(!numParticles) continue;
		
		const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
		const btFluidSortingGrid& grid = fluid->getGrid();
		btFluidParticles& particles = fluid->internalGetParticles();
		btFluidSphSolverPBF::PbfParticles& pbfData = m_pbfData[0];
		
		if( numParticles > pbfData.size() ) pbfData.resize(numParticles);
		
		fluid->insertParticlesIntoGrid();
		
		//Apply forces and predict position
		{
			for(int n = 0; n < numParticles; ++n) pbfData.m_predictedVelocity[n] = particles.m_vel[n] + FL.m_gravity * FG.m_timeStep;
			for(int n = 0; n < numParticles; ++n) 
				pbfData.m_predictedPosition[n] = particles.m_pos[n] + pbfData.m_predictedVelocity[n] * (FG.m_timeStep / FG.m_simulationScale);
		}
		
		//Find neighbors / build neighbor tables
		{
			BT_PROFILE("Find neighbors");
		
			for(int n = 0; n < numParticles; ++n) pbfData.m_neighborTable[n].clear();
			
			for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
			{
				const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
				
				for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
					findNeighborsInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, pbfData);
			}
		}
		
		//Solve positions
		const int MAX_ITERATIONS = 5;
		for(int iteration = 0; iteration < MAX_ITERATIONS; ++iteration)
		{
			//Calculate constraint error C and scaling factor s
			{
				//Calculate density and update distances
				{
					const btScalar poly6ZeroDistance = FG.m_sphRadiusSquared * FG.m_sphRadiusSquared * FG.m_sphRadiusSquared;
					const btScalar initialSphSum = poly6ZeroDistance * FL.m_initialSum;
					for(int n = 0; n < numParticles; ++n) pbfData.m_density[n] = initialSphSum;
					
					for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
					{
						const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
						
						for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
							computeSphSumsInCellSymmetric(FG, FL, multithreadingGroup[cell], grid, particles, pbfData);
					}
			
					for(int n = 0; n < numParticles; ++n) pbfData.m_density[n] *= FL.m_sphParticleMass * FG.m_poly6KernCoeff;
				}
				
				//Calculate scaling factor denominator (denominator of equation 9)
				{
					for(int n = 0; n < numParticles; ++n) pbfData.m_scalingFactorDenominator[n] = btScalar(0.0);
					
					for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
					{
						const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
						
						for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
							calculateScalingFactorDenomInCellSymmetric(FG, FL, multithreadingGroup[cell], grid, particles, pbfData);
					}
					
					const btScalar scalingFactorDenomConstants = (FG.m_spikyKernGradCoeff * FG.m_spikyKernGradCoeff) / (FL.m_restDensity * FL.m_restDensity);
					for(int n = 0; n < numParticles; ++n) pbfData.m_scalingFactorDenominator[n] *= scalingFactorDenomConstants;
				}
				
				//Calculate scaling factor s
				const btScalar EPSILON(1e-5);
				for(int n = 0; n < numParticles; ++n)
				{
					btScalar C = (pbfData.m_density[n] / FL.m_restDensity) - btScalar(1.0);
				
					pbfData.m_scalingFactor[n] = C / (pbfData.m_scalingFactorDenominator[n] + EPSILON);
				}
			}
			
			
			//Calculate delta position
			//(in the paper, collision detection and response is also done here)
			{
				for(int n = 0; n < numParticles; ++n) pbfData.m_deltaPosition[n].setValue(0,0,0);
				
				for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
				{
					const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
					
					for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
						calculateDeltaPositionInCellSymmetric(FG, FL, multithreadingGroup[cell], grid, particles, pbfData);
				}
					
				const btScalar deltaPositionConstants = ( btScalar(-1.0) / FL.m_restDensity ) * FG.m_spikyKernGradCoeff;
				for(int n = 0; n < numParticles; ++n) pbfData.m_deltaPosition[n] *= deltaPositionConstants;
			}
			
			//Update position
			for(int n = 0; n < numParticles; ++n) 
				pbfData.m_predictedPosition[n] += pbfData.m_deltaPosition[n] / (FG.m_simulationScale);// * FG.m_simulationScale);
			
			//Project position into an AABB boundary
			if(FL.m_enableAabbBoundary)
			{
				for(int n = 0; n < numParticles; ++n)
				{
					const btVector3& min = FL.m_aabbBoundaryMin;
					const btVector3& max = FL.m_aabbBoundaryMax;
					
					btVector3 clampedPosition( 	btMax( min.x(), btMin(pbfData.m_predictedPosition[n].x(), max.x()) ),
												btMax( min.y(), btMin(pbfData.m_predictedPosition[n].y(), max.y()) ),
												btMax( min.z(), btMin(pbfData.m_predictedPosition[n].z(), max.z()) ) );
												
					pbfData.m_predictedPosition[n] = clampedPosition;
				}
			}
		}
		
		//Update velocity
		for(int n = 0; n < numParticles; ++n) 
			pbfData.m_nextVelocity[n] = (pbfData.m_predictedPosition[n] - particles.m_pos[n]) * (FG.m_simulationScale / FG.m_timeStep);
		
		//Apply vorticity confinement(not yet implemented) and XSPH viscosity
		{
			for(int n = 0; n < numParticles; ++n) pbfData.m_xsphViscosity[n].setValue(0,0,0);
			
			for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
			{
				const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
				
				for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
					calculateXsphViscosityInCellSymmetric(FG, FL, multithreadingGroup[cell], grid, particles, pbfData);
			}
				
			const btScalar c(1e-8);
			const btScalar xsphViscosityConstants = c * FG.m_poly6KernCoeff;
			for(int n = 0; n < numParticles; ++n) pbfData.m_xsphViscosity[n] *= xsphViscosityConstants;
			
			for(int n = 0; n < numParticles; ++n) pbfData.m_nextVelocity[n] += pbfData.m_xsphViscosity[n];
		}
		
		//Write new velocity and position
		for(int n = 0; n < numParticles; ++n) particles.m_vel[n] = pbfData.m_nextVelocity[n];
		for(int n = 0; n < numParticles; ++n) particles.m_vel_eval[n] = pbfData.m_nextVelocity[n];
		for(int n = 0; n < numParticles; ++n) particles.m_pos[n] = pbfData.m_predictedPosition[n];
	}
}

void btFluidSphSolverPBF::findNeighborsInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
														const btFluidSortingGrid& grid, btFluidParticles& particles,
														btFluidSphSolverPBF::PbfParticles& pbfData)
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
						const btScalar UNUSED_DISTANCE(BT_LARGE_FLOAT);		//Distance is recalculated during each solver iteration
						
						if( !pbfData.m_neighborTable[i].isFilled() ) pbfData.m_neighborTable[i].addNeighbor(n, UNUSED_DISTANCE);
						else if( !pbfData.m_neighborTable[n].isFilled() ) pbfData.m_neighborTable[n].addNeighbor(i, UNUSED_DISTANCE);
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
