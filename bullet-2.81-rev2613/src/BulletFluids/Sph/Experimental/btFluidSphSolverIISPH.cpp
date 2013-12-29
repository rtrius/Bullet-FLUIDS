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
#include "btFluidSphSolverIISPH.h"

#include "BulletFluids/Sph/btFluidSortingGrid.h"

#include "LinearMath/btQuickprof.h"		//BT_PROFILE(name) macro

/*
inline btScalar computePressure(btScalar restDensity, btScalar stiffness, btScalar particleDensity)
{
	const btScalar KAPPA(stiffness);
	const btScalar GAMMA(1.0);	//For WCSPH gamma == 7
	
	return ( (KAPPA * restDensity) / GAMMA ) * ( btPow(particleDensity / restDensity, GAMMA) - btScalar(1.0) );
}
*/

void computeViscosityForceAndDiiNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, int particleIndex, 
													btFluidParticles& particles, btFluidSphSolverIISPH::IiSphParticles& iiSphData)
{
	int i = particleIndex;
	
	for(int j = 0; j < iiSphData.m_neighborTable[i].numNeighbors(); j++) 
	{
		int n = iiSphData.m_neighborTable[i].getNeighborIndex(j);
		btScalar distance = iiSphData.m_neighborTable[i].getDistance(j);
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		
		//Compute viscosity force
		{
			btScalar viscosityScalar = closeness / (iiSphData.m_density[i] * iiSphData.m_density[n]);
			btVector3 viscosityForce = (particles.m_vel[n] - particles.m_vel[i]) * viscosityScalar;
			
			iiSphData.m_viscosityAcceleration[i] += viscosityForce;
			iiSphData.m_viscosityAcceleration[n] += -viscosityForce;
		}
		
		//Compute d_ii
		{
			btScalar spikyKernelGradientScalar = closeness * closeness;
			btVector3 simScaleNormal = (particles.m_pos[i] - particles.m_pos[n]) * (FG.m_simulationScale / distance);
			btVector3 spikyKernelGradient_in = simScaleNormal * spikyKernelGradientScalar;
			
			iiSphData.m_d_ii[i] += spikyKernelGradient_in;
			iiSphData.m_d_ii[n] += -spikyKernelGradient_in;
		}
	}
}
void computeViscosityForceAndDiiInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
										const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverIISPH::IiSphParticles& iiSphData)
{
	BT_PROFILE("computeViscosityForceAndDiiInCellSymmetric()");
	
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computeViscosityForceAndDiiNeighborTableSymmetric(FG, particleIndex, particles, iiSphData);
	}
}

void computeDensityAdvAndAiiNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, int particleIndex, 
													btFluidParticles& particles, btFluidSphSolverIISPH::IiSphParticles& iiSphData)
{
	int i = particleIndex;
	
	for(int j = 0; j < iiSphData.m_neighborTable[i].numNeighbors(); j++) 
	{
		int n = iiSphData.m_neighborTable[i].getNeighborIndex(j);
		btScalar distance = iiSphData.m_neighborTable[i].getDistance(j);
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		
		btScalar spikyKernelGradientScalar = closeness * closeness;
		btVector3 simScaleNormal = (particles.m_pos[i] - particles.m_pos[n]) * (FG.m_simulationScale / distance);
		btVector3 spikyKernelGradient_in = simScaleNormal * spikyKernelGradientScalar;
		
		//Compute density_adv
		{
			btVector3 relativePredictedVelocity_in = (iiSphData.m_predictedVelocity[i] - iiSphData.m_predictedVelocity[n]);
			
			btScalar densityAdv_in = relativePredictedVelocity_in.dot(spikyKernelGradient_in);
			//btScalar densityAdv_ni = (-relativePredictedVelocity_in).dot(-spikyKernelGradient_in);	//densityAdv_in == densityAdv_ni
			
			iiSphData.m_density_adv[i] += densityAdv_in;
			iiSphData.m_density_adv[n] += densityAdv_in;
		}
		
		//Compute a_ii
		{
			btScalar d_in_ni_numerator = -FG.m_timeStep * FG.m_timeStep * FG.m_spikyKernGradCoeff * FL.m_sphParticleMass;
		
			btScalar d_ni_scalar = d_in_ni_numerator / (iiSphData.m_density[i] * iiSphData.m_density[i]);
			btScalar d_in_scalar = d_in_ni_numerator / (iiSphData.m_density[n] * iiSphData.m_density[n]);
			
			btVector3 d_ni = -spikyKernelGradient_in * d_ni_scalar;
			btVector3 d_in = spikyKernelGradient_in * d_in_scalar;
			
			iiSphData.m_a_ii[i] += (iiSphData.m_d_ii[i] - d_ni).dot(spikyKernelGradient_in);
			iiSphData.m_a_ii[n] += (iiSphData.m_d_ii[n] - d_in).dot(-spikyKernelGradient_in);
		}
	}
}
void computeDensityAdvAndAiiInCellSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL,
										int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverIISPH::IiSphParticles& iiSphData)
{
	BT_PROFILE("computeDensityAdvAndAiiInCellSymmetric()");
	
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computeDensityAdvAndAiiNeighborTableSymmetric(FG, FL, particleIndex, particles, iiSphData);
	}
}

void computeDijPjSumNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, int particleIndex, 
													btFluidParticles& particles, btFluidSphSolverIISPH::IiSphParticles& iiSphData)
{
	int i = particleIndex;
	
	for(int j = 0; j < iiSphData.m_neighborTable[i].numNeighbors(); j++) 
	{
		int n = iiSphData.m_neighborTable[i].getNeighborIndex(j);
		btScalar distance = iiSphData.m_neighborTable[i].getDistance(j);
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		
		btScalar spikyKernelGradientScalar = closeness * closeness;
		btVector3 simScaleNormal = (particles.m_pos[i] - particles.m_pos[n]) * (FG.m_simulationScale / distance);
		
		btVector3 spikyKernelGradient = simScaleNormal * spikyKernelGradientScalar;
		
		btScalar scalar_i = iiSphData.m_pressure[n] / (iiSphData.m_density[n] * iiSphData.m_density[n]);
		btScalar scalar_n = iiSphData.m_pressure[i] / (iiSphData.m_density[i] * iiSphData.m_density[i]);
		
		iiSphData.m_d_ij_pj_sum[i] += spikyKernelGradient * scalar_i;
		iiSphData.m_d_ij_pj_sum[n] += -spikyKernelGradient * scalar_n;
	}
}
void computeDijPjSumInCellSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL,
										int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverIISPH::IiSphParticles& iiSphData)
{
	BT_PROFILE("computeDijPjSumInCellSymmetric()");
	
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computeDijPjSumNeighborTableSymmetric(FG, FL, particleIndex, particles, iiSphData);
	}
}


void computeEquation13SumNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, int particleIndex, 
													btFluidParticles& particles, btFluidSphSolverIISPH::IiSphParticles& iiSphData)
{
	int i = particleIndex;
	
	for(int j = 0; j < iiSphData.m_neighborTable[i].numNeighbors(); j++) 
	{
		int n = iiSphData.m_neighborTable[i].getNeighborIndex(j);
		btScalar distance = iiSphData.m_neighborTable[i].getDistance(j);
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		
		btScalar spikyKernelGradientScalar = closeness * closeness;
		btVector3 simScaleNormal = (particles.m_pos[i] - particles.m_pos[n]) * (FG.m_simulationScale / distance);
		btVector3 spikyKernelGradient_in = simScaleNormal * spikyKernelGradientScalar;
		
		btScalar d_in_ni_numerator = -FG.m_timeStep * FG.m_timeStep * FG.m_spikyKernGradCoeff * FL.m_sphParticleMass;
		btScalar d_ni_scalar = d_in_ni_numerator / (iiSphData.m_density[i] * iiSphData.m_density[i]);
		btScalar d_in_scalar = d_in_ni_numerator / (iiSphData.m_density[n] * iiSphData.m_density[n]);
			
		btVector3 d_ni = -spikyKernelGradient_in * d_ni_scalar;
		btVector3 d_in = spikyKernelGradient_in * d_in_scalar;
		
		btVector3 d_ni_pi = d_ni * iiSphData.m_pressure[i];
		btVector3 d_in_pn = d_in * iiSphData.m_pressure[n];
		
		btVector3 i_term = iiSphData.m_d_ij_pj_sum[i] - iiSphData.m_d_ii[n]*iiSphData.m_pressure[n] - (iiSphData.m_d_ij_pj_sum[n] - d_ni_pi);
		btVector3 n_term = iiSphData.m_d_ij_pj_sum[n] - iiSphData.m_d_ii[i]*iiSphData.m_pressure[i] - (iiSphData.m_d_ij_pj_sum[i] - d_in_pn);
		
		iiSphData.m_equation13_sum[i] += i_term.dot(spikyKernelGradient_in);
		iiSphData.m_equation13_sum[n] += n_term.dot(-spikyKernelGradient_in);
	}
}
void computeEquation13SumInCellSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL,
										int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverIISPH::IiSphParticles& iiSphData)
{
	BT_PROFILE("computeEquation13SumInCellSymmetric()");
	
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computeEquation13SumNeighborTableSymmetric(FG, FL, particleIndex, particles, iiSphData);
	}
}

void computePressureForceNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, int particleIndex, 
													btFluidParticles& particles, btFluidSphSolverIISPH::IiSphParticles& iiSphData)
{
	int i = particleIndex;
	
	for(int j = 0; j < iiSphData.m_neighborTable[i].numNeighbors(); j++) 
	{
		int n = iiSphData.m_neighborTable[i].getNeighborIndex(j);
		btScalar distance = iiSphData.m_neighborTable[i].getDistance(j);
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		
		btScalar spikyKernelGradientScalar = closeness * closeness;
		btVector3 simScaleNormal = (particles.m_pos[i] - particles.m_pos[n]) * (FG.m_simulationScale / distance);
		btVector3 spikyKernelGradient = simScaleNormal * spikyKernelGradientScalar;
		
		btScalar pterm_i = iiSphData.m_pressure[i] / (iiSphData.m_density[i] * iiSphData.m_density[i]);
		btScalar pterm_n = iiSphData.m_pressure[n] / (iiSphData.m_density[n] * iiSphData.m_density[n]);
		
		btVector3 pressureAcceleration = spikyKernelGradient * (pterm_i + pterm_n);
		
		//btScalar alternatePressureScalar = btScalar(0.5) * (iiSphData.m_pressure[i] + iiSphData.m_pressure[n]) / (iiSphData.m_density[i] * iiSphData.m_density[n]);
		//btVector3 pressureAcceleration = spikyKernelGradient * alternatePressureScalar;
		
		iiSphData.m_pressureAcceleration[i] += pressureAcceleration;
		iiSphData.m_pressureAcceleration[n] += -pressureAcceleration;
	}
}
void computePressureForceInCellSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL,
										int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
										btFluidSphSolverIISPH::IiSphParticles& iiSphData)
{
	BT_PROFILE("computePressureForceInCellSymmetric()");
	
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computePressureForceNeighborTableSymmetric(FG, FL, particleIndex, particles, iiSphData);
	}
}

void btFluidSphSolverIISPH::updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids)
{
	BT_PROFILE("btFluidSphSolverIISPH::updateGridAndCalculateSphForces()");
	
	//SPH data is discarded/recalculated every frame, so only 1
	//set of arrays are needed if there is no fluid-fluid interaction.
	if( m_iiSphdata.size() != 1 ) m_iiSphdata.resize(1);
	
	for(int fluidIndex = 0; fluidIndex < numFluids; ++fluidIndex)
	{
		btFluidSph* fluid = fluids[fluidIndex];
		int numParticles = fluid->numParticles();
		if(!numParticles) continue;
		
		const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
		const btFluidSortingGrid& grid = fluid->getGrid();
		btFluidParticles& particles = fluid->internalGetParticles();
		btFluidSphSolverIISPH::IiSphParticles& iiSphData = m_iiSphdata[0];
		
		if( numParticles > iiSphData.size() ) iiSphData.resize(numParticles);
		
		fluid->insertParticlesIntoGrid();
		
		//Predict Advection
		{
			BT_PROFILE("Predict Advection");
		
			//Compute current density, and build neighbor tables
			{
				BT_PROFILE("Compute density and get neighbors");
			
				for(int n = 0; n < numParticles; ++n) iiSphData.m_density[n] = FG.m_initialSum;
				for(int n = 0; n < numParticles; ++n) iiSphData.m_neighborTable[n].clear();
				
				{
					BT_PROFILE("compute sums");
					
					for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
					{
						const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
						if( !multithreadingGroup.size() ) continue;
						
						for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
							calculateSumsInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, iiSphData);
					}
				}
				
				for(int n = 0; n < numParticles; ++n) iiSphData.m_density[n] *= FL.m_sphParticleMass * FG.m_poly6KernCoeff;
			}
			
			//Predict next velocity (from viscosity, gravity, surface tension, collision forces)
			{
				BT_PROFILE("Predict next velocity");
			
				for(int n = 0; n < fluid->numParticles(); ++n) iiSphData.m_viscosityAcceleration[n].setValue(0,0,0);
				for(int n = 0; n < fluid->numParticles(); ++n) iiSphData.m_d_ii[n].setValue(0,0,0);
					
				//Compute viscosity force and d_ii
 				{
					for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
					{
						const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
					
						for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
							computeViscosityForceAndDiiInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, iiSphData);
					}
					
					const btScalar viscosityForceConstants = FG.m_viscosityKernLapCoeff * FL.m_viscosity * FL.m_sphParticleMass;
					for(int n = 0; n < fluid->numParticles(); ++n) iiSphData.m_viscosityAcceleration[n] *= viscosityForceConstants;
					
					const btScalar d_ii_constants = -FG.m_timeStep * FG.m_timeStep * FL.m_sphParticleMass * FG.m_spikyKernGradCoeff;
					for(int n = 0; n < fluid->numParticles(); ++n) 
						iiSphData.m_d_ii[n] *= d_ii_constants / (iiSphData.m_density[n] * iiSphData.m_density[n]);
				}
				
				for(int n = 0; n < numParticles; ++n)
				{
					btVector3 predictedAcceleration = iiSphData.m_viscosityAcceleration[n] + FL.m_gravity;
				
					iiSphData.m_predictedVelocity[n] = particles.m_vel[n] + predictedAcceleration * FG.m_timeStep;
				}
			}
			
			//Compute density_adv and a_ii
			{
				BT_PROFILE("Compute density_adv, a_ii");
			
				for(int n = 0; n < numParticles; ++n) iiSphData.m_density_adv[n] = btScalar(0.0);
				for(int n = 0; n < numParticles; ++n) iiSphData.m_a_ii[n] = btScalar(0.0);
				
				for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
				{
					const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
					
					for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
						computeDensityAdvAndAiiInCellSymmetric(FG, FL, multithreadingGroup[cell], grid, particles, iiSphData);
				}
					
				const btScalar density_adv_constants = FL.m_sphParticleMass * FG.m_timeStep * FG.m_spikyKernGradCoeff;
				for(int n = 0; n < numParticles; ++n) 
					iiSphData.m_density_adv[n] = iiSphData.m_density[n] + iiSphData.m_density_adv[n] * density_adv_constants;
					
					
				const btScalar a_ii_constants = FL.m_sphParticleMass * FG.m_spikyKernGradCoeff;
				for(int n = 0; n < numParticles; ++n) iiSphData.m_a_ii[n] *= a_ii_constants;
			}
			
			
			//Compute initial pressure
			//for(int n = 0; n < numParticles; ++n) iiSphData.m_pressure[n] = FL.m_stiffness * (iiSphData.m_density_adv[n] - FL.m_restDensity);
			for(int n = 0; n < numParticles; ++n) iiSphData.m_pressure[n] = btScalar(0.0);
			//for(int n = 0; n < numParticles; ++n) iiSphData.m_pressure[n] = btScalar(0.5) * m_prevPressure[n];
		}
		
		//Pressure Solve
		if(1)
		{
			BT_PROFILE("Solve for pressure");
		
			const int MAX_ITERATIONS = 30;
			for(int iteration = 0; iteration < MAX_ITERATIONS; ++iteration)
			{
				//Loop 1 - compute sum{d_ij * p_j}
				{
					BT_PROFILE("Compute d_ij_pj sum");
				
					for(int n = 0; n < fluid->numParticles(); ++n) iiSphData.m_d_ij_pj_sum[n].setValue(0,0,0);
					
					for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
					{
						const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
						
						for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
							computeDijPjSumInCellSymmetric(FG, FL, multithreadingGroup[cell], grid, particles, iiSphData);
					}
					
					const btScalar d_ij_pj_sum_scalar = -FG.m_timeStep * FG.m_timeStep * FL.m_sphParticleMass * FG.m_spikyKernGradCoeff;
					for(int n = 0; n < fluid->numParticles(); ++n) iiSphData.m_d_ij_pj_sum[n] *= d_ij_pj_sum_scalar;
				}
				
				//Loop 2 - update pressure (equation 13)
				{
					BT_PROFILE("Update pressure");
				
					for(int n = 0; n < fluid->numParticles(); ++n) iiSphData.m_equation13_sum[n] = btScalar(0.0);
					
					for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
					{
						const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
						
						for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
							computeEquation13SumInCellSymmetric(FG, FL, multithreadingGroup[cell], grid, particles, iiSphData);
					}
					
					const btScalar equation13_sum_scalar = FL.m_sphParticleMass * FG.m_spikyKernGradCoeff;
					for(int n = 0; n < fluid->numParticles(); ++n) iiSphData.m_equation13_sum[n] *= equation13_sum_scalar;
					
					//Jacobi Iteration on Pressure
					const bool ALLOW_PRESSURE_FORCE_ATTRACTION = false;
					const btScalar OMEGA(0.5);
					for(int n = 0; n < fluid->numParticles(); ++n) 
					{
						btScalar b = FL.m_restDensity - iiSphData.m_density_adv[n];
						if(!ALLOW_PRESSURE_FORCE_ATTRACTION) b = btMin( btScalar(0.0), b );
						
						iiSphData.m_pressure[n] = (btScalar(1.0) - OMEGA) * iiSphData.m_pressure[n] 
							+ (OMEGA / iiSphData.m_a_ii[n]) * (b - iiSphData.m_equation13_sum[n]);
					}
				}
			}
		}
		
		//Compute Pressure Force
		{
			BT_PROFILE("Compute pressure force");
		
			for(int n = 0; n < fluid->numParticles(); ++n) iiSphData.m_pressureAcceleration[n].setValue(0,0,0);
			
			for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
			{
				const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
					
				for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
					computePressureForceInCellSymmetric(FG, FL, multithreadingGroup[cell], grid, particles, iiSphData);
			}
				
			const btScalar pressureAccelerationScalar = -FL.m_sphParticleMass * FG.m_spikyKernGradCoeff;
			for(int n = 0; n < fluid->numParticles(); ++n) iiSphData.m_pressureAcceleration[n] *= pressureAccelerationScalar;
		}

		//Integrate
		{
			//Apply SPH force to particles
			//Gravity is applied during velocity integration(after btFluidSphSolverIISPH::updateGridAndCalculateSphForces() is called)
			for(int n = 0; n < numParticles; ++n) 
			{
				btVector3 sphAcceleration = iiSphData.m_viscosityAcceleration[n] + iiSphData.m_pressureAcceleration[n];
			
				fluid->applyForce(n, sphAcceleration * FL.m_particleMass);
			}
		}
	}
}

void btFluidSphSolverIISPH::calculateSumsInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
														const btFluidSortingGrid& grid, btFluidParticles& particles,
														btFluidSphSolverIISPH::IiSphParticles& sphData)
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
						sphData.m_density[i] += poly6KernPartialResult;
						sphData.m_density[n] += poly6KernPartialResult;
						
						btScalar distance = btSqrt(distanceSquared);
						distance = (distance < SIMD_EPSILON) ? SIMD_EPSILON : distance;
						
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
