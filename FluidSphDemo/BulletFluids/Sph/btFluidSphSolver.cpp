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

//Link to e.g. 'BulletMultiThreaded.lib' if enabling this
//(must build Bullet with CMake to get BulletMultiThreaded library)
//#define FLUIDS_MULTITHREADED_ENABLED
#ifdef FLUIDS_MULTITHREADED_ENABLED
#include "btParallelFor.h"

const unsigned int NUM_THREADS = 4;
btParallelFor parallelFor("parallelForThreads", NUM_THREADS);

void calculateSumsInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, const btFluidSortingGrid& grid, 
								btFluidParticles& particles, btFluidSphSolverDefault::SphParticles& sphData);	
void calculateForcesInCellSymmetric(const btFluidSphParametersGlobal& FG, const btScalar vterm,
									int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
									btFluidSphSolverDefault::SphParticles& sphData);		
					   
struct PF_ComputePressureData
{
	const btFluidSphParametersGlobal& m_globalParameters; 
	const btAlignedObjectArray<int>& m_gridCellGroup;
	const btFluidSortingGrid& m_grid;
	btFluidParticles& m_particles; 
	btFluidSphSolverDefault::SphParticles& m_sphData;
	
	PF_ComputePressureData(const btFluidSphParametersGlobal& FG, const btAlignedObjectArray<int>& gridCellGroup,
							const btFluidSortingGrid& grid, btFluidParticles& particles, btFluidSphSolverDefault::SphParticles& sphData) 
	: m_globalParameters(FG), m_gridCellGroup(gridCellGroup), m_grid(grid), m_particles(particles), m_sphData(sphData) {}
};
void PF_ComputePressureFunction(void* parameters, int index)
{
	PF_ComputePressureData* data = static_cast<PF_ComputePressureData*>(parameters);
	
	calculateSumsInCellSymmetric(data->m_globalParameters, data->m_gridCellGroup[index], data->m_grid, data->m_particles, data->m_sphData);
}

struct PF_ComputeForceData
{
	const btFluidSphParametersGlobal& m_globalParameters; 
	const btScalar m_vterm;
	const btAlignedObjectArray<int>& m_gridCellGroup;
	const btFluidSortingGrid& m_grid;
	btFluidParticles& m_particles;
	btFluidSphSolverDefault::SphParticles& m_sphData;
	
	PF_ComputeForceData(const btFluidSphParametersGlobal& FG, const btScalar vterm, const btAlignedObjectArray<int>& gridCellGroup, 
						const btFluidSortingGrid& grid, btFluidParticles& particles, btFluidSphSolverDefault::SphParticles& sphData) 
	: m_globalParameters(FG), m_vterm(vterm),  m_gridCellGroup(gridCellGroup), 
	m_grid(grid), m_particles(particles), m_sphData(sphData) {}
};
void PF_ComputeForceFunction(void* parameters, int index)
{
	PF_ComputeForceData* data = static_cast<PF_ComputeForceData*>(parameters);
	
	calculateForcesInCellSymmetric(data->m_globalParameters, data->m_vterm, data->m_gridCellGroup[index], 
									data->m_grid, data->m_particles, data->m_sphData);
}
#endif //FLUIDS_MULTITHREADED_ENABLED


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
void applyBoundaryForceToParticle(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, btScalar simScaleParticleRadius,
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
	
	particles.m_accumulatedForce[i] += acceleration * FL.m_particleMass;
}
void btFluidSphSolver::applyBoundaryForcesSingleFluid(const btFluidSphParametersGlobal& FG, btFluidSph* fluid)
{
	BT_PROFILE("btFluidSphSolver::applyBoundaryForcesSingleFluid()");
	
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	btFluidParticles& particles = fluid->internalGetParticles();
	
	const btScalar simScaleParticleRadius = FL.m_particleRadius * FG.m_simulationScale;

	for(int i = 0; i < particles.size(); ++i) applyBoundaryForceToParticle(FG, FL, simScaleParticleRadius, particles, i);
}


inline void resolveAabbCollision_impulse(const btFluidSphParametersLocal& FL, const btVector3& velocity, 
										const btVector3& normal, btScalar distance, btVector3* impulse)
{
	if( distance < btScalar(0.0) )	//Negative distance indicates penetration
	{
		btScalar penetratingMagnitude = velocity.dot(-normal);
		if( penetratingMagnitude < btScalar(0.0) ) penetratingMagnitude = btScalar(0.0);
		
		btVector3 penetratingVelocity = -normal * penetratingMagnitude;
		btVector3 tangentialVelocity = velocity - penetratingVelocity;
		
		penetratingVelocity *= btScalar(1.0) + FL.m_boundaryRestitution;
		
		btScalar positionError = (-distance) * FL.m_boundaryErp;
		penetratingVelocity += -normal * positionError;
		
		*impulse -= penetratingVelocity;
		*impulse -= tangentialVelocity * FL.m_boundaryFriction;
	}
}
void applyBoundaryImpulseToParticle(btScalar simulationScale, btScalar simScaleParticleRadius, const btFluidSphParametersLocal& FL, 
									btFluidParticles& particles, int particleIndex)
{
	int i = particleIndex;
	
	const btScalar radius = simScaleParticleRadius;
	const btScalar simScale = simulationScale;
	
	const btVector3& boundaryMin = FL.m_aabbBoundaryMin;
	const btVector3& boundaryMax = FL.m_aabbBoundaryMax;
	
	const btVector3& pos = particles.m_pos[i];
	btVector3& vel = particles.m_vel[i];
	btVector3& vel_eval = particles.m_vel_eval[i];
	
	btVector3 impulse(0.f, 0.f, 0.f);
	resolveAabbCollision_impulse( FL, vel, btVector3( 1.0, 0.0, 0.0), ( pos.x() - boundaryMin.x() )*simScale - radius, &impulse );
	resolveAabbCollision_impulse( FL, vel, btVector3(-1.0, 0.0, 0.0), ( boundaryMax.x() - pos.x() )*simScale - radius, &impulse );
	resolveAabbCollision_impulse( FL, vel, btVector3(0.0,  1.0, 0.0), ( pos.y() - boundaryMin.y() )*simScale - radius, &impulse );
	resolveAabbCollision_impulse( FL, vel, btVector3(0.0, -1.0, 0.0), ( boundaryMax.y() - pos.y() )*simScale - radius, &impulse );
	resolveAabbCollision_impulse( FL, vel, btVector3(0.0, 0.0,  1.0), ( pos.z() - boundaryMin.z() )*simScale - radius, &impulse );
	resolveAabbCollision_impulse( FL, vel, btVector3(0.0, 0.0, -1.0), ( boundaryMax.z() - pos.z() )*simScale - radius, &impulse );
	
	//Leapfrog integration
	btVector3 vnext = vel + impulse;
	vel_eval = (vel + vnext) * btScalar(0.5);
	vel = vnext;
}
void btFluidSphSolver::applyBoundaryImpulsesSingleFluid(const btFluidSphParametersGlobal& FG, btFluidSph* fluid)
{
	BT_PROFILE("btFluidSphSolver::applyBoundaryImpulsesSingleFluid()");
	
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	btFluidParticles& particles = fluid->internalGetParticles();
	
	const btScalar simScaleParticleRadius = FL.m_particleRadius * FG.m_simulationScale;

	for(int i = 0; i < particles.size(); ++i)
	{
		applyBoundaryImpulseToParticle(FG.m_simulationScale, simScaleParticleRadius, FL, particles, i);
	}
}

		
void calculateSumsInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, const btFluidSortingGrid& grid, 
								btFluidParticles& particles, btFluidSphSolverDefault::SphParticles& sphData)
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

void btFluidSphSolverDefault::sphComputePressure(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, btFluidSphSolverDefault::SphParticles& sphData)
{
	BT_PROFILE("btFluidSphSolverDefault::sphComputePressure()");
	
	const int numParticles = fluid->numParticles();
	
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	const btFluidSortingGrid& grid = fluid->getGrid();
	btFluidParticles& particles = fluid->internalGetParticles();
	
	{
		BT_PROFILE("sphComputePressure() - reset sums, clear table");
		
		for(int i = 0; i < numParticles; ++i) sphData.m_invDensity[i] = FG.m_initialSum;
		for(int i = 0; i < numParticles; ++i) sphData.m_neighborTable[i].clear();
	}
	
	{
		BT_PROFILE("sphComputePressure() - compute sums");
		
		for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
		{
			const btAlignedObjectArray<int>& currentGroup = grid.internalGetMultithreadingGroup(group);
			if( !currentGroup.size() ) continue;
		
#ifdef FLUIDS_MULTITHREADED_ENABLED
			PF_ComputePressureData PressureData(FG, currentGroup, grid, particles, sphData);
			parallelFor.execute( PF_ComputePressureFunction, &PressureData, 0, currentGroup.size() - 1, 1 );
#else
			for(int cell = 0; cell < currentGroup.size(); ++cell)
				calculateSumsInCellSymmetric(FG, currentGroup[cell], grid, particles, sphData);
#endif		
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
void calculateForcesInCellSymmetric(const btFluidSphParametersGlobal& FG, const btScalar vterm,
									int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
									btFluidSphSolverDefault::SphParticles& sphData)
{
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
	{
		computeForceNeighborTableSymmetric(FG, vterm, i, particles, sphData);
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
		
#ifdef FLUIDS_MULTITHREADED_ENABLED
		PF_ComputeForceData ForceData(FG, vterm, currentGroup, grid, particles, sphData);
		parallelFor.execute( PF_ComputeForceFunction, &ForceData, 0, currentGroup.size() - 1, 1 );
#else
		for(int cell = 0; cell < currentGroup.size(); ++cell)
			calculateForcesInCellSymmetric(FG, vterm, currentGroup[cell], grid, particles, sphData);
#endif		
	}
	
	for(int i = 0; i < particles.size(); ++i)sphData.m_sphForce[i] *= FL.m_sphParticleMass;
}