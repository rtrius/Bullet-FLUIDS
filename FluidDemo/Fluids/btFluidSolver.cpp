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
#include "btFluidSolver.h"

#include "btFluidSortingGrid.h"

//Link to e.g. 'BulletMultiThreaded.lib' if enabling this
//(must build Bullet with CMake to get BulletMultiThreaded library)
//#define FLUIDS_MULTITHREADED_ENABLED
#ifdef FLUIDS_MULTITHREADED_ENABLED
#include "btParallelFor.h"

const unsigned int NUM_THREADS = 4;
btParallelFor parallelFor("parallelForThreads", NUM_THREADS);

void calculateSumsInCellSymmetric(const btFluidParametersGlobal& FG, int gridCellIndex, const btFluidSortingGrid& grid, 
								btFluidParticles& particles, btFluidSolverSph::SphParticles& sphData);	
void calculateForcesInCellSymmetric(const btFluidParametersGlobal& FG, const btScalar vterm,
									int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
									btFluidSolverSph::SphParticles& sphData);
void integrateParticle(const btFluidParametersGlobal& FG, const btFluidParametersLocal& FL,
					   btScalar speedLimitSquared, btScalar R2, int particleIndex, btFluidParticles& particles);			
					   
struct PF_ComputePressureData
{
	const btFluidParametersGlobal& m_globalParameters; 
	const btAlignedObjectArray<int>& m_gridCellGroup;
	const btFluidSortingGrid& m_grid;
	btFluidParticles& m_particles; 
	btFluidSolverSph::SphParticles& m_sphData;
	
	PF_ComputePressureData(const btFluidParametersGlobal& FG, const btAlignedObjectArray<int>& gridCellGroup,
							const btFluidSortingGrid& grid, btFluidParticles& particles, btFluidSolverSph::SphParticles& sphData) 
	: m_globalParameters(FG), m_gridCellGroup(gridCellGroup), m_grid(grid), m_particles(particles), m_sphData(sphData) {}
};
void PF_ComputePressureFunction(void* parameters, int index)
{
	PF_ComputePressureData* data = static_cast<PF_ComputePressureData*>(parameters);
	
	calculateSumsInCellSymmetric(data->m_globalParameters, data->m_gridCellGroup[index], data->m_grid, data->m_particles, data->m_sphData);
}

struct PF_ComputeForceData
{
	const btFluidParametersGlobal& m_globalParameters; 
	const btScalar m_vterm;
	const btAlignedObjectArray<int>& m_gridCellGroup;
	const btFluidSortingGrid& m_grid;
	btFluidParticles& m_particles;
	btFluidSolverSph::SphParticles& m_sphData;
	
	PF_ComputeForceData(const btFluidParametersGlobal& FG, const btScalar vterm, const btAlignedObjectArray<int>& gridCellGroup, 
						const btFluidSortingGrid& grid, btFluidParticles& particles, btFluidSolverSph::SphParticles& sphData) 
	: m_globalParameters(FG), m_vterm(vterm),  m_gridCellGroup(gridCellGroup), 
	m_grid(grid), m_particles(particles), m_sphData(sphData) {}
};
void PF_ComputeForceFunction(void* parameters, int index)
{
	PF_ComputeForceData* data = static_cast<PF_ComputeForceData*>(parameters);
	
	calculateForcesInCellSymmetric(data->m_globalParameters, data->m_vterm, data->m_gridCellGroup[index], 
									data->m_grid, data->m_particles, data->m_sphData);
}

struct PF_IntegrateData
{
	const btFluidParametersGlobal& m_globalParameters; 
	const btFluidParametersLocal& m_localParameters; 
	const btScalar m_speedLimitSquared;
	const btScalar m_sphRadiusSquared;
	btFluidParticles& m_particles;
	
	PF_IntegrateData(const btFluidParametersGlobal& FG, const btFluidParametersLocal& FL, const btScalar speedLimitSquared,
					const btScalar R2, btFluidParticles& particles) 
	: m_globalParameters(FG), m_localParameters(FL), m_speedLimitSquared(speedLimitSquared),
	  m_sphRadiusSquared(R2), m_particles(particles) {}
};
void PF_IntegrateFunction(void* parameters, int index)
{
	PF_IntegrateData* data = static_cast<PF_IntegrateData*>(parameters);
	
	integrateParticle(data->m_globalParameters, data->m_localParameters, data->m_speedLimitSquared, 
					  data->m_sphRadiusSquared, index, data->m_particles);
}

#endif //FLUIDS_MULTITHREADED_ENABLED


inline void resolveAabbCollision(btScalar stiff, btScalar damp, const btVector3& vel_eval,
								 btVector3* acceleration, const btVector3& normal, btScalar penetrationDepth)
{
	const btScalar COLLISION_EPSILON = btScalar(0.00001);

	if(penetrationDepth > COLLISION_EPSILON)
	{
		btScalar adj = stiff * penetrationDepth - damp * normal.dot(vel_eval);
		
		btVector3 collisionAcceleration = normal;
		collisionAcceleration *= adj;	

		*acceleration += collisionAcceleration;
	}
}
inline void resolveAabbCollision_impulse(const btVector3& velocity, const btVector3& normal, btScalar distance, 
										 btVector3* impulse, btVector3* positionProjection)
{
	
	const btScalar COLLISION_EPSILON = btScalar(0.00001);
	
	if(distance < -COLLISION_EPSILON)	//Negative distance indicates penetration
	{
		const btScalar RESTITUTION = btScalar(0.0);	///<Fraction of reflected velocity(bounciness); [0.0, 1.0]; higher values more unstable.
		const btScalar FRICTION = btScalar(0.0);	///<Fraction of tangential velocity removed per frame; [0.0, 1.0]; higher values more unstable.
		const btScalar PROJECTION_FRACTION = btScalar(0.005); ///<Fraction of penetration removed per frame; [0.0, 1.0]; higher values more unstable.
	
		btScalar penetratingMagnitude = velocity.dot(-normal);
		if( penetratingMagnitude < btScalar(0.0) ) penetratingMagnitude = btScalar(0.0);
		
		btVector3 penetratingVelocity = -normal * penetratingMagnitude;
		btVector3 tangentialVelocity = velocity - penetratingVelocity;

		//When calculating restitution, include only the velocity omitted in 
		//the collision(penetration) to avoid introducing energy into the system.
		//
		//Equation 4.58 - "Lagrangian Fluid Dynamics Using Smoothed Particle Hydrodynamics". M. Kelager. 9 January 2006.
		//const btScalar TIME_STEP = btScalar(0.003);
		//*impulse += -penetratingVelocity * ( btScalar(1.0) + RESTITUTION *( -distance / (TIME_STEP * velocity.length()) ) );
		
		*impulse += -penetratingVelocity * ( btScalar(1.0) + RESTITUTION );
		*impulse += -tangentialVelocity * FRICTION;
		
		*positionProjection += normal * (-distance * PROJECTION_FRACTION);
	}
}

void integrateParticle(const btFluidParametersGlobal& FG, const btFluidParametersLocal& FL,
					   btScalar speedLimitSquared, btScalar R2, int particleIndex, btFluidParticles& particles)
{		
	const btScalar ss = FG.m_simulationScale;
	
	const btScalar stiff = FL.m_boundaryStiff;
	const btScalar damp = FL.m_boundaryDamp;
	
	const btVector3& min = FL.m_volumeMin;
	const btVector3& max = FL.m_volumeMax;
	
	int i = particleIndex;
	
	//Compute Acceleration
	btVector3 accel = FL.m_gravity;

	//Apply acceleration to keep particles in the particlesystem's AABB
	const bool FORCE_BOUNDARY = true;	//	if true, ensure that (IMPULSE_BOUNDARY == false)
	if(FORCE_BOUNDARY)
	{
		resolveAabbCollision( stiff, damp, particles.m_vel_eval[i], &accel, 
								btVector3( 1.0, 0.0, 0.0), R2 - ( particles.m_pos[i].x() - min.x() )*ss );
		resolveAabbCollision( stiff, damp, particles.m_vel_eval[i], &accel, 
								btVector3(-1.0, 0.0, 0.0), R2 - ( max.x() - particles.m_pos[i].x() )*ss );
		resolveAabbCollision( stiff, damp, particles.m_vel_eval[i], &accel,
								btVector3(0.0,  1.0, 0.0), R2 - ( particles.m_pos[i].y() - min.y() )*ss );
		resolveAabbCollision( stiff, damp, particles.m_vel_eval[i], &accel, 
								btVector3(0.0, -1.0, 0.0), R2 - ( max.y() - particles.m_pos[i].y() )*ss );
		resolveAabbCollision( stiff, damp, particles.m_vel_eval[i], &accel, 
								btVector3(0.0, 0.0,  1.0), R2 - ( particles.m_pos[i].z() - min.z() )*ss );
		resolveAabbCollision( stiff, damp, particles.m_vel_eval[i], &accel, 
								btVector3(0.0, 0.0, -1.0), R2 - ( max.z() - particles.m_pos[i].z() )*ss );
	}
	
	//Apply forces
	accel += particles.m_accumulatedForce[i] / FL.m_particleMass;
	particles.m_accumulatedForce[i].setValue(0, 0, 0);

	//Integrate velocity
	btVector3 vnext = particles.m_vel[i] + accel * FG.m_timeStep;			//v(t+1/2) = v(t-1/2) + a(t) dt	
	particles.m_vel_eval[i] = (particles.m_vel[i] + vnext) * btScalar(0.5);	//v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
	particles.m_vel[i] = vnext;
	
	//Apply impulses
	btVector3 positionProjection(0.f, 0.f, 0.f);
	const bool IMPULSE_BOUNDARY = false;		//	if true, ensure that (FORCE_BOUNDARY == false)
	if(IMPULSE_BOUNDARY)
	{
		btVector3 impulse(0.f, 0.f, 0.f);
		resolveAabbCollision_impulse( particles.m_vel[i], btVector3(1.0, 0.0, 0.0), 
									 ( particles.m_pos[i].x() - min.x() )*ss - R2, &impulse, &positionProjection );
		resolveAabbCollision_impulse( particles.m_vel[i], btVector3(-1.0, 0.0, 0.0), 
									 ( max.x() - particles.m_pos[i].x() )*ss - R2, &impulse, &positionProjection );
		resolveAabbCollision_impulse( particles.m_vel[i], btVector3(0.0,  1.0, 0.0), 
									 ( particles.m_pos[i].y() - min.y() )*ss - R2, &impulse, &positionProjection );
		resolveAabbCollision_impulse( particles.m_vel[i], btVector3(0.0, -1.0, 0.0), 
									 ( max.y() - particles.m_pos[i].y() )*ss - R2, &impulse, &positionProjection );
		resolveAabbCollision_impulse( particles.m_vel[i], btVector3(0.0, 0.0,  1.0), 
									 ( particles.m_pos[i].z() - min.z() )*ss - R2, &impulse, &positionProjection );
		resolveAabbCollision_impulse( particles.m_vel[i], btVector3(0.0, 0.0, -1.0), 
									 ( max.z() - particles.m_pos[i].z() )*ss - R2, &impulse, &positionProjection );

		btVector3 vnext = particles.m_vel[i] + impulse;
		particles.m_vel_eval[i] = (particles.m_vel[i] + vnext) * btScalar(0.5);
		particles.m_vel[i] = vnext;
	}
	
	//Integrate position
	//p(t+1) = p(t) + v(t+1/2)*dt
	
	if(IMPULSE_BOUNDARY) particles.m_pos[i] += particles.m_vel[i]*(FG.m_timeStep / ss) + positionProjection;
	else particles.m_pos[i] += particles.m_vel[i]*(FG.m_timeStep / ss);
	
}
void btFluidSolver::integrate(const btFluidParametersGlobal& FG, const btFluidParametersLocal& FL, btFluidParticles& particles)
{
	BT_PROFILE("btFluidSolver::integrate()");
	
	btScalar speedLimitSquared = FG.m_speedLimit*FG.m_speedLimit;
	btScalar R2 = btScalar(2.0) * FL.m_particleRadius * FG.m_simulationScale;
	
#ifdef FLUIDS_MULTITHREADED_ENABLED
	PF_IntegrateData IntegrateData(FG, FL, speedLimitSquared, R2, particles);
	parallelFor.execute( PF_IntegrateFunction, &IntegrateData, 0, particles.size() - 1, 256 );
#else
	for(int i = 0; i < particles.size(); ++i)
		integrateParticle(FG, FL, speedLimitSquared, R2, i, particles);
#endif
}

		
void calculateSumsInCellSymmetric(const btFluidParametersGlobal& FG, int gridCellIndex, const btFluidSortingGrid& grid, 
								btFluidParticles& particles, btFluidSolverSph::SphParticles& sphData)
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
						btScalar c_cubed = c * c * c;
						sphData.m_invDensity[i] += c_cubed;
						sphData.m_invDensity[n] += c_cubed;
						
						btScalar distance = btSqrt(distanceSquared);
						if( !sphData.m_neighborTable[i].isFilled() ) sphData.m_neighborTable[i].addNeighbor(n, distance);
						else if( !sphData.m_neighborTable[n].isFilled() ) sphData.m_neighborTable[n].addNeighbor(i, distance);
						else break;
					}
				}
			}
		}
	}
}

void btFluidSolverSph::sphComputePressure(const btFluidParametersGlobal& FG, btFluidSph* fluid, btFluidSolverSph::SphParticles& sphData)
{
	BT_PROFILE("btFluidSolverSph::sphComputePressure()");
	
	const int numParticles = fluid->numParticles();
	
	const btFluidParametersLocal& FL = fluid->getLocalParameters();
	const btFluidSortingGrid& grid = fluid->getGrid();
	btFluidParticles& particles = fluid->internalGetParticles();
	
	{
		BT_PROFILE("sphComputePressure() - reset sums, clear table");
		
		for(int i = 0; i < numParticles; ++i) sphData.m_invDensity[i] = btScalar(0.0);
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
			btScalar density = sphData.m_invDensity[i] * FL.m_particleMass * FG.m_poly6KernCoeff;	
			sphData.m_pressure[i] = (density - FL.m_restDensity) * FL.m_stiffness;
			sphData.m_invDensity[i] = btScalar(1.0) / density;
		}
	}
}

void computeForceNeighborTableSymmetric(const btFluidParametersGlobal& FG, const btScalar vterm, int particleIndex, 
										btFluidParticles& particles, btFluidSolverSph::SphParticles& sphData)
{
	int i = particleIndex;
	
	for(int j = 0; j < sphData.m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = sphData.m_neighborTable[i].getNeighborIndex(j);
		
		btVector3 distance = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
		
		btScalar c = FG.m_sphSmoothRadius - sphData.m_neighborTable[i].getDistance(j);
		btScalar pterm = -0.5f * c * FG.m_spikyKernGradCoeff 
					 * ( sphData.m_pressure[i] + sphData.m_pressure[n]) / sphData.m_neighborTable[i].getDistance(j);
		btScalar dterm = c * sphData.m_invDensity[i] * sphData.m_invDensity[n];

		btVector3 forceAdded( (pterm * distance.x() + vterm * (particles.m_vel_eval[n].x() - particles.m_vel_eval[i].x())) * dterm,
							  (pterm * distance.y() + vterm * (particles.m_vel_eval[n].y() - particles.m_vel_eval[i].y())) * dterm,
							  (pterm * distance.z() + vterm * (particles.m_vel_eval[n].z() - particles.m_vel_eval[i].z())) * dterm );
		
		sphData.m_sphForce[i] += forceAdded;
		sphData.m_sphForce[n] += -forceAdded;
	}
}
void calculateForcesInCellSymmetric(const btFluidParametersGlobal& FG, const btScalar vterm,
									int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
									btFluidSolverSph::SphParticles& sphData)
{
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
	{
		computeForceNeighborTableSymmetric(FG, vterm, i, particles, sphData);
	}
}
void btFluidSolverSph::sphComputeForce(const btFluidParametersGlobal& FG, btFluidSph* fluid, btFluidSolverSph::SphParticles& sphData)
{
	BT_PROFILE("btFluidSolverSph::sphComputeForce()");
	
	const btFluidParametersLocal& FL = fluid->getLocalParameters();
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
	
	for(int i = 0; i < particles.size(); ++i)sphData.m_sphForce[i] *= FL.m_particleMass;
}
