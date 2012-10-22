/* FluidSolver.cpp

	ZLib license
	This software is provided 'as-is', without any express or implied
	warranty. In no event will the authors be held liable for any damages
	arising from the use of this software.
	
	Permission is granted to anyone to use this software for any purpose,
	including commercial applications, and to alter it and redistribute it
	freely, subject to the following restrictions:
	
	1. The origin of this software must not be misrepresented; you must not
	   claim that you wrote the original software. If you use this software
	   in a product, an acknowledgment in the product documentation would be
	   appreciated but is not required.
	2. Altered source versions must be plainly marked as such, and must not be
	   misrepresented as being the original software.
	3. This notice may not be removed or altered from any source distribution.
*/

#include "FluidSolver.h"

#include "FluidSortingGrid.h"

//Link to e.g. 'BulletMultiThreaded.lib' if enabling this
//(must build Bullet with CMake to get BulletMultiThreaded library)
//#define FLUIDS_MULTITHREADED_ENABLED
#ifdef FLUIDS_MULTITHREADED_ENABLED
#include "BulletMultiThreadedSupport.h"

const unsigned int NUM_THREADS = 4;
ParallelFor parallelFor("parallelForThreads", NUM_THREADS);

void calculatePressuresInCellSymmetric(const FluidParametersGlobal &FG, const btScalar gridSearchRadius, int gridCellIndex, 
									   FluidSortingGrid *tempGrid, FluidParticles *fluids, btAlignedObjectArray<btScalar> *sums);	
void calculateForcesInCellSymmetric(const FluidParametersGlobal &FG, const btScalar vterm,
									int gridCellIndex, const FluidSortingGrid *grid, FluidParticles *fluids);
void integrateParticle(const FluidParametersGlobal &FG, const FluidParametersLocal &FL,
					   btScalar speedLimitSquared, btScalar R2,  bool isPlaneGravityEnabled, 
					   int particleIndex, FluidParticles *fluids);			
					   
struct PF_ComputePressureData
{
	const FluidParametersGlobal &m_globalParameters; 
	const btScalar m_gridSearchRadius;
	const btAlignedObjectArray<int> &m_gridCellGroup;
	FluidSortingGrid *m_tempGrid;
	FluidParticles *m_particles; 
	btAlignedObjectArray<btScalar> *m_sums;
	
	PF_ComputePressureData(const FluidParametersGlobal &FG, const btScalar gridSearchRadius,
						   const btAlignedObjectArray<int> &gridCellGroup, FluidSortingGrid *tempGrid,
						   FluidParticles *particles, btAlignedObjectArray<btScalar> *sums) 
	: m_globalParameters(FG), m_gridSearchRadius(gridSearchRadius),
	  m_gridCellGroup(gridCellGroup),m_tempGrid(tempGrid), 
	  m_particles(particles), m_sums(sums) {}
};
void PF_ComputePressureFunction(void *parameters, int index)
{
	PF_ComputePressureData *data = static_cast<PF_ComputePressureData*>(parameters);
	
	calculatePressuresInCellSymmetric(data->m_globalParameters, data->m_gridSearchRadius, 
									  data->m_gridCellGroup[index], data->m_tempGrid, data->m_particles, data->m_sums);
}

struct PF_ComputeForceData
{
	const FluidParametersGlobal &m_globalParameters; 
	const btScalar m_vterm;
	const btAlignedObjectArray<int> &m_gridCellGroup;
	const FluidSortingGrid *m_grid;
	FluidParticles *m_particles;
	
	PF_ComputeForceData(const FluidParametersGlobal &FG, const btScalar vterm,
						const btAlignedObjectArray<int> &gridCellGroup, const FluidSortingGrid *grid, FluidParticles *particles) 
	: m_globalParameters(FG), m_vterm(vterm),
	  m_gridCellGroup(gridCellGroup), m_grid(grid), m_particles(particles) {}
};
void PF_ComputeForceFunction(void *parameters, int index)
{
	PF_ComputeForceData *data = static_cast<PF_ComputeForceData*>(parameters);
	
	calculateForcesInCellSymmetric(data->m_globalParameters, data->m_vterm, 
								   data->m_gridCellGroup[index], data->m_grid, data->m_particles);
}

struct PF_IntegrateData
{
	const FluidParametersGlobal &m_globalParameters; 
	const FluidParametersLocal &m_localParameters; 
	const btScalar m_speedLimitSquared;
	const btScalar m_sphRadiusSquared;
	const bool m_isPlaneGravityEnabled;
	FluidParticles *m_particles;
	
	PF_IntegrateData(const FluidParametersGlobal &FG, const FluidParametersLocal &FL, const btScalar speedLimitSquared,
				   const btScalar R2, const bool isPlaneGravityEnabled, FluidParticles *particles) 
	: m_globalParameters(FG), m_localParameters(FL), m_speedLimitSquared(speedLimitSquared),
	  m_sphRadiusSquared(R2), m_isPlaneGravityEnabled(isPlaneGravityEnabled), m_particles(particles) {}
};
void PF_IntegrateFunction(void *parameters, int index)
{
	PF_IntegrateData *data = static_cast<PF_IntegrateData*>(parameters);
	
	integrateParticle(data->m_globalParameters, data->m_localParameters, data->m_speedLimitSquared, 
					  data->m_sphRadiusSquared, data->m_isPlaneGravityEnabled, index, data->m_particles);
}

#endif //FLUIDS_MULTITHREADED_ENABLED


void FluidSolverGridNeighbor::sphComputePressure(const FluidParametersGlobal &FG, FluidSph *fluid)
{
	BT_PROFILE("sphComputePressure()");
	
	const FluidParametersLocal &FL = fluid->getLocalParameters();
	FluidParticles &particles = fluid->internalGetFluidParticles();
	FluidSortingGrid &grid = fluid->internalGetGrid();
	
	btScalar radius = FG.m_sphSmoothRadius / FG.m_simulationScale;

	for(int i = 0; i < fluid->numParticles(); ++i)
	{
		btScalar sum = 0.0;	
		particles.m_neighborTable[i].clear();

		FluidSortingGrid::FoundCells foundCells;
		grid.findCells(particles.m_pos[i], radius, &foundCells);
		
		for(int cell = 0; cell < FluidSortingGrid::NUM_FOUND_CELLS; cell++) 
		{
			FluidGridIterator &FI = foundCells.m_iterators[cell];
			
			for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
			{
				if(i == n) continue;
				
				btVector3 distance = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
				btScalar distanceSquared = distance.length2();
				
				if(FG.m_sphRadiusSquared > distanceSquared) 
				{
					btScalar c = FG.m_sphRadiusSquared - distanceSquared;
					sum += c * c * c;
					
					if( !particles.m_neighborTable[i].isFilled() ) particles.m_neighborTable[i].addNeighbor( n, btSqrt(distanceSquared) );
				}
			}
		}
		
		btScalar density = sum * FL.m_particleMass * FG.m_poly6KernCoeff;	
		particles.m_pressure[i] = (density - FL.m_restDensity) * FL.m_stiffness;
		particles.m_invDensity[i] = 1.0f / density;
	}
}

void computeForceNeighborTable(const FluidParametersGlobal &FG, const btScalar vterm, int particleIndex, FluidParticles *fluids)
{
	int i = particleIndex;

	btVector3 force(0, 0, 0);
	for(int j = 0; j < fluids->m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = fluids->m_neighborTable[i].getNeighborIndex(j);
		
		btVector3 distance = (fluids->m_pos[i] - fluids->m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
		
		btScalar c = FG.m_sphSmoothRadius - fluids->m_neighborTable[i].getDistance(j);
		btScalar pterm = -0.5f * c * FG.m_spikyKernGradCoeff 
					 * ( fluids->m_pressure[i] + fluids->m_pressure[n]) / fluids->m_neighborTable[i].getDistance(j);
		btScalar dterm = c * fluids->m_invDensity[i] * fluids->m_invDensity[n];

		btVector3 forceAdded( (pterm * distance.x() + vterm * (fluids->m_vel_eval[n].x() - fluids->m_vel_eval[i].x())) * dterm,
							  (pterm * distance.y() + vterm * (fluids->m_vel_eval[n].y() - fluids->m_vel_eval[i].y())) * dterm,
							  (pterm * distance.z() + vterm * (fluids->m_vel_eval[n].z() - fluids->m_vel_eval[i].z())) * dterm );
		force += forceAdded;
	}
	
	fluids->m_sph_force[i] = force;
}

void FluidSolverGridNeighbor::sphComputeForce(const FluidParametersGlobal &FG, FluidSph *fluid)
{
	BT_PROFILE("sphComputeForce()");
	
	const FluidParametersLocal &FL = fluid->getLocalParameters();
	FluidParticles &fluids = fluid->internalGetFluidParticles();
	btScalar vterm = FG.m_viscosityKernLapCoeff * FL.m_viscosity;
	
	for(int i = 0; i < fluids.size(); ++i)
		computeForceNeighborTable(FG, vterm, i, &fluids);
}

inline void resolveAabbCollision(btScalar stiff, btScalar damp, const btVector3 &vel_eval,
								 btVector3 *acceleration, const btVector3 &normal, btScalar penetrationDepth)
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
inline void resolveAabbCollision_impulse(const btVector3 &velocity, const btVector3 &normal, btScalar distance, 
										 btVector3 *impulse, btVector3 *positionProjection)
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

void integrateParticle(const FluidParametersGlobal &FG, const FluidParametersLocal &FL,
					   btScalar speedLimitSquared, btScalar R2, 
					   bool isPlaneGravityEnabled, int particleIndex, FluidParticles *fluids)
{		
	const btScalar ss = FG.m_simulationScale;
	
	const btScalar stiff = FL.m_boundaryStiff;
	const btScalar damp = FL.m_boundaryDamp;
	
	const btVector3 &min = FL.m_volumeMin;
	const btVector3 &max = FL.m_volumeMax;
	
	int i = particleIndex;
	
	//Compute Acceleration		
	btVector3 accel = fluids->m_sph_force[i];
	accel *= FL.m_particleMass;

	//Limit speed
	btScalar speedSquared = accel.length2();
	if(speedSquared > speedLimitSquared) accel *= FG.m_speedLimit / btSqrt(speedSquared);

	//Apply acceleration to keep particles in the FluidSystem's AABB
	const bool FORCE_BOUNDARY = true;	//	if true, ensure that (IMPULSE_BOUNDARY == false)
	if(FORCE_BOUNDARY)
	{
		resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3( 1.0, 0.0, 0.0), R2 - ( fluids->m_pos[i].x() - min.x() )*ss );
		resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3(-1.0, 0.0, 0.0), R2 - ( max.x() - fluids->m_pos[i].x() )*ss );
		resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3(0.0,  1.0, 0.0), R2 - ( fluids->m_pos[i].y() - min.y() )*ss );
		resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3(0.0, -1.0, 0.0), R2 - ( max.y() - fluids->m_pos[i].y() )*ss );
		resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3(0.0, 0.0,  1.0), R2 - ( fluids->m_pos[i].z() - min.z() )*ss );
		resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3(0.0, 0.0, -1.0), R2 - ( max.z() - fluids->m_pos[i].z() )*ss );
	}
	
	//Plane gravity
	if(isPlaneGravityEnabled) accel += FG.m_planeGravity;

	//Point gravity
	if(FG.m_pointGravity > 0.0) 
	{
		btVector3 norm = fluids->m_pos[i] - FG.m_pointGravityPosition;
		norm.normalize();
		norm *= FG.m_pointGravity;
		accel -= norm;
	}
	
	//Apply external forces
	accel += fluids->m_externalAcceleration[i];
	fluids->m_externalAcceleration[i].setValue(0, 0, 0);

	//Integrate velocity
	btVector3 vnext = fluids->m_vel[i] + accel * FG.m_timeStep;			//v(t+1/2) = v(t-1/2) + a(t) dt	
	fluids->m_vel_eval[i] = (fluids->m_vel[i] + vnext) * btScalar(0.5);	//v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
	fluids->m_vel[i] = vnext;
	
	//Apply impulses
	btVector3 positionProjection(0.f, 0.f, 0.f);
	const bool IMPULSE_BOUNDARY = false;		//	if true, ensure that (FORCE_BOUNDARY == false)
	if(IMPULSE_BOUNDARY)
	{
		btVector3 impulse(0.f, 0.f, 0.f);
		resolveAabbCollision_impulse( fluids->m_vel[i], btVector3(1.0, 0.0, 0.0), 
									 ( fluids->m_pos[i].x() - min.x() )*ss - R2, &impulse, &positionProjection );
		resolveAabbCollision_impulse( fluids->m_vel[i], btVector3(-1.0, 0.0, 0.0), 
									 ( max.x() - fluids->m_pos[i].x() )*ss - R2, &impulse, &positionProjection );
		resolveAabbCollision_impulse( fluids->m_vel[i], btVector3(0.0,  1.0, 0.0), 
									 ( fluids->m_pos[i].y() - min.y() )*ss - R2, &impulse, &positionProjection );
		resolveAabbCollision_impulse( fluids->m_vel[i], btVector3(0.0, -1.0, 0.0), 
									 ( max.y() - fluids->m_pos[i].y() )*ss - R2, &impulse, &positionProjection );
		resolveAabbCollision_impulse( fluids->m_vel[i], btVector3(0.0, 0.0,  1.0), 
									 ( fluids->m_pos[i].z() - min.z() )*ss - R2, &impulse, &positionProjection );
		resolveAabbCollision_impulse( fluids->m_vel[i], btVector3(0.0, 0.0, -1.0), 
									 ( max.z() - fluids->m_pos[i].z() )*ss - R2, &impulse, &positionProjection );

		btVector3 vnext = fluids->m_vel[i] + impulse;
		fluids->m_vel_eval[i] = (fluids->m_vel[i] + vnext) * btScalar(0.5);
		fluids->m_vel[i] = vnext;
	}
	
	//Integrate position
	//p(t+1) = p(t) + v(t+1/2)*dt
	
	if(IMPULSE_BOUNDARY) fluids->m_pos[i] += fluids->m_vel[i]*(FG.m_timeStep / ss) + positionProjection;
	else fluids->m_pos[i] += fluids->m_vel[i]*(FG.m_timeStep / ss);
	
}
void FluidSolver::integrate(const FluidParametersGlobal &FG, const FluidParametersLocal &FL, FluidParticles *fluids)
{
	BT_PROFILE("FluidSolver::integrate()");
	
	btScalar speedLimitSquared = FG.m_speedLimit*FG.m_speedLimit;
	btScalar R2 = 2.0f * FG.m_particleRadius;
	
	bool isPlaneGravityEnabled = !FG.m_planeGravity.isZero();
	
#ifdef FLUIDS_MULTITHREADED_ENABLED
	PF_IntegrateData IntegrateData(FG, FL, speedLimitSquared, R2, isPlaneGravityEnabled, fluids);
	parallelFor.execute( PF_IntegrateFunction, &IntegrateData, 0, fluids->size() - 1, 256 );
#else
	for(int i = 0; i < fluids->size(); ++i)
		integrateParticle(FG, FL, speedLimitSquared, R2, isPlaneGravityEnabled, i, fluids);
#endif
}

void calculatePressuresInCellSymmetric(const FluidParametersGlobal &FG, const btScalar gridSearchRadius, int gridCellIndex, 
									   FluidSortingGrid *tempGrid, FluidParticles *fluids, btAlignedObjectArray<btScalar> *sums)
{
#ifdef GRID_CELL_SIZE_2R		//#defined in FluidSortingGrid.h
	FluidGridIterator currentCell = tempGrid->getGridCell(gridCellIndex);
	for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
	{
		//Remove particle, with index i, from grid cell
		tempGrid->internalRemoveFirstParticle(gridCellIndex);
	
		FluidSortingGrid::FoundCells foundCells;
		tempGrid->findCells(fluids->m_pos[i], gridSearchRadius, &foundCells);
		for(int cell = 0; cell < FluidSortingGrid::NUM_FOUND_CELLS; cell++) 
		{
			FluidGridIterator &FI = foundCells.m_iterators[cell];
			
			for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
			{
				//Simulation-scale distance
				btVector3 distanceAsVector = (fluids->m_pos[i] - fluids->m_pos[n]) * FG.m_simulationScale;		
				btScalar distanceSquared = distanceAsVector.length2();
				
				if(FG.m_sphRadiusSquared > distanceSquared) 
				{
					btScalar c = FG.m_sphRadiusSquared - distanceSquared;
					btScalar c_cubed = c * c * c;
					(*sums)[i] += c_cubed;
					(*sums)[n] += c_cubed;
					
					btScalar distance = btSqrt(distanceSquared);
					if( !fluids->m_neighborTable[i].isFilled() ) fluids->m_neighborTable[i].addNeighbor(n, distance);
					else if( !fluids->m_neighborTable[n].isFilled() ) fluids->m_neighborTable[n].addNeighbor(i, distance);
					else break;
				}
			}
		}
	}
#else
	FluidGridIterator currentCell = tempGrid->getGridCell(gridCellIndex);
	if(currentCell.m_firstIndex <= currentCell.m_lastIndex)	//if cell is not empty
	{
		FluidSortingGrid::FoundCells foundCells;
		tempGrid->findCells(fluids->m_pos[currentCell.m_firstIndex], gridSearchRadius, &foundCells);
		//tempGrid->internalRemoveAllParticles(gridCellIndex);
		
		for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
		{
			//Remove particle, with index i, from grid cell
			tempGrid->internalRemoveFirstParticle(gridCellIndex);	//tempGrid cell
			++foundCells.m_iterators[0].m_firstIndex;				//local cell; currentCell == foundCells.m_iterators[0]
			
			for(int cell = 0; cell < FluidSortingGrid::NUM_FOUND_CELLS; cell++) 
			{
				FluidGridIterator &FI = foundCells.m_iterators[cell];
				
				for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
				{
					//Simulation-scale distance
					btVector3 distanceAsVector = (fluids->m_pos[i] - fluids->m_pos[n]) * FG.m_simulationScale;		
					btScalar distanceSquared = distanceAsVector.length2();
					
					if(FG.m_sphRadiusSquared > distanceSquared) 
					{
						btScalar c = FG.m_sphRadiusSquared - distanceSquared;
						btScalar c_cubed = c * c * c;
						(*sums)[i] += c_cubed;
						(*sums)[n] += c_cubed;
						
						btScalar distance = btSqrt(distanceSquared);
						if( !fluids->m_neighborTable[i].isFilled() ) fluids->m_neighborTable[i].addNeighbor(n, distance);
						else if( !fluids->m_neighborTable[n].isFilled() ) fluids->m_neighborTable[n].addNeighbor(i, distance);
						else break;
					}
				}
			}
		}
	}
#endif
}

//Remove particles from grid cells as their interactions are calculated
void FluidSolverReducedGridNeighbor::sphComputePressureGridReduce(const FluidParametersGlobal &FG, FluidSph *fluid)
{
	BT_PROFILE("sphComputePressureGridReduce()");
	
	const int numParticles = fluid->numParticles();
	
	const FluidParametersLocal &FL = fluid->getLocalParameters();
	FluidParticles &particles = fluid->internalGetFluidParticles();
	FluidSortingGrid &grid = fluid->internalGetGrid();
	
	static FluidSortingGrid tempGrid;
	static btAlignedObjectArray<btScalar> sums;
	{
		BT_PROFILE("sphComputePressureGridReduce() - copy grid, reset sums, clear table");
		
		tempGrid = grid;
		
		sums.resize(numParticles);
		for(int i = 0; i < numParticles; ++i) sums[i] = btScalar(0.0);
		for(int i = 0; i < numParticles; ++i) particles.m_neighborTable[i].clear();
	}
	
	{
		BT_PROFILE("sphComputePressureGridReduce() - compute sums");
		btScalar gridSearchRadius = FG.m_sphSmoothRadius / FG.m_simulationScale;
		
		for(int group = 0; group < FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS; ++group)
		{
			const btAlignedObjectArray<int> &currentGroup = grid.internalGetCellProcessingGroup(group);
			if( !currentGroup.size() ) continue;
		
#ifdef FLUIDS_MULTITHREADED_ENABLED
			PF_ComputePressureData PressureData(FG, gridSearchRadius, currentGroup, &tempGrid, &particles, &sums);
			parallelFor.execute( PF_ComputePressureFunction, &PressureData, 0, currentGroup.size() - 1, 1 );
#else
			for(int cell = 0; cell < currentGroup.size(); ++cell)
				calculatePressuresInCellSymmetric(FG, gridSearchRadius, currentGroup[cell], &tempGrid, &particles, &sums);
#endif		
		}
	}
	
	{
		BT_PROFILE("sphComputePressureGridReduce() - compute pressure/density");
	
		for(int i = 0; i < numParticles; ++i)
		{
			btScalar density = sums[i] * FL.m_particleMass * FG.m_poly6KernCoeff;	
			particles.m_pressure[i] = (density - FL.m_restDensity) * FL.m_stiffness;
			particles.m_invDensity[i] = btScalar(1.0) / density;
		}
	}
}

void computeForceNeighborTableSymmetric(const FluidParametersGlobal &FG, const btScalar vterm, int particleIndex, FluidParticles *fluids)
{
	int i = particleIndex;

	for(int j = 0; j < fluids->m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = fluids->m_neighborTable[i].getNeighborIndex(j);
		
		btVector3 distance = (fluids->m_pos[i] - fluids->m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
		
		btScalar c = FG.m_sphSmoothRadius - fluids->m_neighborTable[i].getDistance(j);
		btScalar pterm = -0.5f * c * FG.m_spikyKernGradCoeff 
					 * ( fluids->m_pressure[i] + fluids->m_pressure[n]) / fluids->m_neighborTable[i].getDistance(j);
		btScalar dterm = c * fluids->m_invDensity[i] * fluids->m_invDensity[n];

		btVector3 forceAdded( (pterm * distance.x() + vterm * (fluids->m_vel_eval[n].x() - fluids->m_vel_eval[i].x())) * dterm,
							  (pterm * distance.y() + vterm * (fluids->m_vel_eval[n].y() - fluids->m_vel_eval[i].y())) * dterm,
							  (pterm * distance.z() + vterm * (fluids->m_vel_eval[n].z() - fluids->m_vel_eval[i].z())) * dterm );
		
		fluids->m_sph_force[i] += forceAdded;
		fluids->m_sph_force[n] += -forceAdded;
	}
}
void calculateForcesInCellSymmetric(const FluidParametersGlobal &FG, const btScalar vterm,
									int gridCellIndex, const FluidSortingGrid *grid, FluidParticles *fluids)
{
	FluidGridIterator currentCell = grid->getGridCell(gridCellIndex);
	for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
	{
		computeForceNeighborTableSymmetric(FG, vterm, i, fluids);
	}
}
void FluidSolverReducedGridNeighbor::sphComputeForceReduce(const FluidParametersGlobal &FG, FluidSph *fluid)
{
	BT_PROFILE("FluidSolverReducedGridNeighbor::sphComputeForceReduce()");
	
	const FluidParametersLocal &FL = fluid->getLocalParameters();
	FluidParticles &particles = fluid->internalGetFluidParticles();
	FluidSortingGrid &grid = fluid->internalGetGrid();
	
	btScalar vterm = FG.m_viscosityKernLapCoeff * FL.m_viscosity;
	
	for(int i = 0; i < particles.size(); ++i)particles.m_sph_force[i].setValue(0, 0, 0);
	
	for(int group = 0; group < FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS; ++group)
	{
		const btAlignedObjectArray<int> &currentGroup = grid.internalGetCellProcessingGroup(group);
		if( !currentGroup.size() ) continue;
		
#ifdef FLUIDS_MULTITHREADED_ENABLED
		PF_ComputeForceData ForceData(FG, vterm, currentGroup, &grid, &particles);
		parallelFor.execute( PF_ComputeForceFunction, &ForceData, 0, currentGroup.size() - 1, 1 );
#else
		for(int cell = 0; cell < currentGroup.size(); ++cell)
			calculateForcesInCellSymmetric(FG, vterm, currentGroup[cell], &grid, &particles);
#endif		
	}
}
