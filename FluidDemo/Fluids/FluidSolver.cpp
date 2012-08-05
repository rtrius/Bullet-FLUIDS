/** FluidSolver.cpp

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

#include "grid.h"
#include "hashgrid.h"

//Link to e.g. 'BulletMultiThreaded.lib' if enabling this
//(must build Bullet with CMake to get BulletMultiThreaded library)
//#define FLUIDS_MULTITHREADED_ENABLED
#ifdef FLUIDS_MULTITHREADED_ENABLED
#include "BulletMultiThreadedSupport.h"

const unsigned int NUM_THREADS = 4;
ParallelFor parallelFor("parallelForThreads", NUM_THREADS);

void calculatePressuresInCellSymmetric(const FluidParametersGlobal &FG, const btScalar gridSearchRadius, int gridCellIndex, 
									   FluidGrid *tempGrid, FluidParticles *fluids, btAlignedObjectArray<btScalar> *sums);	
void calculateForcesInCellSymmetric(const FluidParametersGlobal &FG, const btScalar vterm, bool isLinkedList, 
									int gridCellIndex, const FluidGrid *grid, FluidParticles *fluids);
void integrateParticle(const FluidParametersGlobal &FG, const FluidParametersLocal &FL,
					   btScalar speedLimitSquared, btScalar R2,  bool isPlaneGravityEnabled, 
					   int particleIndex, FluidParticles *fluids);			
					   
struct PF_ComputePressureData
{
	const FluidParametersGlobal &m_globalParameters; 
	const btScalar m_gridSearchRadius;
	const btAlignedObjectArray<int> &m_gridCellGroup;
	FluidGrid *m_tempGrid;
	FluidParticles *m_particles; 
	btAlignedObjectArray<btScalar> *m_sums;
	
	PF_ComputePressureData(const FluidParametersGlobal &FG, const btScalar gridSearchRadius,
						   const btAlignedObjectArray<int> &gridCellGroup, FluidGrid *tempGrid,
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
	const bool m_isLinkedList;
	const btAlignedObjectArray<int> &m_gridCellGroup;
	const FluidGrid *m_grid;
	FluidParticles *m_particles;
	
	PF_ComputeForceData(const FluidParametersGlobal &FG, const btScalar vterm, const bool isLinkedList,
						const btAlignedObjectArray<int> &gridCellGroup, const FluidGrid *grid, FluidParticles *particles) 
	: m_globalParameters(FG), m_vterm(vterm), m_isLinkedList(isLinkedList),
	  m_gridCellGroup(gridCellGroup), m_grid(grid), m_particles(particles) {}
};
void PF_ComputeForceFunction(void *parameters, int index)
{
	PF_ComputeForceData *data = static_cast<PF_ComputeForceData*>(parameters);
	
	calculateForcesInCellSymmetric(data->m_globalParameters, data->m_vterm, data->m_isLinkedList, 
								   data->m_gridCellGroup[index], data->m_grid, data->m_particles);
}

struct PF_IntegrateData
{
	const FluidParametersGlobal &m_globalParameters; 
	const FluidParametersLocal &m_localParameters; 
	const btScalar m_speedLimitSquared;
	const btScalar m_R2;
	const bool m_isPlaneGravityEnabled;
	FluidParticles *m_particles;
	
	PF_IntegrateData(const FluidParametersGlobal &FG, const FluidParametersLocal &FL, const btScalar speedLimitSquared,
				   const btScalar R2, const bool isPlaneGravityEnabled, FluidParticles *particles) 
	: m_globalParameters(FG), m_localParameters(FL), m_speedLimitSquared(speedLimitSquared),
	  m_R2(R2), m_isPlaneGravityEnabled(isPlaneGravityEnabled), m_particles(particles) {}
};
void PF_IntegrateFunction(void *parameters, int index)
{
	PF_IntegrateData *data = static_cast<PF_IntegrateData*>(parameters);
	
	integrateParticle(data->m_globalParameters, data->m_localParameters, data->m_speedLimitSquared, 
					  data->m_R2, data->m_isPlaneGravityEnabled, index, data->m_particles);
}

#endif //FLUIDS_MULTITHREADED_ENABLED


void FluidSolverGridNeighbor::sphComputePressure(const FluidParametersGlobal &FG, FluidSph *fluid)
{
	BT_PROFILE("sphComputePressure()");
	
	const FluidParametersLocal &FL = fluid->getLocalParameters();
	FluidParticles &particles = fluid->internalGetFluidParticles();
	FluidGrid *grid = fluid->internalGetGrid();
	
	const bool isLinkedList = (grid->getGridType() == FT_LinkedList);
	btScalar radius = FG.sph_smoothradius / FG.sph_simscale;

	for(int i = 0; i < fluid->numParticles(); ++i)
	{
		btScalar sum = 0.0;	
		particles.m_neighborTable[i].clear();

		FindCellsResult foundCells;
		grid->findCells(particles.m_pos[i], radius, &foundCells);
		
		for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; cell++) 
		{
			FluidGridIterator &FI = foundCells.m_iterators[cell];
			
			for( int n = FI.m_firstIndex; FluidGridIterator::isIndexValid(n, FI.m_lastIndex); 
					 n = FluidGridIterator::getNextIndex(n, isLinkedList, particles.m_nextFluidIndex) )
			{
				if(i == n) continue;
				
				btVector3 distance = (particles.m_pos[i] - particles.m_pos[n]) * FG.sph_simscale;		//Simulation-scale distance
				btScalar distanceSquared = distance.length2();
				
				if(FG.m_R2 > distanceSquared) 
				{
					btScalar c = FG.m_R2 - distanceSquared;
					sum += c * c * c;
					
					if( !particles.m_neighborTable[i].isFilled() ) particles.m_neighborTable[i].addNeighbor( n, btSqrt(distanceSquared) );
				}
			}
		}
		
		btScalar tempDensity = sum * FL.m_particleMass * FG.m_Poly6Kern;	
		particles.m_pressure[i] = (tempDensity - FL.m_restDensity) * FL.m_intstiff;
		particles.m_density[i] = 1.0f / tempDensity;
	}
}

void computeForceNeighborTable(const FluidParametersGlobal &FG, const btScalar vterm, int particleIndex, FluidParticles *fluids)
{
	int i = particleIndex;

	btVector3 force(0, 0, 0);
	for(int j = 0; j < fluids->m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = fluids->m_neighborTable[i].getNeighborIndex(j);
		
		btVector3 distance = (fluids->m_pos[i] - fluids->m_pos[n]) * FG.sph_simscale;		//Simulation-scale distance
		
		btScalar c = FG.sph_smoothradius - fluids->m_neighborTable[i].getDistance(j);
		btScalar pterm = -0.5f * c * FG.m_SpikyKern 
					 * ( fluids->m_pressure[i] + fluids->m_pressure[n]) / fluids->m_neighborTable[i].getDistance(j);
		btScalar dterm = c * fluids->m_density[i] * fluids->m_density[n];

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
	btScalar vterm = FG.m_LapKern * FL.m_viscosity;
	
	for(int i = 0; i < fluids.size(); ++i)
		computeForceNeighborTable(FG, vterm, i, &fluids);
}

inline void resolveAabbCollision(btScalar stiff, btScalar damp, const btVector3 &vel_eval,
								 btVector3 *acceleration, const btVector3 &normal, btScalar depthOfPenetration)
{
	const btScalar COLLISION_EPSILON = 0.00001f;

	if(depthOfPenetration > COLLISION_EPSILON)
	{
		btScalar adj = stiff * depthOfPenetration - damp * normal.dot(vel_eval);
		
		btVector3 collisionAcceleration = normal;
		collisionAcceleration *= adj;	

		*acceleration += collisionAcceleration;
	}
}
void integrateParticle(const FluidParametersGlobal &FG, const FluidParametersLocal &FL,
					   btScalar speedLimitSquared, btScalar R2, 
					   bool isPlaneGravityEnabled, int particleIndex, FluidParticles *fluids)
{		
	const btScalar speedLimit = FG.sph_limit;
	const btScalar ss = FG.sph_simscale;
	
	const btScalar stiff = FL.m_extstiff;
	const btScalar damp = FL.m_extdamp;
	
	const btVector3 &min = FL.m_volumeMin;
	const btVector3 &max = FL.m_volumeMax;
	
	int i = particleIndex;
	
	//Compute Acceleration		
	btVector3 accel = fluids->m_sph_force[i];
	accel *= FL.m_particleMass;

	//Limit speed
	btScalar speedSquared = accel.length2();
	if(speedSquared > speedLimitSquared) accel *= speedLimit / btSqrt(speedSquared);

	//Apply acceleration to keep particles in the FluidSystem's AABB
	resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3( 1.0, 0.0, 0.0), R2 - ( fluids->m_pos[i].x() - min.x() )*ss );
	resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3(-1.0, 0.0, 0.0), R2 - ( max.x() - fluids->m_pos[i].x() )*ss );
	resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3(0.0,  1.0, 0.0), R2 - ( fluids->m_pos[i].y() - min.y() )*ss );
	resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3(0.0, -1.0, 0.0), R2 - ( max.y() - fluids->m_pos[i].y() )*ss );
	resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3(0.0, 0.0,  1.0), R2 - ( fluids->m_pos[i].z() - min.z() )*ss );
	resolveAabbCollision( stiff, damp, fluids->m_vel_eval[i], &accel, btVector3(0.0, 0.0, -1.0), R2 - ( max.z() - fluids->m_pos[i].z() )*ss );

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

	// Leapfrog Integration ----------------------------
	btVector3 vnext = accel;							
	vnext *= FG.m_timeStep;
	vnext += fluids->m_vel[i];						// v(t+1/2) = v(t-1/2) + a(t) dt
	fluids->m_vel_eval[i] = fluids->m_vel[i];
	fluids->m_vel_eval[i] += vnext;
	fluids->m_vel_eval[i] *= 0.5;					// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
	fluids->m_vel[i] = vnext;
	vnext *= FG.m_timeStep / ss;
	fluids->m_pos[i] += vnext;						// p(t+1) = p(t) + v(t+1/2) dt


	// Euler integration -------------------------------
	// accel += m_Gravity;
	//accel *= FP.m_timeStep;
	//fluids->vel[i] += accel;				// v(t+1) = v(t) + a(t) dt
	//fluids->m_vel_eval[i] += accel;
	//fluids->m_vel_eval[i] *= FP.m_timeStep/d;
	//fluids->m_pos[i] += fluids->m_vel_eval[i];
	//fluids->m_vel_eval[i] = fluids->vel[i]; 
}
void FluidSolver::integrate(const FluidParametersGlobal &FG, const FluidParametersLocal &FL, FluidParticles *fluids)
{
	BT_PROFILE("FluidSolver::integrate()");
	
	btScalar speedLimitSquared = FG.sph_limit*FG.sph_limit;
	btScalar R2 = 2.0f * FG.sph_pradius;
	
	bool isPlaneGravityEnabled = !FG.m_planeGravity.isZero();
	
#ifdef FLUIDS_MULTITHREADED_ENABLED
	PF_IntegrateData IntegrateData(FG, FL, speedLimitSquared, R2, isPlaneGravityEnabled, fluids);
	parallelFor.execute( PF_IntegrateFunction, &IntegrateData, 0, fluids->size() - 1, 256 );
#else
	for(int i = 0; i < fluids->size(); ++i)
		integrateParticle(FG, FL, speedLimitSquared, R2, isPlaneGravityEnabled, i, fluids);
#endif
}

void FluidSolverGridNeighbor::sphComputeForceGrid(const FluidParametersGlobal &FG, FluidSph *fluid)
{
	BT_PROFILE("FluidSolverGridNeighbor::sphComputeForceGrid()");

	const FluidParametersLocal &FL = fluid->getLocalParameters();
	FluidParticles &particles = fluid->internalGetFluidParticles();
	FluidGrid *grid = fluid->internalGetGrid();
	
	const bool isLinkedList = (grid->getGridType() == FT_LinkedList);
	btScalar radius = FG.sph_smoothradius / FG.sph_simscale;
	
	btScalar vterm = FG.m_LapKern * FL.m_viscosity;
		
	for(int i = 0; i < particles.size(); ++i)
	{
		btVector3 force(0, 0, 0);

		FindCellsResult foundCells;
		grid->findCells(particles.m_pos[i], radius, &foundCells);
		for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; cell++) 
		{
			FluidGridIterator &FI = foundCells.m_iterators[cell];
			
			for( int n = FI.m_firstIndex; FluidGridIterator::isIndexValid(n, FI.m_lastIndex); 
					 n = FluidGridIterator::getNextIndex(n, isLinkedList, particles.m_nextFluidIndex) )
			{
				if(i == n)continue; 
				
				btVector3 distance = (particles.m_pos[i] - particles.m_pos[n]) * FG.sph_simscale;		//Simulation-scale distance
				btScalar distanceSquared = distance.length2();
				
				if(FG.m_R2 > distanceSquared) 
				{
					btScalar r = btSqrt(distanceSquared);
					
					btScalar c = FG.sph_smoothradius - r;
					btScalar pterm = -0.5f * c * FG.m_SpikyKern * ( particles.m_pressure[i] + particles.m_pressure[n]) / r;
					btScalar dterm = c * particles.m_density[i] * particles.m_density[n];
					
					btVector3 forceAdded( (pterm * distance.x() + vterm * (particles.m_vel_eval[n].x() - particles.m_vel_eval[i].x())) * dterm,
										  (pterm * distance.y() + vterm * (particles.m_vel_eval[n].y() - particles.m_vel_eval[i].y())) * dterm,
										  (pterm * distance.z() + vterm * (particles.m_vel_eval[n].z() - particles.m_vel_eval[i].z())) * dterm );
					force += forceAdded;
				}
			}
		}
		
		particles.m_sph_force[i] = force;
	}
}


void calculatePressuresInCellSymmetric(const FluidParametersGlobal &FG, const btScalar gridSearchRadius, int gridCellIndex, 
									   FluidGrid *tempGrid, FluidParticles *fluids, btAlignedObjectArray<btScalar> *sums)
{
	const bool isLinkedList = (tempGrid->getGridType() == FT_LinkedList);

	FluidGridIterator currentCell = tempGrid->getGridCell(gridCellIndex);
	for( int i = currentCell.m_firstIndex; FluidGridIterator::isIndexValid(i, currentCell.m_lastIndex); 
			 i = FluidGridIterator::getNextIndex(i, isLinkedList, fluids->m_nextFluidIndex) )
	{
		FindCellsResult foundCells;
		tempGrid->findCells(fluids->m_pos[i], gridSearchRadius, &foundCells);
		for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; ++cell) 
		{
			FluidGridIterator &FI = foundCells.m_iterators[cell];
			
			for( int n = FI.m_firstIndex; FluidGridIterator::isIndexValid(n, FI.m_lastIndex); 
					 n = FluidGridIterator::getNextIndex(n, isLinkedList, fluids->m_nextFluidIndex) )
			{
				if(i == n) continue;
				
				//Simulation-scale distance
				btVector3 distanceAsVector = (fluids->m_pos[i] - fluids->m_pos[n]) * FG.sph_simscale;		
				btScalar distanceSquared = distanceAsVector.length2();
				
				if(FG.m_R2 > distanceSquared) 
				{
					btScalar c = FG.m_R2 - distanceSquared;
					btScalar c_cubed = c * c * c;
					(*sums)[i] += c_cubed;
					(*sums)[n] += c_cubed;
					
					btScalar distance = btSqrt(distanceSquared);
					if( !fluids->m_neighborTable[i].isFilled() ) fluids->m_neighborTable[i].addNeighbor(n, distance);
					else if( !fluids->m_neighborTable[n].isFilled() ) fluids->m_neighborTable[n].addNeighbor(i, distance);
				}
			}
		}
		
		//Remove particle from grid cell
		tempGrid->removeFirstParticle(gridCellIndex, fluids->m_nextFluidIndex);
	}
}

//Remove particles from grid cells as their interactions are calculated
void FluidSolverReducedGridNeighbor::sphComputePressureGridReduce(const FluidParametersGlobal &FG, FluidSph *fluid)
{
	BT_PROFILE("sphComputePressureGridReduce()");
	
	const int numParticles = fluid->numParticles();
	
	const FluidParametersLocal &FL = fluid->getLocalParameters();
	FluidParticles &particles = fluid->internalGetFluidParticles();
	FluidGrid *grid = fluid->internalGetGrid();
	const bool isLinkedList = (grid->getGridType() == FT_LinkedList);
	
	
	FluidGrid *tempGrid = 0;
	static Grid tempStaticGrid;
	static HashGrid tempHashGrid;
	static btAlignedObjectArray<btScalar> sums;
	{
		BT_PROFILE("sphComputePressureGridReduce() - copy grid, reset sums, clear table");
		if(isLinkedList)
		{
			tempStaticGrid = *static_cast<Grid*>(grid);
			tempGrid = static_cast<FluidGrid*>(&tempStaticGrid);
		}
		else 
		{
			tempHashGrid = *static_cast<HashGrid*>(grid);
			tempGrid = static_cast<FluidGrid*>(&tempHashGrid);
		}
		
		sums.resize(numParticles);
		for(int i = 0; i < numParticles; ++i) sums[i] = 0.f;
		for(int i = 0; i < numParticles; ++i) particles.m_neighborTable[i].clear();
	}
	
	{
		BT_PROFILE("sphComputePressureGridReduce() - compute sums");
		btScalar gridSearchRadius = FG.sph_smoothradius / FG.sph_simscale;
		
		for(int group = 0; group < NUM_CELL_PROCESSING_GROUPS; ++group)
		{
			const btAlignedObjectArray<int> &currentGroup = grid->getCellProcessingGroup(group);
		
#ifdef FLUIDS_MULTITHREADED_ENABLED
			PF_ComputePressureData PressureData(FG, gridSearchRadius, currentGroup, tempGrid, &particles, &sums);
			parallelFor.execute( PF_ComputePressureFunction, &PressureData, 0, currentGroup.size() - 1, 1 );
#else
			for(int cell = 0; cell < currentGroup.size(); ++cell)
				calculatePressuresInCellSymmetric(FG, gridSearchRadius, currentGroup[cell], tempGrid, &particles, &sums);
#endif		
		}
	}
	
	{
		BT_PROFILE("sphComputePressureGridReduce() - compute pressure/density");
	
		for(int i = 0; i < numParticles; ++i)
		{
			btScalar tempDensity = sums[i] * FL.m_particleMass * FG.m_Poly6Kern;	
			particles.m_pressure[i] = (tempDensity - FL.m_restDensity) * FL.m_intstiff;
			particles.m_density[i] = 1.0f / tempDensity;
		}
	}
}

void computeForceNeighborTableSymmetric(const FluidParametersGlobal &FG, const btScalar vterm, int particleIndex, FluidParticles *fluids)
{
	int i = particleIndex;

	for(int j = 0; j < fluids->m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = fluids->m_neighborTable[i].getNeighborIndex(j);
		
		btVector3 distance = (fluids->m_pos[i] - fluids->m_pos[n]) * FG.sph_simscale;		//Simulation-scale distance
		
		btScalar c = FG.sph_smoothradius - fluids->m_neighborTable[i].getDistance(j);
		btScalar pterm = -0.5f * c * FG.m_SpikyKern 
					 * ( fluids->m_pressure[i] + fluids->m_pressure[n]) / fluids->m_neighborTable[i].getDistance(j);
		btScalar dterm = c * fluids->m_density[i] * fluids->m_density[n];

		btVector3 forceAdded( (pterm * distance.x() + vterm * (fluids->m_vel_eval[n].x() - fluids->m_vel_eval[i].x())) * dterm,
							  (pterm * distance.y() + vterm * (fluids->m_vel_eval[n].y() - fluids->m_vel_eval[i].y())) * dterm,
							  (pterm * distance.z() + vterm * (fluids->m_vel_eval[n].z() - fluids->m_vel_eval[i].z())) * dterm );
		
		fluids->m_sph_force[i] += forceAdded;
		fluids->m_sph_force[n] += -forceAdded;
	}
}
void calculateForcesInCellSymmetric(const FluidParametersGlobal &FG, const btScalar vterm, bool isLinkedList, 
									int gridCellIndex, const FluidGrid *grid, FluidParticles *fluids)
{
	FluidGridIterator currentCell = grid->getGridCell(gridCellIndex);
	for( int i = currentCell.m_firstIndex; FluidGridIterator::isIndexValid(i, currentCell.m_lastIndex); 
			 i = FluidGridIterator::getNextIndex(i, isLinkedList, fluids->m_nextFluidIndex) )
	{
		computeForceNeighborTableSymmetric(FG, vterm, i, fluids);
	}
}
void FluidSolverReducedGridNeighbor::sphComputeForceReduce(const FluidParametersGlobal &FG, FluidSph *fluid)
{
	BT_PROFILE("FluidSolverReducedGridNeighbor::sphComputeForceReduce()");
	
	const FluidParametersLocal &FL = fluid->getLocalParameters();
	FluidParticles &particles = fluid->internalGetFluidParticles();
	FluidGrid *grid = fluid->internalGetGrid();
	const bool isLinkedList = (grid->getGridType() == FT_LinkedList);
	
	btScalar vterm = FG.m_LapKern * FL.m_viscosity;
	
	for(int i = 0; i < particles.size(); ++i)particles.m_sph_force[i].setValue(0, 0, 0);
	
	for(int group = 0; group < NUM_CELL_PROCESSING_GROUPS; ++group)
	{
		const btAlignedObjectArray<int> &currentGroup = grid->getCellProcessingGroup(group);
		
#ifdef FLUIDS_MULTITHREADED_ENABLED
		PF_ComputeForceData ForceData(FG, vterm, isLinkedList, currentGroup, grid, &particles);
		parallelFor.execute( PF_ComputeForceFunction, &ForceData, 0, currentGroup.size() - 1, 1 );
#else
		for(int cell = 0; cell < currentGroup.size(); ++cell)
			calculateForcesInCellSymmetric(FG, vterm, isLinkedList, currentGroup[cell], grid, &particles);
#endif		
	}
}
