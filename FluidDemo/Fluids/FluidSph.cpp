/*
  FLUIDS v.1 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

  ZLib license
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
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

#include "FluidSph.h"

#include <cmath>
#include <cstdio>

//Link to e.g. 'BulletMultiThreaded.lib' if enabling this
//(must build Bullet with CMake to get BulletMultiThreaded library)
//#define FLUIDS_MULTITHREADED_ENABLED	//Experimental -- may crash intermittently
#ifdef FLUIDS_MULTITHREADED_ENABLED
#include "BulletMultiThreadedSupport.h"

const unsigned int NUM_THREADS = 4;
ParallelFor parallelFor("parallelForThreads", NUM_THREADS);

void processCell(const FluidParametersGlobal &FG, const btScalar gridSearchRadius, int gridCellIndex, 
				 Grid *tempGrid, FluidParticles *fluids, btAlignedObjectArray<btScalar> *sums);					   
void computeForceNeighborTable(const FluidParametersGlobal &FG, const btScalar vterm, int particleIndex, FluidParticles *fluids);				   
void integrateParticle(const FluidParametersGlobal &FG, const FluidParametersLocal &FL,
					   btScalar speedLimitSquared, btScalar R2,  bool isPlaneGravityEnabled, 
					   int particleIndex, FluidParticles *fluids);			
					   
struct PF_ComputePressureData
{
	const FluidParametersGlobal &m_globalParameters; 
	const btScalar m_gridSearchRadius;
	const btAlignedObjectArray<int> &m_gridCellGroup;
	Grid *m_tempGrid;
	FluidParticles *m_particles; 
	btAlignedObjectArray<btScalar> *m_sums;
	
	PF_ComputePressureData(const FluidParametersGlobal &FG, const btScalar gridSearchRadius,
						   const btAlignedObjectArray<int> &gridCellGroup, Grid *tempGrid,
						   FluidParticles *particles, btAlignedObjectArray<btScalar> *sums) 
	: m_globalParameters(FG), m_gridSearchRadius(gridSearchRadius),
	  m_gridCellGroup(gridCellGroup),m_tempGrid(tempGrid), 
	  m_particles(particles), m_sums(sums) {}
};
void PF_ComputePressureFunction(void *parameters, int index)
{
	PF_ComputePressureData *data = static_cast<PF_ComputePressureData*>(parameters);
	
	processCell(data->m_globalParameters, data->m_gridSearchRadius, 
				data->m_gridCellGroup[index], data->m_tempGrid, data->m_particles, data->m_sums);
}

struct PF_ComputeForceData
{
	const FluidParametersGlobal &m_globalParameters; 
	const btScalar m_vterm;
	FluidParticles *m_particles;
	
	PF_ComputeForceData(const FluidParametersGlobal &FG, const btScalar vterm, FluidParticles *particles) 
	: m_globalParameters(FG), m_vterm(vterm), m_particles(particles) {}
};
void PF_ComputeForceFunction(void *parameters, int index)
{
	PF_ComputeForceData *data = static_cast<PF_ComputeForceData*>(parameters);
	
	computeForceNeighborTable(data->m_globalParameters, data->m_vterm, index, data->m_particles);
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

#endif	//FLUIDS_MULTITHREADED_ENABLED

void FluidSph::initialize(const FluidParametersGlobal &FG, int maxNumParticles, const btVector3 &volumeMin, const btVector3 &volumeMax)
{
	reset(maxNumParticles);	
	
	m_localParameters.m_volumeMin = volumeMin;
	m_localParameters.m_volumeMax = volumeMax;

	btScalar simCellSize = FG.sph_smoothradius * 2.0;					//Grid cell size (2r)	
	m_grid.setup(volumeMin, volumeMax,  FG.sph_simscale, simCellSize, 1.0);
	
#ifdef USE_HASHGRID
	m_hashgrid.setup(FG.sph_simscale, simCellSize);
#endif
}

void FluidSph::stepSimulation(const FluidParametersGlobal &FG)
{
	BT_PROFILE("FluidSph::stepSimulation()");

	removeMarkedFluids();
	
	//CPU Branch
	{
#ifdef USE_HASHGRID
		grid_insertParticles_HASHGRID();
		
		sph_computePressureGrid_HASHGRID(FG);
		
		grid_insertParticles();		//	External interface is not implemented for hashgrid(e.g. marching cubes)
#else
		grid_insertParticles();
		
		//sph_computePressureGrid(FG);
		sph_computePressureGridReduce(FG);
#endif	
		sph_computeForceGridNC(FG);
		
		integrate(FG);
	}
}

void FluidSph::reset(int maxNumParticles)
{
	clear();
	
	m_particles.resize(0);
	m_particles.setMaxParticles(maxNumParticles);
	
	setDefaultParameters();
}

void FluidSph::clear()
{
	m_particles.resize(0);
	
	m_removedFluidIndicies.resize(0);
	
	m_grid.clear();
}

//------------------------------------------------------ SPH Setup 
//
//  Range = +/- 10.0 * 0.006 (r) =	   0.12			m (= 120 mm = 4.7 inch)
//  Container Volume (Vc) =			   0.001728		m^3
//  Rest Density (D) =				   1000.0		kg / m^3
//  Particle Mass (Pm) =			   0.00020543	kg		(mass = vol * density)
//  Number of Particles (N) =		   4000.0
//  Water Mass (M) =				   0.821		kg (= 821 grams)
//  Water Volume (V) =				   0.000821     m^3 (= 3.4 cups, .21 gals)
//  Smoothing Radius (R) =             0.02			m (= 20 mm = ~3/4 inch)
//  Particle Radius (Pr) =			   0.00366		m (= 4 mm  = ~1/8 inch)
//  Particle Volume (Pv) =			   2.054e-7		m^3	(= .268 milliliters)
//  Rest Distance (Pd) =			   0.0059		m
//
//  Given: D, Pm, N
//    Pv = Pm / D			0.00020543 kg / 1000 kg/m^3 = 2.054e-7 m^3	
//    Pv = 4/3*pi*Pr^3    cuberoot( 2.054e-7 m^3 * 3/(4pi) ) = 0.00366 m
//     M = Pm * N			0.00020543 kg * 4000.0 = 0.821 kg		
//     V =  M / D              0.821 kg / 1000 kg/m^3 = 0.000821 m^3
//     V = Pv * N			 2.054e-7 m^3 * 4000 = 0.000821 m^3
//    Pd = cuberoot(Pm/D)    cuberoot(0.00020543/1000) = 0.0059 m 
//
// Ideal grid cell size (gs) = 2 * smoothing radius = 0.02*2 = 0.04
// Ideal domain size = k*gs/d = k*0.02*2/0.005 = k*8 = {8, 16, 24, 32, 40, 48, ..}
//    (k = number of cells, gs = cell size, d = simulation scale)

btScalar FluidSph::getValue(btScalar x, btScalar y, btScalar z) const
{
	btScalar sum = 0.0;
	
	const btScalar searchRadius = m_grid.getParameters().m_gridCellSize / 2.0f;
	const btScalar R2 = 1.8*1.8;
	//const btScalar R2 = 0.8*0.8;		//	marching cubes rendering test

	GridCellIndicies GC;
	m_grid.findCells( btVector3(x,y,z), searchRadius, &GC );
	
	for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; ++cell) 
	{
		for(int n = m_grid.getLastParticleIndex( GC.m_indicies[cell] ); n != INVALID_PARTICLE_INDEX; n = m_particles.m_nextFluidIndex[n])
		{
			const btVector3 &position = m_particles.m_pos[n];
			btScalar dx = x - position.x();
			btScalar dy = y - position.y();
			btScalar dz = z - position.z();
			btScalar dsq = dx*dx+dy*dy+dz*dz;
				
			if(dsq < R2) sum += R2 / dsq;
		}
	}
	
	return sum;	
}	
btVector3 FluidSph::getGradient(btScalar x, btScalar y, btScalar z) const
{
	btVector3 norm(0,0,0);
	
	const btScalar searchRadius = m_grid.getParameters().m_gridCellSize / 2.0f;
	const btScalar R2 = searchRadius*searchRadius;
	
	GridCellIndicies GC;
	m_grid.findCells( btVector3(x,y,z), searchRadius, &GC );
	for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; cell++)
	{
		for(int n = m_grid.getLastParticleIndex( GC.m_indicies[cell] ); n != INVALID_PARTICLE_INDEX; n = m_particles.m_nextFluidIndex[n])
		{
			const btVector3 &position = m_particles.m_pos[n];
			btScalar dx = x - position.x();
			btScalar dy = y - position.y();
			btScalar dz = z - position.z();
			btScalar dsq = dx*dx+dy*dy+dz*dz;
			
			if(0 < dsq && dsq < R2) 
			{
				dsq = 2.0*R2 / (dsq*dsq);
				
				btVector3 particleNorm(dx * dsq, dy * dsq, dz * dsq);
				norm += particleNorm;
			}
		}
	}
	
	norm.normalize();
	return norm;
}

//for btAlignedObjectArray<int>::heapSort()/quickSort()
struct DescendingSortPredicate { inline bool operator() (const int &a, const int &b) const { return !(a < b); } };
void FluidSph::removeMarkedFluids()
{
	//Since removing elements from the array invalidates(higher) indicies,
	//elements should be removed in descending order.
	m_removedFluidIndicies.heapSort( DescendingSortPredicate() );
	//m_removedFluidIndicies.quickSort( DescendingSortPredicate() );	//	crashes (issue with btAlignedObjectArray<int>::quickSort())
	
	//Remove duplicates; rearrange array such that all unique indicies are
	//in the range [0, uniqueSize) and in descending order.
	int uniqueSize = 0;
	if( m_removedFluidIndicies.size() ) 
	{
		uniqueSize = 1;
		for(int i = 1; i < m_removedFluidIndicies.size(); ++i)
		{
			if( m_removedFluidIndicies[i] != m_removedFluidIndicies[i-1] )
			{
				m_removedFluidIndicies[uniqueSize] = m_removedFluidIndicies[i];
				++uniqueSize;
			}
		}
	}
	
	//
	for(int i = 0; i < uniqueSize; ++i) m_particles.removeParticle( m_removedFluidIndicies[i] );
	
	m_removedFluidIndicies.resize(0);
}

void FluidSph::grid_insertParticles()
{	
	BT_PROFILE("FluidSph::grid_insertParticles()");

	//Reset particles
	for(int i = 0; i < numParticles(); ++i) m_particles.m_nextFluidIndex[i] = INVALID_PARTICLE_INDEX;	

	//Reset grid
	m_grid.clear();
	
	//Insert particles into grid
	for(int i = 0; i < numParticles(); ++i) 
		m_grid.insertParticle(m_particles.m_pos[i], i, &m_particles.m_nextFluidIndex[i]);
}


//Compute Pressures - Using spatial grid, and also create neighbor table
void FluidSph::sph_computePressureGrid(const FluidParametersGlobal &FG)
{
	BT_PROFILE("FluidSph::sph_computePressureGrid()");
	
	btScalar radius = FG.sph_smoothradius / FG.sph_simscale;

	for(int i = 0; i < numParticles(); ++i)
	{
		btScalar sum = 0.0;	
		m_particles.m_neighborTable[i].clear();

		GridCellIndicies GC;
		m_grid.findCells(m_particles.m_pos[i], radius, &GC);
		for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; cell++) 
		{
			for(int n = m_grid.getLastParticleIndex( GC.m_indicies[cell] );	n != INVALID_PARTICLE_INDEX; n = m_particles.m_nextFluidIndex[n])
			{
				if(i == n) continue;
				
				btVector3 distance = (m_particles.m_pos[i] - m_particles.m_pos[n]) * FG.sph_simscale;		//Simulation-scale distance
				btScalar distanceSquared = distance.length2();
				
				if(FG.m_R2 > distanceSquared) 
				{
					btScalar c = FG.m_R2 - distanceSquared;
					sum += c * c * c;
					
					if( !m_particles.m_neighborTable[i].isFilled() ) m_particles.m_neighborTable[i].addNeighbor( n, btSqrt(distanceSquared) );
				}
			}
		}
		
		btScalar tempDensity = sum * m_localParameters.m_particleMass * FG.m_Poly6Kern;	
		m_particles.m_pressure[i] = (tempDensity - m_localParameters.m_restDensity) * m_localParameters.m_intstiff;
		m_particles.m_density[i] = 1.0f / tempDensity;
	}
}
#ifdef USE_HASHGRID
void FluidSph::sph_computePressureGrid_HASHGRID(const FluidParametersGlobal &FG)
{
	BT_PROFILE("FluidSph::sph_computePressureGrid - hashgrid()");
	
	btScalar radius = FG.sph_smoothradius / FG.sph_simscale;

	for(int i = 0; i < numParticles(); ++i)
	{
		btScalar sum = 0.0;	
		m_particles.m_neighborTable[i].clear();

		HashGridQueryResult HG;
		m_hashgrid.findCells(m_particles.m_pos[i], radius, &HG);
		for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; ++cell) 
		{
			if(!HG.m_cells[cell]) continue;
			
			for(int n = HG.m_cells[cell]->m_firstIndex; n <= HG.m_cells[cell]->m_lastIndex; ++n)
			{
				if(i == n)continue; 
				
				//Simulation-scale distance
				btVector3 distance = (m_particles.m_pos[i] - m_particles.m_pos[n]) * FG.sph_simscale;
				btScalar distanceSquared = distance.length2();
				
				if(FG.m_R2 > distanceSquared) 
				{
					btScalar c = FG.m_R2 - distanceSquared;
					sum += c * c * c;
					
					if( !m_particles.m_neighborTable[i].isFilled() ) m_particles.m_neighborTable[i].addNeighbor( n, btSqrt(distanceSquared) );
				}
			}
		}
		
		btScalar tempDensity = sum * m_localParameters.m_particleMass * FG.m_Poly6Kern;	
		m_particles.m_pressure[i] = (tempDensity - m_localParameters.m_restDensity) * m_localParameters.m_intstiff;
		m_particles.m_density[i] = 1.0f / tempDensity;
	}
}
#endif

void getCellProcessingGroups(const Grid &G, btAlignedObjectArray<int> *gridCellIndicies)
{
	//Although a single particle only accesses 2^3 grid cells,
	//the particles within a single grid cell may access up to 
	//3^3 grid cells overall.
	//
	//In order to simultaneously calculate pressure on a per grid 
	//cell basis, with particles being removed from a grid cell
	//as it is processed, the cells must be split into 27
	//sequentially processed groups.

	BT_PROFILE("getCellProcessingGroups()");
	
	const GridParameters &GP = G.getParameters();
	
	for(int i = 0; i < 27; ++i) gridCellIndicies[i].resize(0);
	
	int resX = GP.m_resolutionX;
	int resY = GP.m_resolutionY;
	int resZ = GP.m_resolutionZ;
	for(int cell = 0; cell < GP.m_numCells; ++cell)
	{
		if( G.getLastParticleIndex(cell) == INVALID_PARTICLE_INDEX ) continue;
	
		char cellX = cell % resX;
		char cellZ = cell / (resX*resY);
		char cellY = (cell - cellZ*resX*resY) / resX;

		//Convert range from [-128, 127] to [1, 256]	//	for HashGrid
		//cellX += 129;
		//cellY += 129;
		//cellZ += 129;

		//Convert range from [0, n] to [1, n+1]
		++cellX;
		++cellY;
		++cellZ;
		
		char group = 0;
		if(cellX % 3 == 0) group += 0;
		else if(cellX % 2 == 0) group += 1;
		else group += 2;
		if(cellY % 3 == 0) group += 0;
		else if(cellY % 2 == 0) group += 3;
		else group += 6;
		if(cellZ % 3 == 0) group += 0;
		else if(cellZ % 2 == 0) group += 9;
		else group += 18;
		
		gridCellIndicies[group].push_back(cell);
	}
}
void processCell(const FluidParametersGlobal &FG, const btScalar gridSearchRadius, int gridCellIndex, 
				 Grid *tempGrid, FluidParticles *fluids, btAlignedObjectArray<btScalar> *sums)
{
	for(int i = tempGrid->getLastParticleIndex(gridCellIndex); i != INVALID_PARTICLE_INDEX; i = fluids->m_nextFluidIndex[i])
	{
		GridCellIndicies GC;
		tempGrid->findCells(fluids->m_pos[i], gridSearchRadius, &GC);
		for(int foundCell = 0; foundCell < RESULTS_PER_GRID_SEARCH; ++foundCell) 
		{
			for( int n = tempGrid->getLastParticleIndex(GC.m_indicies[foundCell]); 
				 n != INVALID_PARTICLE_INDEX; n = fluids->m_nextFluidIndex[n] )
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
					if( !fluids->m_neighborTable[n].isFilled() ) fluids->m_neighborTable[n].addNeighbor(i, distance);
				}
			}
		}
		
		//Remove particle from grid cell
		int lastParticleIndex = tempGrid->getLastParticleIndex(gridCellIndex);
		if(lastParticleIndex != INVALID_PARTICLE_INDEX)
		{
			int nextIndex = fluids->m_nextFluidIndex[lastParticleIndex];
			tempGrid->setLastParticleIndex(gridCellIndex, nextIndex);
		}
	}
}
//Remove particles from grid cells as their interactions are calculated
void FluidSph::sph_computePressureGridReduce(const FluidParametersGlobal &FG)
{
	BT_PROFILE("FluidSph::sph_computePressureGridReduce()");
	
	static Grid tempGrid;
	static btAlignedObjectArray<btScalar> sums;
	{
		BT_PROFILE("sph_computePressureGridReduce() - copy grid, reset sums, clear table");
		tempGrid = m_grid;
		
		sums.resize( numParticles() );
		for(int i = 0; i < numParticles(); ++i) sums[i] = 0.f;
		for(int i = 0; i < numParticles(); ++i) m_particles.m_neighborTable[i].clear();
	}
	
	static btAlignedObjectArray<int> gridCellIndicies[27];
	getCellProcessingGroups(tempGrid, gridCellIndicies);
	
	{
		BT_PROFILE("sph_computePressureGridReduce() - compute sums");
		btScalar gridSearchRadius = FG.sph_smoothradius / FG.sph_simscale;
		
		for(int group = 0; group < 27; ++group)
		{
#ifdef FLUIDS_MULTITHREADED_ENABLED
			PF_ComputePressureData PressureData(FG, gridSearchRadius, gridCellIndicies[group], &tempGrid, &m_particles, &sums);
			parallelFor.execute( PF_ComputePressureFunction, &PressureData, 0, gridCellIndicies[group].size() - 1, 1 );
#else
			for(int cell = 0; cell < gridCellIndicies[group].size(); ++cell)
				processCell(FG, gridSearchRadius, gridCellIndicies[group][cell], &tempGrid, &m_particles, &sums);
#endif		
		}
	}
	
	{
		BT_PROFILE("sph_computePressureGridReduce() - compute pressure/density");
	
		for(int i = 0; i < numParticles(); ++i)
		{
			btScalar tempDensity = sums[i] * m_localParameters.m_particleMass * FG.m_Poly6Kern;	
			m_particles.m_pressure[i] = (tempDensity - m_localParameters.m_restDensity) * m_localParameters.m_intstiff;
			m_particles.m_density[i] = 1.0f / tempDensity;
		}
	}
}

//Compute Forces - Using spatial grid. Faster.
void FluidSph::sph_computeForceGrid(const FluidParametersGlobal &FG)
{
	BT_PROFILE("FluidSph::sph_computeForceGrid()");

	btScalar radius = FG.sph_smoothradius / FG.sph_simscale;
	btScalar vterm = FG.m_LapKern * m_localParameters.m_viscosity;
		
	for(int i = 0; i < numParticles(); ++i)
	{
		btVector3 force(0, 0, 0);

		GridCellIndicies GC;
		m_grid.findCells(m_particles.m_pos[i], radius, &GC);
		for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; cell++) 
		{
			for(int n = m_grid.getLastParticleIndex( GC.m_indicies[cell] ); n != INVALID_PARTICLE_INDEX; n = m_particles.m_nextFluidIndex[n])
			{
				if(i == n)continue; 
				
				btVector3 distance = (m_particles.m_pos[i] - m_particles.m_pos[n]) * FG.sph_simscale;		//Simulation-scale distance
				btScalar distanceSquared = distance.length2();
				
				if(FG.m_R2 > distanceSquared) 
				{
					btScalar r = btSqrt(distanceSquared);
					
					btScalar c = FG.sph_smoothradius - r;
					btScalar pterm = -0.5f * c * FG.m_SpikyKern * ( m_particles.m_pressure[i] + m_particles.m_pressure[n]) / r;
					btScalar dterm = c * m_particles.m_density[i] * m_particles.m_density[n];
					
					btVector3 forceAdded( (pterm * distance.x() + vterm * (m_particles.m_vel_eval[n].x() - m_particles.m_vel_eval[i].x())) * dterm,
										  (pterm * distance.y() + vterm * (m_particles.m_vel_eval[n].y() - m_particles.m_vel_eval[i].y())) * dterm,
										  (pterm * distance.z() + vterm * (m_particles.m_vel_eval[n].z() - m_particles.m_vel_eval[i].z())) * dterm );
					force += forceAdded;
				}
			}
		}
		
		m_particles.m_sph_force[i] = force;
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
//Compute Forces - Using spatial grid with saved neighbor table. Fastest.
void FluidSph::sph_computeForceGridNC(const FluidParametersGlobal &FG)
{
	BT_PROFILE("FluidSph::sph_computeForceGridNC()");

	btScalar vterm = FG.m_LapKern * m_localParameters.m_viscosity;
	
#ifdef FLUIDS_MULTITHREADED_ENABLED
	PF_ComputeForceData ForceData(FG, vterm, &m_particles);
	parallelFor.execute( PF_ComputeForceFunction, &ForceData, 0, numParticles() - 1, 256 );
#else
	for(int i = 0; i < numParticles(); ++i)
		computeForceNeighborTable(FG, vterm, i, &m_particles);
#endif	
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

	//CCD_TEST
	fluids->m_prev_pos[i] = fluids->m_pos[i];
	
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
void FluidSph::integrate(const FluidParametersGlobal &FG)
{
	BT_PROFILE("FluidSph::integrate()");
	
	btScalar speedLimitSquared = FG.sph_limit*FG.sph_limit;
	btScalar R2 = 2.0f * FG.sph_pradius;
	
	bool isPlaneGravityEnabled = !FG.m_planeGravity.isZero();
	
#ifdef FLUIDS_MULTITHREADED_ENABLED
	PF_IntegrateData IntegrateData(FG, m_localParameters, speedLimitSquared, R2, isPlaneGravityEnabled, &m_particles);
	parallelFor.execute( PF_IntegrateFunction, &IntegrateData, 0, numParticles() - 1, 256 );
#else
	for(int i = 0; i < numParticles(); ++i)
		integrateParticle(FG, m_localParameters, speedLimitSquared, R2, isPlaneGravityEnabled, i, &m_particles);
#endif	
}

////////////////////////////////////////////////////////////////////////////////
/// struct FluidEmitter
////////////////////////////////////////////////////////////////////////////////
void FluidEmitter::emit(FluidSph *fluid, int numParticles, btScalar spacing)
{
	int x = static_cast<int>( btSqrt(static_cast<btScalar>(numParticles)) );
	
	for(int i = 0; i < numParticles; i++) 
	{
		btScalar ang_rand = ( static_cast<btScalar>(rand()*2.0/RAND_MAX) - 1.0 ) * m_yawSpread;
		btScalar tilt_rand = ( static_cast<btScalar>(rand()*2.0/RAND_MAX) - 1.0 ) * m_pitchSpread;
		
		//Modified - set y to vertical axis
		btVector3 dir( 	btCos((m_yaw + ang_rand) * SIMD_RADS_PER_DEG) * btSin((m_pitch + tilt_rand) * SIMD_RADS_PER_DEG) * m_velocity,
						btCos((m_pitch + tilt_rand) * SIMD_RADS_PER_DEG) * m_velocity,
						btSin((m_yaw + ang_rand) * SIMD_RADS_PER_DEG) * btSin((m_pitch + tilt_rand) * SIMD_RADS_PER_DEG) * m_velocity );
		
		btVector3 position( spacing*(i/x), spacing*(i%x), 0 );
		position += m_position;
		
		int index = fluid->addParticle(position);
		fluid->setVelocity(index, dir);
	}
}

void FluidEmitter::addVolume(FluidSph *fluid, const btVector3 &min, const btVector3 &max, btScalar spacing)
{
	for(btScalar z = max.z(); z >= min.z(); z -= spacing) 
		for(btScalar y = min.y(); y <= max.y(); y += spacing) 
			for(btScalar x = min.x(); x <= max.x(); x += spacing) 
			{
				fluid->addParticle( btVector3(x,y,z) );
			}
}


