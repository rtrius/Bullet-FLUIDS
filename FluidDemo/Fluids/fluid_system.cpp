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

#include "fluid_system.h"

#include <cmath>
#include <cstdio>

#define FLUIDS_OPENCL_ENABLED
#ifdef FLUIDS_OPENCL_ENABLED
	#include "OpenCL_support/fluids_opencl_support.h"
#endif

void FluidSystem::initialize(int maxNumParticles, const Vector3DF &volumeMin, const Vector3DF &volumeMax)
{
	m_particles.reserve(maxNumParticles);
	m_neighborTable.reserve(maxNumParticles);
	
	sph_setup();
	
	reset(maxNumParticles);	
	
	m_parameters.m_volumeMin = volumeMin;
	m_parameters.m_volumeMax = volumeMax;

	float simCellSize = m_parameters.sph_smoothradius * 2.0;					//Grid cell size (2r)	
	m_grid.setup(volumeMin, volumeMax,  m_parameters.sph_simscale, simCellSize, 1.0);
}

void FluidSystem::reset(int maxNumParticles)
{
#ifndef USE_BTVECTOR3_IN_FLUIDS
	m_particles.clear();
	m_neighborTable.clear();
#else
	m_particles.resize(0);
	m_neighborTable.clear();
#endif
	
	m_maxParticles = maxNumParticles;
	m_particles.reserve(maxNumParticles);
	m_neighborTable.reserve(maxNumParticles);
	
	m_parameters.m_timeStep = 0.003; //  0.001;			// .001 = for point grav

	//Reset parameters
	m_parameters.sph_visc = 0.2;
	m_parameters.sph_intstiff = 0.50;
	m_parameters.sph_extstiff = 20000.0;
	m_parameters.sph_smoothradius = 0.01;
	
	m_parameters.m_pointGravity = 0.0;
	m_parameters.m_pointGravityPosition.setValue(0, 0, 50);
	m_parameters.m_planeGravity.setValue(0, -9.8, 0);
	
	//
	m_grid.clear();
}

/*
int FluidSystem::addPoint()
{
	m_particles.push_back( Fluid() );
	m_neighborTable.push_back( Neighbors() );
	
	Fluid *f = &m_particles[ m_particles.size() - 1 ];
	
	//f->pos.setValue(0,0,0);	//Should be set by caller of addPoint();
	f->vel.setValue(0,0,0);
	f->vel_eval.setValue(0,0,0);
	f->sph_force.setValue(0,0,0);
	f->externalAcceleration.setValue(0,0,0);
	f->nextFluidIndex = INVALID_PARTICLE_INDEX;
	f->pressure = 0;
	f->density = 0;
	
	int particleIndex = m_particles.size() - 1;
	return particleIndex;
}
*/

int FluidSystem::addPointReuse(const Vector3DF &position)
{
	int particleIndex;
	Fluid *f;
	if( m_particles.size() < m_maxParticles )
	{
		m_particles.push_back( Fluid() );
		m_neighborTable.push_back( Neighbors() );
		
		particleIndex = m_particles.size() - 1;
		f = &m_particles[ m_particles.size() - 1 ];
	}
	else
	{
		particleIndex = (m_particles.size()-1) * rand() / RAND_MAX;		//Random index
		f = &m_particles[particleIndex];
	}
	
	f->pos = position;
	f->prev_pos = position;		//CCD_TEST
	f->vel.setValue(0,0,0);
	f->vel_eval.setValue(0,0,0);
	f->sph_force.setValue(0,0,0);
	f->externalAcceleration.setValue(0,0,0);
	f->nextFluidIndex = INVALID_PARTICLE_INDEX;
	f->pressure = 0;
	f->density = 0;
	
	return particleIndex;
}

void FluidSystem::removeMarkedFluids()
{
	//Descending order sort since FluidSystem::removeFluid() swaps elements -- indicies will change
	std::sort( m_removedFluidIndicies.begin(), m_removedFluidIndicies.end(), descendingSort );
	
	//Remove duplicates
	std::vector<int>::iterator end = std::unique( m_removedFluidIndicies.begin(), m_removedFluidIndicies.end() );
	
	//printf("removeMarkedFluids()\n");
	for(std::vector<int>::iterator i = m_removedFluidIndicies.begin(); i != end; ++i) 
	{
		//printf("m_removedFluidIndicies[]: %d \n", *i);
		removeFluid(*i);
	}
	//printf("\n");
	
	m_removedFluidIndicies.clear();
}
void FluidSystem::removeFluid(int index)
{
	int lastIndex = m_particles.size() - 1;
	
	if(index < lastIndex) 
	{
		m_particles[index] = m_particles[lastIndex];
		m_neighborTable[index] = m_neighborTable[lastIndex];
	}
	m_particles.pop_back();
	m_neighborTable.pop_back();
}

void FluidSystem::stepSimulation()
{
	BT_PROFILE("FluidSystem::stepSimulation()");

	removeMarkedFluids();
	
	if(m_useOpenCL) 					//GPU Branch
	{
		if( !m_particles.size() ) return;
	
		GridParameters GP = m_grid.getParameters();
		int *gridCells = m_grid.getCellsPointer();
		int *gridCellsNumFluids = m_grid.getCellsNumFluidsPointer();
		
		FluidParameters_float FP(m_parameters);
		int numFluidParticles = m_particles.size();
		Fluid *fluids = &m_particles[0];
		Neighbors *neighbors = &m_neighborTable.front();

#ifdef FLUIDS_OPENCL_ENABLED
		static FluidSystem_OpenCL *pCL_System = 0;
		if(!pCL_System)
		{
			pCL_System = new FluidSystem_OpenCL;
			pCL_System->initialize();
			//pCL_System->deactivate();
		}
		
		pCL_System->stepSimulation(&FP, &GP,
								   numFluidParticles, fluids, neighbors,
								   GP.m_numCells, gridCells, gridCellsNumFluids);
#endif
	} 
	else 								//CPU Branch
	{
		const bool USING_GRID = true;
		if(USING_GRID)
		{
			grid_insertParticles();
			sph_computePressureGrid();	
			sph_computeForceGridNC();
		}
		else
		{
			//Slow method - O(n^2)
			sph_computePressureSlow();
			sph_computeForceSlow();
		}
		
		advance();
	}
}

inline void resolveAabbCollision(float stiff, float damp, const Vector3DF &vel_eval,
							 	 Vector3DF *acceleration, const Vector3DF &normal, float depthOfPenetration)
{
	const float COLLISION_EPSILON = 0.00001f;

	if(depthOfPenetration > COLLISION_EPSILON)
	{
		double adj = stiff * depthOfPenetration - damp * normal.dot(vel_eval);
		
		Vector3DF collisionAcceleration = normal;
		collisionAcceleration *= adj;	

		*acceleration += collisionAcceleration;
	}
}
void FluidSystem::advance()
{
	BT_PROFILE("FluidSystem::advance()");
	
	float speedLimit = m_parameters.sph_limit;
	float speedLimitSquared = speedLimit*speedLimit;
	
	float stiff = m_parameters.sph_extstiff;
	float damp = m_parameters.sph_extdamp;
	float radius = m_parameters.sph_pradius;
	float R2 = 2.0f * radius;
	float ss = m_parameters.sph_simscale;
	
	const Vector3DF &min = m_parameters.m_volumeMin;
	const Vector3DF &max = m_parameters.m_volumeMax;
	
	bool isPlaneGravityEnabled = ( m_parameters.m_planeGravity.x() != 0.0 
								|| m_parameters.m_planeGravity.y() != 0.0 
								|| m_parameters.m_planeGravity.z() != 0.0 );
	
	for(int i = 0; i < m_particles.size(); ++i)
	{
		Fluid *p = &m_particles[i];

		//CCD_TEST
		p->prev_pos = p->pos;
		
		//Compute Acceleration		
		Vector3DF accel = p->sph_force;
		accel *= m_parameters.sph_pmass;

		//Limit speed
		float speedSquared = accel.length2();
		if(speedSquared > speedLimitSquared) accel *= speedLimit / sqrt(speedSquared);
	
		//Apply acceleration to keep particles in the FluidSystem's AABB
		resolveAabbCollision( stiff, damp, p->vel_eval, &accel, Vector3DF( 1.0, 0.0, 0.0), R2 - ( p->pos.x() - min.x() )*ss );
		resolveAabbCollision( stiff, damp, p->vel_eval, &accel, Vector3DF(-1.0, 0.0, 0.0), R2 - ( max.x() - p->pos.x() )*ss );
		resolveAabbCollision( stiff, damp, p->vel_eval, &accel, Vector3DF(0.0,  1.0, 0.0), R2 - ( p->pos.y() - min.y() )*ss );
		resolveAabbCollision( stiff, damp, p->vel_eval, &accel, Vector3DF(0.0, -1.0, 0.0), R2 - ( max.y() - p->pos.y() )*ss );
		resolveAabbCollision( stiff, damp, p->vel_eval, &accel, Vector3DF(0.0, 0.0,  1.0), R2 - ( p->pos.z() - min.z() )*ss );
		resolveAabbCollision( stiff, damp, p->vel_eval, &accel, Vector3DF(0.0, 0.0, -1.0), R2 - ( max.z() - p->pos.z() )*ss );
	
		//Plane gravity
		if(isPlaneGravityEnabled) accel += m_parameters.m_planeGravity;

		//Point gravity
		if(m_parameters.m_pointGravity > 0.0) 
		{
			Vector3DF norm( p->pos.x() - m_parameters.m_pointGravityPosition.x(),
							p->pos.y() - m_parameters.m_pointGravityPosition.y(),
							p->pos.z() - m_parameters.m_pointGravityPosition.z() );
			norm.normalize();
			norm *= m_parameters.m_pointGravity;
			accel -= norm;
		}
		
		//Apply external forces
		accel += p->externalAcceleration;
		p->externalAcceleration.setValue(0, 0, 0);

		// Leapfrog Integration ----------------------------
		Vector3DF vnext = accel;							
		vnext *= m_parameters.m_timeStep;
		vnext += p->vel;						// v(t+1/2) = v(t-1/2) + a(t) dt
		p->vel_eval = p->vel;
		p->vel_eval += vnext;
		p->vel_eval *= 0.5;						// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
		p->vel = vnext;
		vnext *= m_parameters.m_timeStep/ss;
		p->pos += vnext;						// p(t+1) = p(t) + v(t+1/2) dt


		// Euler integration -------------------------------
		// accel += m_Gravity;
		//accel *= m_parameters.m_timeStep;
		//p->vel += accel;				// v(t+1) = v(t) + a(t) dt
		//p->vel_eval += accel;
		//p->vel_eval *= m_parameters.m_timeStep/d;
		//p->pos += p->vel_eval;
		//p->vel_eval = p->vel; 
	}
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

void FluidSystem::sph_setup()
{
	m_parameters.sph_simscale =		0.004;			// unit size
	m_parameters.sph_visc =			0.2;			// pascal-second (Pa.s) = 1 kg m^-1 s^-1  (see wikipedia page on viscosity)
	m_parameters.sph_restdensity =	600.0;			// kg / m^3
	m_parameters.sph_pmass =		0.00020543;		// kg
	m_parameters.sph_pradius =		0.004;			// m
	m_parameters.sph_smoothradius =	0.01;			// m 
	//m_parameters.sph_intstiff =		1.00;		//	overwritten by FluidSystem::reset()
	m_parameters.sph_extstiff =		10000.0;
	m_parameters.sph_extdamp =		256.0;
	m_parameters.sph_limit =		200.0;			// m / s

	sph_computeKernels();
}

void FluidSystem::sph_computeKernels()
{
	m_parameters.sph_pdist = pow(m_parameters.sph_pmass/m_parameters.sph_restdensity, 1.0/3.0);
	m_parameters.m_R2 = m_parameters.sph_smoothradius * m_parameters.sph_smoothradius;
	
	//Wpoly6 kernel (denominator part) - 2003 Muller, p.4
	m_parameters.m_Poly6Kern = 315.0f / (64.0f * 3.141592 * pow( m_parameters.sph_smoothradius, 9) );	
	
	//Laplacian of viscocity (denominator): PI h^6
	m_parameters.m_SpikyKern = -45.0f / (3.141592 * pow( m_parameters.sph_smoothradius, 6) );	
	
	m_parameters.m_LapKern = 45.0f / (3.141592 * pow( m_parameters.sph_smoothradius, 6) );
}

//Compute Pressures - Very slow yet simple. O(n^2)
void FluidSystem::sph_computePressureSlow()
{
	for(int i = 0; i < m_particles.size(); ++i)
	{
		Fluid *p = &m_particles[i];

		double sum = 0.0;
		for(int n = 0; n < m_particles.size(); ++n)
		{
			Fluid *q = &m_particles[n];
			if(p == q) continue;
			
			double dx = (p->pos.x() - q->pos.x()) * m_parameters.sph_simscale;		//Simulation-scale distance
			double dy = (p->pos.y() - q->pos.y()) * m_parameters.sph_simscale;
			double dz = (p->pos.z() - q->pos.z()) * m_parameters.sph_simscale;
			double dsq = dx*dx + dy*dy + dz*dz;
			if(m_parameters.m_R2 > dsq) 
			{
				double c = m_parameters.m_R2 - dsq;
				sum += c * c * c;
			}
		}	
		p->density = sum * m_parameters.sph_pmass * m_parameters.m_Poly6Kern;	
		p->pressure = (p->density - m_parameters.sph_restdensity) * m_parameters.sph_intstiff;
		p->density = 1.0f / p->density;
	}
}

//Compute Pressures - Using spatial grid, and also create neighbor table
void FluidSystem::sph_computePressureGrid()
{
	BT_PROFILE("FluidSystem::sph_computePressureGrid()");
	
	float radius = m_parameters.sph_smoothradius / m_parameters.sph_simscale;

	for(int i = 0; i < m_particles.size(); ++i)
	{
		Fluid* p = &m_particles[i];

		float sum = 0.0;	
		m_neighborTable[i].clear();

		GridCellIndicies GC;
		m_grid.findCells(p->pos, radius, &GC);
		for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; cell++) 
		{
			int pndx = m_grid.getLastParticleIndex( GC.m_indicies[cell] );				
			while(pndx != INVALID_PARTICLE_INDEX) 
			{					
				Fluid* pcurr = &m_particles[pndx];					
				if(pcurr == p) 
				{
					pndx = pcurr->nextFluidIndex; 
					continue; 
				}
				
				float dx = (p->pos.x() - pcurr->pos.x()) * m_parameters.sph_simscale;		//Simulation-scale distance
				float dy = (p->pos.y() - pcurr->pos.y()) * m_parameters.sph_simscale;
				float dz = (p->pos.z() - pcurr->pos.z()) * m_parameters.sph_simscale;
				float dsq = dx*dx + dy*dy + dz*dz;
				if(m_parameters.m_R2 > dsq) 
				{
					float c = m_parameters.m_R2 - dsq;
					sum += c * c * c;
					
					if( !m_neighborTable[i].isFilled() ) m_neighborTable[i].addNeighbor( pndx, sqrt(dsq) );
				}
				
				pndx = pcurr->nextFluidIndex;
			}
		}
		p->density = sum * m_parameters.sph_pmass * m_parameters.m_Poly6Kern;	
		p->pressure = (p->density - m_parameters.sph_restdensity) * m_parameters.sph_intstiff;		
		p->density = 1.0f / p->density;		
	}
}

//Compute Forces - Very slow, but simple. O(n^2)
void FluidSystem::sph_computeForceSlow()
{
	double vterm = m_parameters.m_LapKern * m_parameters.sph_visc;

	for(int i = 0; i < m_particles.size(); ++i)
	{
		Fluid *p = &m_particles[i];

		Vector3DF force(0, 0, 0);
		
		for(int n = 0; n < m_particles.size(); ++n)
		{
			Fluid *q = &m_particles[n];
			if(p == q) continue;
			
			double dx = (p->pos.x() - q->pos.x()) * m_parameters.sph_simscale;			//Simulation-scale distance
			double dy = (p->pos.y() - q->pos.y()) * m_parameters.sph_simscale;
			double dz = (p->pos.z() - q->pos.z()) * m_parameters.sph_simscale;
			double dsq = dx*dx + dy*dy + dz*dz;
			if(m_parameters.m_R2 > dsq) 
			{
				double r = sqrt(dsq);
				double c = m_parameters.sph_smoothradius - r;
				double pterm = -0.5f * c * m_parameters.m_SpikyKern * (p->pressure + q->pressure) / r;
				double dterm = c * p->density * q->density;
				
				Vector3DF forceAdded( (pterm * dx + vterm * (q->vel_eval.x() - p->vel_eval.x())) * dterm,
									  (pterm * dy + vterm * (q->vel_eval.y() - p->vel_eval.y())) * dterm,
									  (pterm * dz + vterm * (q->vel_eval.z() - p->vel_eval.z())) * dterm );
				force += forceAdded;
			}
		}		
		
		p->sph_force = force;		
	}
}

//Compute Forces - Using spatial grid. Faster.
void FluidSystem::sph_computeForceGrid()
{
	BT_PROFILE("FluidSystem::sph_computeForceGrid()");

	float radius = m_parameters.sph_smoothradius / m_parameters.sph_simscale;
	double vterm =	m_parameters.m_LapKern * m_parameters.sph_visc;
		
	for(int i = 0; i < m_particles.size(); ++i)
	{
		Fluid *p = &m_particles[i];

		Vector3DF force(0, 0, 0);

		GridCellIndicies GC;
		m_grid.findCells(p->pos, radius, &GC);
		for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; cell++) 
		{
			int pndx = m_grid.getLastParticleIndex( GC.m_indicies[cell] );				
			while(pndx != INVALID_PARTICLE_INDEX) 
			{					
				Fluid *pcurr = &m_particles[pndx];
				if(pcurr == p) {pndx = pcurr->nextFluidIndex; continue; }
			
				double dx = (p->pos.x() - pcurr->pos.x()) * m_parameters.sph_simscale;		//Simulation-scale distance
				double dy = (p->pos.y() - pcurr->pos.y()) * m_parameters.sph_simscale;
				double dz = (p->pos.z() - pcurr->pos.z()) * m_parameters.sph_simscale;
				double dsq = dx*dx + dy*dy + dz*dz;
				if(m_parameters.m_R2 > dsq) 
				{
					double r = sqrt(dsq);
					double c = m_parameters.sph_smoothradius - r;
					double pterm = -0.5f * c * m_parameters.m_SpikyKern * ( p->pressure + pcurr->pressure) / r;
					double dterm = c * p->density * pcurr->density;
					
					Vector3DF forceAdded( (pterm * dx + vterm * (pcurr->vel_eval.x() - p->vel_eval.x())) * dterm,
										  (pterm * dy + vterm * (pcurr->vel_eval.y() - p->vel_eval.y())) * dterm,
										  (pterm * dz + vterm * (pcurr->vel_eval.z() - p->vel_eval.z())) * dterm );
					force += forceAdded;
				}
				pndx = pcurr->nextFluidIndex;
			}
		}
		p->sph_force = force;
	}
}

//Compute Forces - Using spatial grid with saved neighbor table. Fastest.
void FluidSystem::sph_computeForceGridNC()
{
	BT_PROFILE("FluidSystem::sph_computeForceGridNC()");

	float vterm = m_parameters.m_LapKern * m_parameters.sph_visc;
	
	for(int i = 0; i < m_particles.size(); ++i) 
	{
		Fluid *p = &m_particles[i];
		
		Vector3DF force(0, 0, 0);
		for(int j = 0; j < m_neighborTable[i].numNeighbors(); j++ ) 
		{
			Fluid *pcurr = &m_particles[ m_neighborTable[i].getNeighborIndex(j) ];
			float dx = (p->pos.x() - pcurr->pos.x()) * m_parameters.sph_simscale;		//Simulation-scale distance
			float dy = (p->pos.y() - pcurr->pos.y()) * m_parameters.sph_simscale;
			float dz = (p->pos.z() - pcurr->pos.z()) * m_parameters.sph_simscale;
			float c = m_parameters.sph_smoothradius - m_neighborTable[i].getDistance(j);
			float pterm = -0.5f * c * m_parameters.m_SpikyKern * ( p->pressure + pcurr->pressure) / m_neighborTable[i].getDistance(j);
			float dterm = c * p->density * pcurr->density;
			
			Vector3DF forceAdded( (pterm * dx + vterm * (pcurr->vel_eval.x() - p->vel_eval.x())) * dterm,
								  (pterm * dy + vterm * (pcurr->vel_eval.y() - p->vel_eval.y())) * dterm,
								  (pterm * dz + vterm * (pcurr->vel_eval.z() - p->vel_eval.z())) * dterm );
			force += forceAdded;
		}
		p->sph_force = force;
	}
}

float FluidSystem::getValue(float x, float y, float z) const
{
	float sum = 0.0;
	
	const float searchRadius = m_grid.getParameters().m_gridCellSize / 2.0f;
	const float R2 = 1.8*1.8;
	//const float R2 = 0.8*0.8;		//	marching cubes rendering test

	GridCellIndicies GC;
	m_grid.findCells( Vector3DF(x,y,z), searchRadius, &GC );
	
	for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; ++cell) 
	{
		int pndx = m_grid.getLastParticleIndex( GC.m_indicies[cell] );
		while(pndx != INVALID_PARTICLE_INDEX) 
		{
			const Fluid &F = m_particles[pndx];
			float dx = x - F.pos.x();
			float dy = y - F.pos.y();
			float dz = z - F.pos.z();
			float dsq = dx*dx+dy*dy+dz*dz;
				
			if(dsq < R2) sum += R2 / dsq;
			
			pndx = F.nextFluidIndex;
		}
	}
	
	return sum;	
}	
Vector3DF FluidSystem::getGradient(float x, float y, float z) const
{
	Vector3DF norm(0,0,0);
	
	const float searchRadius = m_grid.getParameters().m_gridCellSize / 2.0f;
	const float R2 = searchRadius*searchRadius;
	
	GridCellIndicies GC;
	m_grid.findCells( Vector3DF(x,y,z), searchRadius, &GC );
	for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; cell++)
	{
		int pndx = m_grid.getLastParticleIndex( GC.m_indicies[cell] );
		while ( pndx != INVALID_PARTICLE_INDEX ) 
		{					
			const Fluid &F = m_particles[pndx];
			float dx = x - F.pos.x();
			float dy = y - F.pos.y();
			float dz = z - F.pos.z();
			float dsq = dx*dx+dy*dy+dz*dz;
			
			if(0 < dsq && dsq < R2) 
			{
				dsq = 2.0*R2 / (dsq*dsq);
				
				Vector3DF particleNorm(dx * dsq, dy * dsq, dz * dsq);
				norm += particleNorm;
			}
			
			pndx = F.nextFluidIndex;
		}
	}
	
	norm.normalize();
	return norm;
}

void FluidSystem::grid_insertParticles()
{	
	BT_PROFILE("FluidSystem::grid_insertParticles()");

	//Reset particles
	for(int i = 0; i < m_particles.size(); ++i) m_particles[i].nextFluidIndex = INVALID_PARTICLE_INDEX;	

	//Reset grid
	m_grid.clear();
	
	//Insert particles into grid
	for(int particleIndex = 0; particleIndex < m_particles.size(); ++particleIndex) 
	{
		Fluid *p = &m_particles[particleIndex];
		m_grid.insertParticle(p, particleIndex);
	}
}


////////////////////////////////////////////////////////////////////////////////
/// struct FluidEmitter
////////////////////////////////////////////////////////////////////////////////
void FluidEmitter::emit(FluidSystem *fluidSystem, int numParticles, float spacing)
{	
	const float DEGtoRAD = 3.141592/180.0;
	
	int x = static_cast<int>( sqrt(static_cast<float>(numParticles)) );
	
	for(int i = 0; i < numParticles; i++) 
	{
		float ang_rand = ( static_cast<float>(rand()*2.0/RAND_MAX) - 1.0 ) * m_yawSpread;
		float tilt_rand = ( static_cast<float>(rand()*2.0/RAND_MAX) - 1.0 ) * m_pitchSpread;
		
		//Modified - set y to vertical axis
		Vector3DF dir( 	cos((m_yaw + ang_rand) * DEGtoRAD) * sin((m_pitch + tilt_rand) * DEGtoRAD) * m_velocity,
						cos((m_pitch + tilt_rand) * DEGtoRAD) * m_velocity,
						sin((m_yaw + ang_rand) * DEGtoRAD) * sin((m_pitch + tilt_rand) * DEGtoRAD) * m_velocity );
		
		Vector3DF pos( spacing*(i/x), spacing*(i%x), 0 );
		pos += m_position;
		
		Fluid *f = fluidSystem->addFluid(pos);
		f->vel = dir;
		f->vel_eval = dir;
	}
}

void FluidEmitter::addVolume(FluidSystem *fluidSystem, const Vector3DF &min, const Vector3DF &max, float spacing)
{
	Vector3DF d(max.x() - min.x(), max.y() - min.y(), max.z() - min.z());
	
	for(float z = max.z(); z >= min.z(); z -= spacing) 
		for(float y = min.y(); y <= max.y(); y += spacing) 
			for(float x = min.x(); x <= max.x(); x += spacing) 
			{
				Fluid *f = fluidSystem->addFluid( Vector3DF(x,y,z) );
			}
}


