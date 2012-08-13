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

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

#include "FluidStaticGrid.h"
#include "FluidSortingGrid.h"

FluidSph::FluidSph(const FluidParametersGlobal &FG, const btVector3 &volumeMin, const btVector3 &volumeMax, 
				   FluidGrid::Type gridType, int maxNumParticles)
{
	m_grid = 0;

	setMaxParticles(maxNumParticles);
	configureGridAndAabb(FG, volumeMin, volumeMax, gridType);
}

void FluidSph::configureGridAndAabb(const FluidParametersGlobal &FG, const btVector3 &volumeMin, const btVector3 &volumeMax, 
									FluidGrid::Type gridType)
{
	m_localParameters.m_volumeMin = volumeMin;
	m_localParameters.m_volumeMax = volumeMax;

	btScalar simCellSize = FG.m_sphSmoothRadius * btScalar(2.0);	//Grid cell size (2r)
	
	if(m_grid)delete m_grid;
	
	if(gridType == FluidGrid::FT_IndexRange) m_grid = new FluidSortingGrid(FG.m_simulationScale, simCellSize);
	else m_grid = new FluidStaticGrid( volumeMin, volumeMax, FG.m_simulationScale, simCellSize, btScalar(1.0) );
}
void FluidSph::getCurrentAabb(const FluidParametersGlobal &FG, btVector3 *out_min, btVector3 *out_max) const
{
	m_grid->getPointAabb(out_min, out_max);

	btScalar particleRadius = FG.m_particleRadius / FG.m_simulationScale;
	btVector3 radius(particleRadius, particleRadius, particleRadius);
	
	*out_min -= radius;
	*out_max += radius;
}

void FluidSph::setMaxParticles(int maxNumParticles)
{
	if( maxNumParticles < m_particles.size() )m_particles.resize(maxNumParticles);
	m_particles.setMaxParticles(maxNumParticles);
}

void FluidSph::removeAllParticles()
{
	m_particles.resize(0);
	
	m_removedFluidIndicies.resize(0);
	
	m_grid->clear();
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

btScalar FluidSph::getValue(btScalar x, btScalar y, btScalar z) const
{
	btScalar sum = 0.0;
	
	const bool isLinkedList = m_grid->isLinkedListGrid();
	const btScalar searchRadius = m_grid->getCellSize() / btScalar(2.0);
	const btScalar R2 = btScalar(1.8) * btScalar(1.8);
	//const btScalar R2 = btScalar(0.8) * btScalar(0.8);		//	marching cubes rendering test
	
	FluidGrid::FoundCells foundCells;
	m_grid->findCells( btVector3(x,y,z), searchRadius, &foundCells );
		
	for(int cell = 0; cell < FluidGrid::NUM_FOUND_CELLS; cell++) 
	{
		FluidGridIterator &FI = foundCells.m_iterators[cell];
			
		for( int n = FI.m_firstIndex; FluidGridIterator::isIndexValid(n, FI.m_lastIndex); 
				 n = FluidGridIterator::getNextIndex(n, isLinkedList, m_particles.m_nextFluidIndex) )
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
	
	const bool isLinkedList = m_grid->isLinkedListGrid();
	const btScalar searchRadius = m_grid->getCellSize() / btScalar(2.0);
	const btScalar R2 = searchRadius*searchRadius;
	
	FluidGrid::FoundCells foundCells;
	m_grid->findCells( btVector3(x,y,z), searchRadius, &foundCells );
	
	for(int cell = 0; cell < FluidGrid::NUM_FOUND_CELLS; cell++)
	{
		FluidGridIterator &FI = foundCells.m_iterators[cell];
			
		for( int n = FI.m_firstIndex; FluidGridIterator::isIndexValid(n, FI.m_lastIndex); 
				 n = FluidGridIterator::getNextIndex(n, isLinkedList, m_particles.m_nextFluidIndex) )
		{
			const btVector3 &position = m_particles.m_pos[n];
			btScalar dx = x - position.x();
			btScalar dy = y - position.y();
			btScalar dz = z - position.z();
			btScalar dsq = dx*dx+dy*dy+dz*dz;
			
			if(0 < dsq && dsq < R2) 
			{
				dsq = btScalar(2.0)*R2 / (dsq*dsq);
				
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
void FluidSph::removeMarkedParticles()
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

void FluidSph::insertParticlesIntoGrid()
{	
	BT_PROFILE("FluidSph::insertParticlesIntoGrid()");

	//Reset particles
	for(int i = 0; i < numParticles(); ++i) m_particles.m_nextFluidIndex[i] = INVALID_PARTICLE_INDEX;
	
	//
	m_grid->clear();
	m_grid->insertParticles(&m_particles);
}


// /////////////////////////////////////////////////////////////////////////////
// struct FluidEmitter
// /////////////////////////////////////////////////////////////////////////////
void FluidEmitter::emit(FluidSph *fluid, int numParticles, btScalar spacing)
{
	int x = static_cast<int>( btSqrt(static_cast<btScalar>(numParticles)) );
	
	for(int i = 0; i < numParticles; i++) 
	{
		btScalar ang_rand = ( static_cast<btScalar>(rand()*btScalar(2.0)/RAND_MAX) - btScalar(1.0) ) * m_yawSpread;
		btScalar tilt_rand = ( static_cast<btScalar>(rand()*btScalar(2.0)/RAND_MAX) - btScalar(1.0) ) * m_pitchSpread;
		
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


