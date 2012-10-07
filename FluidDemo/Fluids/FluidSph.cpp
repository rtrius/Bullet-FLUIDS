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
#include "LinearMath/btAabbUtil2.h"		//TestPointAgainstAabb2()

#include "FluidSortingGrid.h"

FluidSph::FluidSph(const FluidParametersGlobal &FG, const btVector3 &volumeMin, const btVector3 &volumeMax, int maxNumParticles)
{
	setMaxParticles(maxNumParticles);
	configureGridAndAabb(FG, volumeMin, volumeMax);
}

void FluidSph::configureGridAndAabb(const FluidParametersGlobal &FG, const btVector3 &volumeMin, const btVector3 &volumeMax)
{
	m_localParameters.m_volumeMin = volumeMin;
	m_localParameters.m_volumeMax = volumeMax;

	m_grid.setup(FG.m_simulationScale, FG.m_sphSmoothRadius);
}
void FluidSph::getCurrentAabb(const FluidParametersGlobal &FG, btVector3 *out_min, btVector3 *out_max) const
{
	m_grid.getPointAabb(out_min, out_max);

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
	
	m_grid.clear();
}

btScalar FluidSph::getValue(btScalar x, btScalar y, btScalar z) const
{
	btScalar sum = 0.0;
	
	const btScalar searchRadius = m_grid.getCellSize() / btScalar(2.0);
	const btScalar R2 = btScalar(1.8) * btScalar(1.8);
	
	FluidSortingGrid::FoundCells foundCells;
	m_grid.findCells( btVector3(x,y,z), searchRadius, &foundCells );
		
	for(int cell = 0; cell < FluidSortingGrid::NUM_FOUND_CELLS; cell++) 
	{
		FluidGridIterator &FI = foundCells.m_iterators[cell];
			
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
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
	
	const btScalar searchRadius = m_grid.getCellSize() / btScalar(2.0);
	const btScalar R2 = searchRadius*searchRadius;
	
	FluidSortingGrid::FoundCells foundCells;
	m_grid.findCells( btVector3(x,y,z), searchRadius, &foundCells );
	
	for(int cell = 0; cell < FluidSortingGrid::NUM_FOUND_CELLS; cell++)
	{
		FluidGridIterator &FI = foundCells.m_iterators[cell];
			
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
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

//Assumes that out_unique is already sorted.
//Removes duplicates; rearranges array such that all unique values are in the range [0, uniqueSize).
void makeUnique(btAlignedObjectArray<int> *out_unique)
{
	int uniqueSize = 0;
	if( out_unique->size() ) 
	{
		uniqueSize = 1;
		for(int i = 1; i < out_unique->size(); ++i)
		{
			if( (*out_unique)[i] != (*out_unique)[i-1] )
			{
				(*out_unique)[uniqueSize] = (*out_unique)[i];
				++uniqueSize;
			}
		}
	}
	
	out_unique->resize(uniqueSize);
}
struct AscendingSortPredicate { inline bool operator() (const int &a, const int &b) const { return (a < b); } };
void FluidSph::removeMarkedParticles()
{
	//makeUnique() assumes that the array is sorted
	//m_removedFluidIndicies.heapSort( AscendingSortPredicate() );
	m_removedFluidIndicies.quickSort( AscendingSortPredicate() );
	
	//Remove duplicate indicies
	makeUnique(&m_removedFluidIndicies);
	
	//Since removing elements from the array invalidates(higher) indicies,
	//elements should be removed in descending order.
	for(int i = m_removedFluidIndicies.size() - 1; i >= 0; --i) m_particles.removeParticle( m_removedFluidIndicies[i] );
	
	m_removedFluidIndicies.resize(0);
}

void FluidSph::insertParticlesIntoGrid()
{	
	BT_PROFILE("FluidSph::insertParticlesIntoGrid()");
	
	//
	m_grid.clear();
	m_grid.insertParticles(&m_particles);
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


// /////////////////////////////////////////////////////////////////////////////
// struct FluidAbsorber
// /////////////////////////////////////////////////////////////////////////////
void FluidAbsorber::absorb(FluidSph *fluid)
{
	const FluidSortingGrid &grid = fluid->getGrid();
	
	btAlignedObjectArray<int> gridCellIndicies;
	grid.getGridCellIndiciesInAabb(m_min, m_max, &gridCellIndicies);
	
	for(int i = 0; i < gridCellIndicies.size(); ++i)
	{
		FluidGridIterator FI = grid.getGridCell( gridCellIndicies[i] );
		
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
		{
			if( TestPointAgainstAabb2( m_min, m_max, fluid->getPosition(n) ) ) fluid->markParticleForRemoval(n);
		}
	}
}

