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

#include "btFluidSph.h"

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro
#include "LinearMath/btAabbUtil2.h"		//TestPointAgainstAabb2()
#include "LinearMath/btRandom.h"		//GEN_rand(), GEN_RAND_MAX

#include "btFluidSortingGrid.h"
#include "btFluidSphCollisionShape.h"

btFluidSph::btFluidSph(const btFluidSphParametersGlobal& FG, int maxNumParticles)
{
	m_overrideSolver = 0;
	m_overrideParameters = 0;

	setMaxParticles(maxNumParticles);
	setGridCellSize(FG);
	
	//btCollisionObject
	{
		m_worldTransform.setIdentity();
		m_internalType = CO_USER_TYPE;	// replace later with CO_FLUID_SPH
		
		void* ptr = btAlignedAlloc( sizeof(btFluidSphCollisionShape), 16 );
		m_collisionShape = new(ptr) btFluidSphCollisionShape(this);
		m_collisionShape->setMargin( btScalar(0.25) );	//Arbitrary value
		
		m_rootCollisionShape = m_collisionShape;
	}
}
btFluidSph::~btFluidSph()
{
	//btCollisionObject
	{
		m_collisionShape->~btCollisionShape();
		btAlignedFree(m_collisionShape);
	}
}

void btFluidSph::setGridCellSize(const btFluidSphParametersGlobal& FG)
{
	m_grid.setCellSize(FG.m_simulationScale, FG.m_sphSmoothRadius);
}

void btFluidSph::setMaxParticles(int maxNumParticles)
{
	if( maxNumParticles < m_particles.size() )m_particles.resize(maxNumParticles);
	m_particles.setMaxParticles(maxNumParticles);
}

void btFluidSph::removeAllParticles()
{
	m_particles.resize(0);
	
	m_removedFluidIndicies.resize(0);
	
	m_grid.clear();
}

btScalar btFluidSph::getCombinedPosition(btScalar x, btScalar y, btScalar z) const
{
	const btScalar worldSphRadius = m_grid.getCellSize();	//Grid cell size == sph interaction radius, at world scale
	const btScalar R2 = worldSphRadius * worldSphRadius;
	
	btFluidSortingGrid::FoundCells foundCells;
	m_grid.findCells( btVector3(x,y,z), foundCells );
		
	btScalar sum = 0.0;
	for(int cell = 0; cell < btFluidSortingGrid::NUM_FOUND_CELLS; cell++) 
	{
		btFluidGridIterator& FI = foundCells.m_iterators[cell];
			
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
		{
			const btVector3& position = m_particles.m_pos[n];
			btScalar dx = x - position.x();
			btScalar dy = y - position.y();
			btScalar dz = z - position.z();
			btScalar distanceSquared = dx*dx + dy*dy + dz*dz;
				
			if(distanceSquared < R2) sum += R2 / distanceSquared;
		}
	}
	
	return sum;
}	
btVector3 btFluidSph::getGradient(btScalar x, btScalar y, btScalar z) const
{
	const btScalar worldSphRadius = m_grid.getCellSize();	//Grid cell size == sph interaction radius, at world scale
	const btScalar R2 = worldSphRadius*worldSphRadius;
	
	btFluidSortingGrid::FoundCells foundCells;
	m_grid.findCells( btVector3(x,y,z), foundCells );
	
	btVector3 normal(0,0,0);
	for(int cell = 0; cell < btFluidSortingGrid::NUM_FOUND_CELLS; cell++)
	{
		btFluidGridIterator& FI = foundCells.m_iterators[cell];
			
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
		{
			const btVector3& position = m_particles.m_pos[n];
			btScalar dx = x - position.x();
			btScalar dy = y - position.y();
			btScalar dz = z - position.z();
			btScalar distanceSquared = dx*dx + dy*dy + dz*dz;
			
			if( btScalar(0.0) < distanceSquared && distanceSquared < R2 ) 
			{
				distanceSquared = btScalar(2.0)*R2 / (distanceSquared*distanceSquared);
				
				btVector3 particleNorm(dx * distanceSquared, dy * distanceSquared, dz * distanceSquared);
				normal += particleNorm;
			}
		}
	}
	
	normal.normalize();
	return normal;
}

//Assumes that out_unique is already sorted.
//Removes duplicates; rearranges array such that all unique values are in the range [0, uniqueSize).
void makeUniqueInt(btAlignedObjectArray<int>& out_unique)
{
	int uniqueSize = 0;
	if( out_unique.size() ) 
	{
		uniqueSize = 1;
		for(int i = 1; i < out_unique.size(); ++i)
		{
			if( out_unique[i] != out_unique[i-1] )
			{
				out_unique[uniqueSize] = out_unique[i];
				++uniqueSize;
			}
		}
	}
	
	out_unique.resize(uniqueSize);
}
struct AscendingSortPredicate { inline bool operator() (const int& a, const int& b) const { return (a < b); } };
void btFluidSph::removeMarkedParticles()
{
	//makeUnique() assumes that the array is sorted
	//m_removedFluidIndicies.heapSort( AscendingSortPredicate() );
	m_removedFluidIndicies.quickSort( AscendingSortPredicate() );
	
	//Remove duplicate indicies
	makeUniqueInt(m_removedFluidIndicies);
	
	//Since removing elements from the array invalidates(higher) indicies,
	//elements should be removed in descending order.
	for(int i = m_removedFluidIndicies.size() - 1; i >= 0; --i) m_particles.removeParticle( m_removedFluidIndicies[i] );
	
	m_removedFluidIndicies.resize(0);
}

void btFluidSph::insertParticlesIntoGrid()
{	
	BT_PROFILE("btFluidSph::insertParticlesIntoGrid()");
	
	//
	m_grid.clear();
	m_grid.insertParticles(m_particles);
}


// /////////////////////////////////////////////////////////////////////////////
// struct btFluidEmitter
// /////////////////////////////////////////////////////////////////////////////
void btFluidEmitter::emit(btFluidSph* fluid, int numParticles, btScalar spacing)
{
	int x = static_cast<int>( btSqrt(static_cast<btScalar>(numParticles)) );
	
	for(int i = 0; i < numParticles; i++) 
	{
		btScalar ang_rand = ( static_cast<btScalar>(GEN_rand()*btScalar(2.0)/GEN_RAND_MAX) - btScalar(1.0) ) * m_yawSpread;
		btScalar tilt_rand = ( static_cast<btScalar>(GEN_rand()*btScalar(2.0)/GEN_RAND_MAX) - btScalar(1.0) ) * m_pitchSpread;
		
		btVector3 dir( 	btCos((m_yaw + ang_rand) * SIMD_RADS_PER_DEG) * btSin((m_pitch + tilt_rand) * SIMD_RADS_PER_DEG) * m_velocity,
						btCos((m_pitch + tilt_rand) * SIMD_RADS_PER_DEG) * m_velocity,
						btSin((m_yaw + ang_rand) * SIMD_RADS_PER_DEG) * btSin((m_pitch + tilt_rand) * SIMD_RADS_PER_DEG) * m_velocity );
		
		btVector3 position( spacing*(i/x), spacing*(i%x), 0 );
		position += m_position;
		
		int index = fluid->addParticle(position);
		
		if( index != fluid->numParticles() ) fluid->setVelocity(index, dir);
		else if(m_useRandomIfAllParticlesAllocated)
		{
			index = ( fluid->numParticles() - 1 ) * GEN_rand() / GEN_RAND_MAX;		//Random index
		
			fluid->setPosition(index, position);
			fluid->setVelocity(index, dir);
		}
	}
}
void btFluidEmitter::addVolume(btFluidSph* fluid, const btVector3& min, const btVector3& max, btScalar spacing)
{
	for(btScalar z = max.z(); z >= min.z(); z -= spacing) 
		for(btScalar y = min.y(); y <= max.y(); y += spacing) 
			for(btScalar x = min.x(); x <= max.x(); x += spacing) 
			{
				fluid->addParticle( btVector3(x,y,z) );
			}
}


// /////////////////////////////////////////////////////////////////////////////
// struct btFluidAbsorber
// /////////////////////////////////////////////////////////////////////////////
struct btFluidAbsorberCallback : public btFluidSortingGrid::AabbCallback
{
	btFluidSph* m_fluidSph;

	const btVector3& m_min;
	const btVector3& m_max;

	btFluidAbsorberCallback(btFluidSph* fluidSph, const btVector3& min, const btVector3& max) 
	: m_fluidSph(fluidSph), m_min(min), m_max(max) {}
	
	virtual bool processParticles(const btFluidGridIterator FI, const btVector3& aabbMin, const btVector3& aabbMax)
	{
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
		{
			if( TestPointAgainstAabb2( m_min, m_max, m_fluidSph->getPosition(n) ) ) m_fluidSph->markParticleForRemoval(n);
		}
	
		return true;
	}
};
void btFluidAbsorber::absorb(btFluidSph* fluid)
{
	const btFluidSortingGrid& grid = fluid->getGrid();
	
	btFluidAbsorberCallback absorber(fluid, m_min, m_max);
	grid.forEachGridCell(m_min, m_max, absorber);
}
