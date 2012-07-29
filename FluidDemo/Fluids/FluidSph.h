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
#ifndef FLUID_SPH_H
#define FLUID_SPH_H

#include "FluidParticles.h"
#include "FluidParameters.h"
#include "FluidGrid.h"

class FluidSph
{
	FluidParametersLocal	m_localParameters;
	
	FluidGrid				*m_grid;
	
	FluidParticles 			m_particles;
	
	btAlignedObjectArray<int> m_removedFluidIndicies;

public:
	FluidSph(const FluidParametersGlobal &FG, const btVector3 &volumeMin, const btVector3 &volumeMax, 
			 FluidGridType gridType, int maxNumParticles);
	~FluidSph() { if(m_grid) delete m_grid; }
	
	int	numParticles() const { return m_particles.size(); }
	void setMaxParticles(int maxNumParticles);	//Removes particles if( maxNumParticles < numParticles() )
	int addParticle(const btVector3 &position) { return m_particles.addParticle(position); }
	
	///Removal occurs on the next call to stepSimulation()
	void markParticleForRemoval(int index) { m_removedFluidIndicies.push_back(index); }
	
	void removeAllParticles();
	void removeMarkedParticles();	//Automatically called during FluidWorld::stepSimulation(); invalidates grid
	void insertParticlesIntoGrid();
	
	//
	void setPosition(int index, const btVector3 &position) { m_particles.m_pos[index] = position; }
	void setVelocity(int index, const btVector3 &velocity)
	{
		m_particles.m_vel[index] = velocity;
		m_particles.m_vel_eval[index] = velocity;
	}
	void applyAcceleration(int index, const btVector3 &acceleration) { m_particles.m_externalAcceleration[index] += acceleration; }
	
	const btVector3& getPosition(int index) const { return m_particles.m_pos[index]; }
	const btVector3& getVelocity(int index) const { return m_particles.m_vel[index]; }
	const btVector3& getEvalVelocity(int index) const { return m_particles.m_vel_eval[index]; }
	
	//
	const FluidGrid* getGrid() const { return m_grid; }
	const btAlignedObjectArray<int>& getNextFluidIndicies() const { return m_particles.m_nextFluidIndex; }
	void configureGridAndAabb(const FluidParametersGlobal &FG, const btVector3 &volumeMin, const btVector3 &volumeMax, FluidGridType gridType);
	
	//Parameters
	const FluidParametersLocal& getLocalParameters() const { return m_localParameters; }
	void setLocalParameters(const FluidParametersLocal &FP) { m_localParameters = FP; }
	btScalar getEmitterSpacing(const FluidParametersGlobal &FG) const { return m_localParameters.m_particleDist / FG.sph_simscale; }
	
	//Metablobs	
	btScalar getValue(btScalar x, btScalar y, btScalar z) const;
	btVector3 getGradient(btScalar x, btScalar y, btScalar z) const;

	///Internal functions; do not use.
	FluidParticles& internalGetFluidParticles() { return m_particles; }
	FluidGrid* internalGetGrid() { return m_grid; }
};

struct FluidEmitter
{
	btVector3 m_position;

	btScalar m_velocity;
	
	btScalar m_yaw;
	btScalar m_pitch;
	
	btScalar m_yawSpread;
	btScalar m_pitchSpread;
	
	FluidEmitter() : m_position(0,0,0), m_yaw(0), m_pitch(0), 
					 m_velocity(0), m_yawSpread(0), m_pitchSpread(0) {}
	
	void emit(FluidSph *fluid, int numParticles, btScalar spacing);

	static void addVolume(FluidSph *fluid, const btVector3 &min, const btVector3 &max, btScalar spacing);
};

inline bool isInsideAabb(const btVector3 &min, const btVector3 &max, const btVector3 &point)
{
	return (   min.x() <= point.x() && point.x() <= max.x()
			&& min.y() <= point.y() && point.y() <= max.y()
			&& min.z() <= point.z() && point.z() <= max.z() );
}
struct FluidAbsorber
{
	btVector3 m_min;
	btVector3 m_max;
	
	//int m_maxParticlesRemoved;
	//	add velocity limit / max particles removed, etc.?
	
	FluidAbsorber() : m_min(0,0,0), m_max(0,0,0) {}
	
	void absorb(FluidSph *fluid)
	{
		const FluidGrid *grid = fluid->getGrid();
		const bool isLinkedList = (grid->getGridType() == FT_LinkedList);
		
		btAlignedObjectArray<int> gridCellIndicies;
		grid->getGridCellIndiciesInAabb(m_min, m_max, &gridCellIndicies);
		
		for(int i = 0; i < gridCellIndicies.size(); ++i)
		{
			FluidGridIterator FI = grid->getGridCell( gridCellIndicies[i] );
			
			for( int n = FI.m_firstIndex; FluidGridIterator::isIndexValid(n, FI.m_lastIndex);
					 n = FluidGridIterator::getNextIndex(n, isLinkedList, fluid->getNextFluidIndicies()) )
			{
				if( isInsideAabb( m_min, m_max, fluid->getPosition(n) ) ) fluid->markParticleForRemoval(n);
			}
		}
	}
};

#endif


