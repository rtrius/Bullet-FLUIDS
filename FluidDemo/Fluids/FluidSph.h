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

///@brief Main fluid class. Coordinates a set of FluidParticles, a FluidParametersLocal, and a FluidGrid.
class FluidSph
{
	FluidParametersLocal	m_localParameters;
	
	FluidGrid				*m_grid;
	
	FluidParticles 			m_particles;
	
	btAlignedObjectArray<int> m_removedFluidIndicies;

public:
	///See FluidSph::configureGridAndAabb().
	FluidSph(const FluidParametersGlobal &FG, const btVector3 &volumeMin, const btVector3 &volumeMax, 
			 FluidGrid::Type gridType, int maxNumParticles);
	~FluidSph();
	
	int	numParticles() const { return m_particles.size(); }
	int getMaxParticles() const { return m_particles.getMaxParticles(); }
	void setMaxParticles(int maxNumParticles);	///<Removes particles if( maxNumParticles < numParticles() ).
	
	///Returns a particle index; creates a new particle if numParticles() < getMaxParticles(), returns a random index otherwise.
	int addParticle(const btVector3 &position) { return m_particles.addParticle(position); }
	
	///Duplicate indicies are ignored, so a particle may be marked twice without any issues.
	void markParticleForRemoval(int index) { m_removedFluidIndicies.push_back(index); }
	
	void removeAllParticles();
	void removeMarkedParticles();	///<Automatically called during FluidWorld::stepSimulation(); invalidates grid.
	void insertParticlesIntoGrid(); ///<Automatically called during FluidWorld::stepSimulation(); updates the grid.
	
	//
	void setPosition(int index, const btVector3 &position) { m_particles.m_pos[index] = position; }
	
	///Sets both velocities; getVelocity() and getEvalVelocity().
	void setVelocity(int index, const btVector3 &velocity) 
	{
		m_particles.m_vel[index] = velocity;
		m_particles.m_vel_eval[index] = velocity;
	}
	
	///Accumulates a simulation scale acceleration that is applied, and then set to 0 during FluidWorld::stepSimulation().
	void applyAcceleration(int index, const btVector3 &acceleration) { m_particles.m_externalAcceleration[index] += acceleration; }
	
	const btVector3& getPosition(int index) const { return m_particles.m_pos[index]; }
	const btVector3& getVelocity(int index) const { return m_particles.m_vel[index]; }			///<Returns m_vel of FluidParticles.
	const btVector3& getEvalVelocity(int index) const { return m_particles.m_vel_eval[index]; } ///<Returns m_vel_eval of FluidParticles.
	
	//
	const FluidGrid* getGrid() const { return m_grid; }
	const btAlignedObjectArray<int>& getNextFluidIndicies() const { return m_particles.m_nextFluidIndex; } ///<See FluidGridIterator.
	
	///@param FG Reference returned by FluidWorld::getGlobalParameters().
	///@param volumeMin, volumeMax AABB defining the extent to which particles may move.
	///@param gridType FluidGrid::FT_LinkedList for FluidStaticGrid, FluidGrid::FT_IndexRange for FluidSortingGrid.
	void configureGridAndAabb(const FluidParametersGlobal &FG, const btVector3 &volumeMin, const btVector3 &volumeMax, FluidGrid::Type gridType);
	void getCurrentAabb(const FluidParametersGlobal &FG, btVector3 *out_min, btVector3 *out_max) const;
	
	
	//Parameters
	const FluidParametersLocal& getLocalParameters() const { return m_localParameters; }
	void setLocalParameters(const FluidParametersLocal &FP) { m_localParameters = FP; }
	btScalar getEmitterSpacing(const FluidParametersGlobal &FG) const { return m_localParameters.m_particleDist / FG.m_simulationScale; }
	
	//Metablobs	
	btScalar getValue(btScalar x, btScalar y, btScalar z) const;
	btVector3 getGradient(btScalar x, btScalar y, btScalar z) const;

	FluidParticles& internalGetFluidParticles() { return m_particles; }
	FluidGrid* internalGetGrid() { return m_grid; }
};

///@brief Adds particles to a FluidSph.
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

///@brief Marks particles from a FluidSph for removal; see FluidSph::removeMarkedParticles().
struct FluidAbsorber
{
	btVector3 m_min;
	btVector3 m_max;
	
	//int m_maxParticlesRemoved;
	//	add velocity limit / max particles removed, etc.?
	
	FluidAbsorber() : m_min(0,0,0), m_max(0,0,0) {}
	
	void absorb(FluidSph *fluid);
};

#endif


