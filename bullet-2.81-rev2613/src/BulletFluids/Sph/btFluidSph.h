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

#ifndef BT_FLUID_SPH_H
#define BT_FLUID_SPH_H

#include "BulletCollision/CollisionDispatch/btCollisionObject.h"

#include "btFluidParticles.h"
#include "btFluidSphParameters.h"
#include "btFluidSortingGrid.h"

///Describes a single contact between a btFluidSph particle and a btCollisionObject or btRigidBody.
struct btFluidSphRigidContact
{
	int m_fluidParticleIndex;
	
	btVector3 m_normalOnObject;
	btVector3 m_hitPointWorldOnObject;
	btScalar m_distance;
};

///Contains all btFluidSphRigidContact between a btFluidSph and a btCollisionObject.
struct btFluidSphRigidContactGroup
{
	const btCollisionObject* m_object;
	btAlignedObjectArray<btFluidSphRigidContact> m_contacts;
	
	void addContact(const btFluidSphRigidContact &contact) { m_contacts.push_back(contact); }
	int numContacts() const { return m_contacts.size(); }
};

class btFluidSphSolver;

///@brief Main fluid class. Coordinates a set of btFluidParticles with material definition and grid broadphase.
class btFluidSph : public btCollisionObject
{
protected:
	btFluidSphParametersLocal	m_localParameters;
	
	btFluidSortingGrid		m_grid;
	
	btFluidParticles 		m_particles;
	
	btAlignedObjectArray<int> m_removedFluidIndicies;

	btAlignedObjectArray<const btCollisionObject*> m_intersectingRigidAabb;	///<Contains btCollisionObject/btRigidBody(not btSoftbody)
	btAlignedObjectArray<btFluidSphRigidContactGroup> m_rigidContacts;
	
	//If either override is set, the fluid is passed separately to the solver(btFluidSph-btFluidSph interaction is disabled)
	btFluidSphSolver* m_overrideSolver;
	btFluidSphParametersGlobal* m_overrideParameters;
	
public:
	///@param FG Reference returned by btFluidRigidDynamicsWorld::getGlobalParameters().
	btFluidSph(const btFluidSphParametersGlobal& FG, int maxNumParticles);
	virtual ~btFluidSph();
	
	int	numParticles() const { return m_particles.size(); }
	int getMaxParticles() const { return m_particles.getMaxParticles(); }
	void setMaxParticles(int maxNumParticles);	///<Removes particles if( maxNumParticles < numParticles() ).
	
	///Returns a particle index; creates a new particle if numParticles() < getMaxParticles(), returns numParticles() otherwise.
	int addParticle(const btVector3& position) { return m_particles.addParticle(position); }
	
	///Duplicate indicies are ignored, so a particle may be marked twice without any issues.
	void markParticleForRemoval(int index) { m_removedFluidIndicies.push_back(index); }
	
	void removeAllParticles();
	void removeMarkedParticles();	///<Automatically called during btFluidRigidDynamicsWorld::stepSimulation(); invalidates grid.
	void insertParticlesIntoGrid(); ///<Automatically called during btFluidRigidDynamicsWorld::stepSimulation(); updates the grid.
	
	///Avoid placing particles at the same position; particles with same position and velocity will experience identical SPH forces and not seperate.
	void setPosition(int index, const btVector3& position) { m_particles.m_pos[index] = position; }
	
	///Sets both velocities; getVelocity() and getEvalVelocity().
	void setVelocity(int index, const btVector3& velocity) 
	{
		m_particles.m_vel[index] = velocity;
		m_particles.m_vel_eval[index] = velocity;
	}
	
	///Accumulates a simulation scale force that is applied, and then set to 0 during btFluidRigidDynamicsWorld::stepSimulation().
	void applyForce(int index, const btVector3& force) { m_particles.m_accumulatedForce[index] += force; }
	
	const btVector3& getPosition(int index) const { return m_particles.m_pos[index]; }
	const btVector3& getVelocity(int index) const { return m_particles.m_vel[index]; }			///<Returns the 'current+(1/2)*timestep' velocity.
	const btVector3& getEvalVelocity(int index) const { return m_particles.m_vel_eval[index]; } ///<Returns the current velocity.
	
	void setParticleUserPointer(int index, void* userPointer) { m_particles.m_userPointer[index] = userPointer; }
	void* getParticleUserPointer(int index) const { return m_particles.m_userPointer[index]; }
	//
	const btFluidSortingGrid& getGrid() const { return m_grid; }
	
	///@param FG Reference returned by btFluidRigidDynamicsWorld::getGlobalParameters().
	void setGridCellSize(const btFluidSphParametersGlobal& FG);
	
	//Parameters
	const btFluidSphParametersLocal& getLocalParameters() const { return m_localParameters; }
	btFluidSphParametersLocal& getLocalParameters() { return m_localParameters; }
	void setLocalParameters(const btFluidSphParametersLocal& FP) { m_localParameters = FP; }
	btScalar getEmitterSpacing(const btFluidSphParametersGlobal& FG) const { return m_localParameters.m_particleDist / FG.m_simulationScale; }
	
	///If solver is not 0, then it will be used instead of the solver specified by btFluidRigidDynamicsWorld::getFluidSolver()
	void setOverrideSolver(btFluidSphSolver* solver) { m_overrideSolver = solver; }
	btFluidSphSolver* getOverrideSolver() const { return m_overrideSolver; }
	
	///If parameters is not 0, then it will be used instead of the parameters specified by btFluidRigidDynamicsWorld::getGlobalParameters()
	void setOverrideParameters(btFluidSphParametersGlobal* parameters) { m_overrideParameters = parameters; }
	btFluidSphParametersGlobal* getOverrideParameters() const { return m_overrideParameters; }
	
	//Metablobs	
	btScalar getValue(btScalar x, btScalar y, btScalar z) const;
	btVector3 getGradient(btScalar x, btScalar y, btScalar z) const;

	const btFluidParticles& getParticles() const { return m_particles; }
	btFluidParticles& internalGetParticles() { return m_particles; }
	btFluidSortingGrid& internalGetGrid() { return m_grid; }
	
	//btFluidSph-Rigid collisions
	void internalClearRigidContacts()
	{
		m_intersectingRigidAabb.clear();
		m_rigidContacts.clear();
	}
	btAlignedObjectArray<const btCollisionObject*>& internalGetIntersectingRigidAabbs() { return m_intersectingRigidAabb; }
	btAlignedObjectArray<btFluidSphRigidContactGroup>& internalGetRigidContacts() { return m_rigidContacts; }
	
	const btAlignedObjectArray<const btCollisionObject*>& getIntersectingRigidAabbs() const { return m_intersectingRigidAabb; }
	const btAlignedObjectArray<btFluidSphRigidContactGroup>& getRigidContacts() const { return m_rigidContacts; }
	
	//btCollisionObject
	virtual void setCollisionShape(btCollisionShape *collisionShape) { btAssert(0); }
	
	virtual void getAabb(btVector3& aabbMin, btVector3& aabbMax) const
	{
		m_grid.getPointAabb(aabbMin, aabbMax);
		
		btScalar radius = m_localParameters.m_particleRadius;
		btVector3 extent(radius, radius, radius);
		
		aabbMin -= extent;
		aabbMax += extent;
	}

	static const btFluidSph* upcast(const btCollisionObject* colObj)
	{
		// replace later with CO_FLUID_SPH
		return (colObj->getInternalType() == CO_USER_TYPE) ? (const btFluidSph*)colObj : 0;
	}
	static btFluidSph* upcast(btCollisionObject* colObj)
	{
		// replace later with CO_FLUID_SPH
		return (colObj->getInternalType() == CO_USER_TYPE) ? (btFluidSph*)colObj : 0;
	}
};

///@brief Adds particles to a btFluidSph.
struct btFluidEmitter
{
	btVector3 m_position;

	btScalar m_velocity;
	
	btScalar m_yaw;
	btScalar m_pitch;
	
	btScalar m_yawSpread;
	btScalar m_pitchSpread;
	
	bool m_useRandomIfAllParticlesAllocated;
	
	btFluidEmitter() : m_position(0,0,0), m_yaw(0), m_pitch(0), 
					 m_velocity(0), m_yawSpread(0), m_pitchSpread(0),
					 m_useRandomIfAllParticlesAllocated(true) {}
	
	void emit(btFluidSph* fluid, int numParticles, btScalar spacing);

	static void addVolume(btFluidSph* fluid, const btVector3& min, const btVector3& max, btScalar spacing);
};

///@brief Marks particles from a btFluidSph for removal; see btFluidSph::removeMarkedParticles().
struct btFluidAbsorber
{
	btVector3 m_min;
	btVector3 m_max;
	
	//int m_maxParticlesRemoved;
	//	add velocity limit / max particles removed, etc.?
	
	btFluidAbsorber() : m_min(0,0,0), m_max(0,0,0) {}
	
	void absorb(btFluidSph* fluid);
};

#endif


