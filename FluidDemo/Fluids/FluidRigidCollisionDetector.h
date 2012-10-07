/* 	FluidRigidCollisionDetector.h
	Copyright (C) 2012 Jackson Lee

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
#ifndef FLUID_RIGID_COLLISION_DETECTOR_H
#define FLUID_RIGID_COLLISION_DETECTOR_H

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"

#include "BulletCollision/CollisionDispatch/btCollisionWorld.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"

#include "FluidSph.h"

class btCollisionWorld;

///Describes a contact between a FluidSph particle and a btCollisionObject or btRigidBody.
struct FluidRigidContact
{
	int m_fluidParticleIndex;
	
	btCollisionObject *m_object;
	btVector3 m_normalOnObject;
	btVector3 m_hitPointWorldOnObject;
	btScalar m_distance;
};

///Contains all FluidRigidContact for a single FluidSph. 
struct FluidRigidContactGroup
{
	FluidSph *m_fluid;
	btAlignedObjectArray<FluidRigidContact> m_contacts;
	
	FluidRigidContactGroup(FluidSph *fluid) : m_fluid(fluid) {}
};

///Detects and stores collisions between FluidSph and btCollisionObject/btRigidBody.
class FluidRigidCollisionDetector
{
	btAlignedObjectArray<FluidRigidContactGroup> m_contactGroups;

public:
	void detectCollisions(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids, btCollisionWorld *world)
	{
		BT_PROFILE("FluidRigidCollisionDetector::detectCollisions()");
		
		m_contactGroups.clear();
	
		for(int i = 0; i < fluids->size(); ++i)
			detectCollisionsSingleFluid( FG, (*fluids)[i], world );
	}
	
	const btAlignedObjectArray<FluidRigidContactGroup>& getContactGroups() const { return m_contactGroups; }
	btAlignedObjectArray<FluidRigidContactGroup>& internalGetContactGroups() { return m_contactGroups; }
	
private:
	//Collide individual fluid particles against several btCollsionObject(s) using Bullet's broadphase
	void detectCollisionsSingleFluid(const FluidParametersGlobal &FG, FluidSph *fluid, btCollisionWorld *world);
	
	//Collide individual btCollisionObjects against several fluid particles using FluidSortingGrid broadphase
	void detectCollisionsSingleFluid2(const FluidParametersGlobal &FG, FluidSph *fluid, btCollisionWorld *world);
	
	//Experimental
	void detectCollisionsSingleFluidCcd(const FluidParametersGlobal &FG, FluidSph *fluid, btCollisionWorld *world);
};


#endif
