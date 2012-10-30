/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. 
   If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_FLUID_RIGID_COLLISION_DETECTOR_H
#define BT_FLUID_RIGID_COLLISION_DETECTOR_H

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"

#include "BulletCollision/CollisionDispatch/btCollisionWorld.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"

#include "btFluidSph.h"

class btCollisionWorld;

///Describes a contact between a btFluidSph particle and a btCollisionObject or btRigidBody.
struct btFluidRigidContact
{
	int m_fluidParticleIndex;
	
	btCollisionObject* m_object;
	btVector3 m_normalOnObject;
	btVector3 m_hitPointWorldOnObject;
	btScalar m_distance;
};

///Contains all btFluidRigidContact for a single btFluidSph. 
struct btFluidRigidContactGroup
{
	btFluidSph* m_fluid;
	btAlignedObjectArray<btFluidRigidContact> m_contacts;
	
	btFluidRigidContactGroup(btFluidSph* fluid) : m_fluid(fluid) {}
};

///Detects and stores collisions between btFluidSph and btCollisionObject/btRigidBody.
class btFluidRigidCollisionDetector
{
	btAlignedObjectArray<btFluidRigidContactGroup> m_contactGroups;

public:
	void detectCollisions(const btFluidParametersGlobal& FG, btAlignedObjectArray<btFluidSph*>* fluids, btCollisionWorld* world)
	{
		BT_PROFILE("btFluidRigidCollisionDetector::detectCollisions()");
		
		m_contactGroups.clear();
	
		for(int i = 0; i < fluids->size(); ++i)
			detectCollisionsSingleFluid( FG, (*fluids)[i], world );
	}
	
	const btAlignedObjectArray<btFluidRigidContactGroup>& getContactGroups() const { return m_contactGroups; }
	btAlignedObjectArray<btFluidRigidContactGroup>& internalGetContactGroups() { return m_contactGroups; }
	
private:
	//Collide individual fluid particles against several btCollsionObject(s) using Bullet's broadphase
	void detectCollisionsSingleFluid(const btFluidParametersGlobal& FG, btFluidSph* fluid, btCollisionWorld* world);
	
	//Collide individual btCollisionObjects against several fluid particles using btFluidSortingGrid broadphase
	void detectCollisionsSingleFluid2(const btFluidParametersGlobal& FG, btFluidSph* fluid, btCollisionWorld* world);
	
	//Experimental
	void detectCollisionsSingleFluidCcd(const btFluidParametersGlobal& FG, btFluidSph* fluid, btCollisionWorld* world);
};


#endif
