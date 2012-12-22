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
#include "btFluidSphRigidCollisionDetector.h"

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btAabbUtil2.h"		//TestPointAgainstAabb2()
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionDispatch/btManifoldResult.h"
#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"
#include "BulletCollision/BroadphaseCollision/btDispatcher.h"

#include "btFluidSph.h"


struct btFluidSphRigidContactResult : public btManifoldResult
{
	const btCollisionObject* m_particleObject;
	int m_fluidParticleIndex;

	btFluidSphRigidContactGroup& m_rigidContactGroup;

	btFluidSphRigidContactResult(const btCollisionObjectWrapper* obj0Wrap, const btCollisionObjectWrapper* obj1Wrap,
								btFluidSphRigidContactGroup& rigidContactGroup, 
								const btCollisionObject* particleObject, int particleIndex)
	: btManifoldResult(obj0Wrap, obj1Wrap), m_rigidContactGroup(rigidContactGroup), 
	m_particleObject(particleObject), m_fluidParticleIndex(particleIndex) {}

	virtual void addContactPoint(const btVector3& normalOnBInWorld, const btVector3& pointBInWorld, btScalar distance)
	{
		bool isSwapped = m_manifoldPtr->getBody0() != m_body0Wrap->getCollisionObject();
		
		const btCollisionObjectWrapper* obj0Wrap = (isSwapped) ? m_body1Wrap : m_body0Wrap;
		const btCollisionObjectWrapper* obj1Wrap = (isSwapped) ? m_body0Wrap : m_body1Wrap;
		
		const btCollisionObject* colObj0 = obj0Wrap->m_collisionObject;
		const btCollisionObject* colObj1 = obj1Wrap->m_collisionObject;
		
		btFluidSphRigidContact m_contact;
		m_contact.m_fluidParticleIndex = m_fluidParticleIndex;
		m_contact.m_distance = distance;
		
		if(m_particleObject == colObj0 && colObj1)		
		{
			m_contact.m_normalOnObject = normalOnBInWorld;
			m_contact.m_hitPointWorldOnObject = pointBInWorld;
		}
		else if(m_particleObject == colObj1 && colObj0)
		{
			//This branch may never be reached
			btVector3 pointAInWorld = pointBInWorld + normalOnBInWorld * distance;
		
			m_contact.m_normalOnObject = -normalOnBInWorld;
			m_contact.m_hitPointWorldOnObject = pointAInWorld;
		}
		
		m_rigidContactGroup.addContact(m_contact);
	}
};



///Adds contacts to a btFluidSphRigidContactGroup
struct btFluidSphRigidNarrowphaseCallback : public btFluidSortingGrid::AabbCallback
{
	//All members must be set
	btDispatcher* m_dispatcher;
	const btDispatcherInfo* m_dispatchInfo;
	
	btFluidSphRigidContactGroup* m_contactGroup;
	
	const btFluidSph* m_fluid;
	
	btCollisionObject* m_particleObject;
	const btCollisionObject* m_rigidObject;
	
	btVector3 m_expandedRigidAabbMin;
	btVector3 m_expandedRigidAabbMax;

	btFluidSphRigidNarrowphaseCallback() {}
	
	virtual bool processParticles(const btFluidGridIterator FI, const btVector3& aabbMin, const btVector3& aabbMax)
	{
		btTransform& particleTransform = m_particleObject->getWorldTransform();
			
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
		{
			const btVector3& fluidPos = m_fluid->getPosition(n);
			
			if( TestPointAgainstAabb2(m_expandedRigidAabbMin, m_expandedRigidAabbMax, fluidPos) )
			{
				particleTransform.setOrigin(fluidPos);
				
				btCollisionObjectWrapper particleWrap( 0, m_particleObject->getCollisionShape(), m_particleObject, particleTransform );
				btCollisionObjectWrapper rigidWrap( 0, m_rigidObject->getCollisionShape(), m_rigidObject, m_rigidObject->getWorldTransform() );
				
				btCollisionAlgorithm* algorithm = m_dispatcher->findAlgorithm(&particleWrap, &rigidWrap);
				if(algorithm)
				{
					btFluidSphRigidContactResult result(&particleWrap, &rigidWrap, *m_contactGroup, m_particleObject, n);
					
					{
						//BT_PROFILE("algorithm->processCollision()");
						algorithm->processCollision(&particleWrap, &rigidWrap, *m_dispatchInfo, &result);
					}
					algorithm->~btCollisionAlgorithm();
					m_dispatcher->freeCollisionAlgorithm(algorithm);
				}
					
			}
		}
		return true;
	}
};

void btFluidSphRigidCollisionDetector::performNarrowphase(btDispatcher* dispatcher, const btDispatcherInfo& dispatchInfo, btFluidSph* fluid)
{
	BT_PROFILE("FluidSphRigid - performNarrowphase()");
	
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	const btFluidSortingGrid& grid = fluid->getGrid();
	btAlignedObjectArray<btFluidSphRigidContactGroup>& rigidContacts = fluid->internalGetRigidContacts();
	
	btSphereShape particleShape(FL.m_particleRadius);
	
	btCollisionObject particleObject;
	particleObject.setCollisionShape(&particleShape);
		
	btTransform& particleTransform = particleObject.getWorldTransform();
	particleTransform.setRotation( btQuaternion::getIdentity() );
	
	const btAlignedObjectArray<const btCollisionObject*>& intersectingRigidAabbs = fluid->internalGetIntersectingRigidAabbs();
	for(int i = 0; i < intersectingRigidAabbs.size(); ++i)
	{
		const btCollisionObject* rigidObject = intersectingRigidAabbs[i];
		
		btFluidSphRigidContactGroup contactGroup;
		contactGroup.m_object = rigidObject;
		
		//Add particle radius to rigid AABB to avoid calculating particle AABB; use point-AABB test instead of AABB-AABB
		const btVector3 fluidRadius(FL.m_particleRadius, FL.m_particleRadius, FL.m_particleRadius);
		const btVector3 rigidMin = rigidObject->getBroadphaseHandle()->m_aabbMin - fluidRadius;
		const btVector3 rigidMax = rigidObject->getBroadphaseHandle()->m_aabbMax + fluidRadius;
		
		btFluidSphRigidNarrowphaseCallback particleRigidCollider;
		particleRigidCollider.m_dispatcher = dispatcher;
		particleRigidCollider.m_dispatchInfo = &dispatchInfo;
		particleRigidCollider.m_contactGroup = &contactGroup;
		particleRigidCollider.m_fluid = fluid;
		particleRigidCollider.m_particleObject = &particleObject;
		particleRigidCollider.m_rigidObject = rigidObject;
		particleRigidCollider.m_expandedRigidAabbMin = rigidMin;
		particleRigidCollider.m_expandedRigidAabbMax = rigidMax;
		
		grid.forEachGridCell(rigidMin, rigidMax, particleRigidCollider);
		
		if( contactGroup.numContacts() ) rigidContacts.push_back(contactGroup);
	}
}
