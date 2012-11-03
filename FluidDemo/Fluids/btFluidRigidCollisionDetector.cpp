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
#include "btFluidRigidCollisionDetector.h"

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btAabbUtil2.h"		//TestAabbAgainstAabb2()
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionDispatch/btManifoldResult.h"
#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"
#include "BulletCollision/BroadphaseCollision/btDispatcher.h"

#include "btFluidSph.h"


struct btFluidRigidContactResult : public btManifoldResult
{
	const btCollisionObject* m_particleObject;
	int m_fluidParticleIndex;

	btFluidRigidContactGroup& m_rigidContactGroup;

	btFluidRigidContactResult(const btCollisionObjectWrapper* obj0Wrap, const btCollisionObjectWrapper* obj1Wrap,
								btFluidRigidContactGroup& rigidContactGroup, 
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
		
		btFluidRigidContact m_contact;
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

void btFluidRigidCollisionDetector::detectCollisionsSingleFluid(btDispatcher* dispatcher, const btDispatcherInfo& dispatchInfo, btFluidSph* fluid)
{
	BT_PROFILE("detectCollisionsSingleFluid()");
	
	const btFluidParametersLocal& FL = fluid->getLocalParameters();
	const btFluidSortingGrid& grid = fluid->getGrid();
	btAlignedObjectArray<btFluidRigidContactGroup>& rigidContacts = fluid->internalGetRigidContacts();
	
	btSphereShape particleShape(FL.m_particleRadius);
	
	btCollisionObject particleObject;
	particleObject.setCollisionShape(&particleShape);
		
	btTransform& particleTransform = particleObject.getWorldTransform();
	particleTransform.setRotation( btQuaternion::getIdentity() );
	
	btAlignedObjectArray<int> gridCellIndicies;
	
	const btAlignedObjectArray<const btCollisionObject*>& intersectingRigidAabbs = fluid->internalGetIntersectingRigidAabbs();
	for(int i = 0; i < intersectingRigidAabbs.size(); ++i)
	{
		gridCellIndicies.clear();
	
		const btCollisionObject* rigidObject = intersectingRigidAabbs[i];
		btCollisionObjectWrapper rigidWrap( 0, rigidObject->getCollisionShape(), rigidObject, rigidObject->getWorldTransform() );
		btFluidRigidContactGroup contactGroup;
		contactGroup.m_object = rigidObject;
		
		const btVector3& rigidMin = rigidObject->getBroadphaseHandle()->m_aabbMin;
		const btVector3& rigidMax = rigidObject->getBroadphaseHandle()->m_aabbMax;
		
		grid.getGridCellIndiciesInAabb(rigidMin, rigidMax, &gridCellIndicies);
		for(int j = 0; j < gridCellIndicies.size(); ++j)
		{
			btFluidGridIterator FI = grid.getGridCell( gridCellIndicies[j] );
			
			for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
			{
				const btVector3& fluidPos = fluid->getPosition(n);
				
				btVector3 fluidMin( fluidPos.x() - FL.m_particleRadius, fluidPos.y() - FL.m_particleRadius, fluidPos.z() - FL.m_particleRadius );
				btVector3 fluidMax( fluidPos.x() + FL.m_particleRadius, fluidPos.y() + FL.m_particleRadius, fluidPos.z() + FL.m_particleRadius );
				
				if( TestAabbAgainstAabb2(fluidMin, fluidMax, rigidMin, rigidMax) )
				{
					particleTransform.setOrigin(fluidPos);
					
					btCollisionObjectWrapper particleWrap( 0, particleObject.getCollisionShape(), 
															&particleObject, particleObject.getWorldTransform() );

					btCollisionAlgorithm* algorithm = dispatcher->findAlgorithm(&particleWrap, &rigidWrap);
					if(algorithm)
					{
						btFluidRigidContactResult result(&particleWrap, &rigidWrap, contactGroup, &particleObject, n);
						algorithm->processCollision(&particleWrap, &rigidWrap, dispatchInfo, &result);

						algorithm->~btCollisionAlgorithm();
						dispatcher->freeCollisionAlgorithm(algorithm);
					}
					
				}
			}
		}
		
		if( contactGroup.numContacts() ) rigidContacts.push_back(contactGroup);
	}
}
