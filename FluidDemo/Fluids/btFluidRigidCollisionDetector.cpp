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


struct btParticleContactResultMulti : public btCollisionWorld::ContactResultCallback
{
	static const int MAX_COLLISIONS = 4;	//Max number of collisions with btCollisionObject(s)

	btFluidRigidContact m_contacts[MAX_COLLISIONS];
	int m_numCollisions;
	
	btCollisionObject* m_particleObject;

	btParticleContactResultMulti(btCollisionObject* particleObject, int particleIndex) : m_particleObject(particleObject), m_numCollisions(0)
	{
		for(int i = 0; i < MAX_COLLISIONS; ++i) m_contacts[i].m_fluidParticleIndex = particleIndex; 
	}
	
	virtual	btScalar addSingleResult(btManifoldPoint& cp, const btCollisionObjectWrapper* colObj0Wrap, int partId0, int index0,
														  const btCollisionObjectWrapper* colObj1Wrap, int partId1, int index1)
	{
		const btCollisionObject* colObj0 = colObj0Wrap->m_collisionObject;
		const btCollisionObject* colObj1 = colObj1Wrap->m_collisionObject;
		
		if(m_numCollisions < MAX_COLLISIONS)
		{
			//Assume 0 == A, 1 == B
			if(m_particleObject == colObj0 && colObj1)
			{
				m_contacts[m_numCollisions].m_object = const_cast<btCollisionObject*>(colObj1);
				m_contacts[m_numCollisions].m_normalOnObject = cp.m_normalWorldOnB;
				m_contacts[m_numCollisions].m_hitPointWorldOnObject = cp.getPositionWorldOnB();
				m_contacts[m_numCollisions].m_distance = cp.getDistance();
				++m_numCollisions;
			}
			else if(m_particleObject == colObj1 && colObj0)
			{
				//	this branch is never reached?
				
				m_contacts[m_numCollisions].m_object = const_cast<btCollisionObject*>(colObj0);
				m_contacts[m_numCollisions].m_normalOnObject = -cp.m_normalWorldOnB;
				m_contacts[m_numCollisions].m_hitPointWorldOnObject = cp.getPositionWorldOnA();
				m_contacts[m_numCollisions].m_distance = cp.getDistance();
				++m_numCollisions;
			}
		}
		
		//Value returned from btCollisionWorld::ContactResultCallback::addSingleResult() appears to be unused
		const btScalar UNUSED = btScalar(1.0);
		return UNUSED;
	}
};
void btFluidRigidCollisionDetector::detectCollisionsSingleFluid(const btFluidParametersGlobal& FG, btFluidSph* fluid, btCollisionWorld* world)
{
	//m_particleRadius is at simscale; divide by simscale to transform into world scale
	btSphereShape particleShape( FG.m_particleRadius / FG.m_simulationScale );
	
	btCollisionObject particleObject;
	particleObject.setCollisionShape(&particleShape);
		
	btTransform& particleTransform = particleObject.getWorldTransform();
	particleTransform.setRotation( btQuaternion::getIdentity() );
	
	
	btFluidRigidContactGroup resultContacts(fluid);
	for(int i = 0; i < fluid->numParticles(); ++i)
	{
		particleTransform.setOrigin( fluid->getPosition(i) );
			
		btParticleContactResultMulti result(&particleObject, i);
		world->contactTest(&particleObject, result);
			
		for(int n = 0; n < result.m_numCollisions; ++n)
			resultContacts.m_contacts.push_back(result.m_contacts[n]);
	}
	
	if( resultContacts.m_contacts.size() ) m_contactGroups.push_back(resultContacts);
}

struct btParticleContactResult : public btCollisionWorld::ContactResultCallback
{
	btFluidRigidContact m_contact;
	btCollisionObject* m_particleObject;

	btParticleContactResult(btCollisionObject* particleObject, int particleIndex) : m_particleObject(particleObject)
	{
		m_contact.m_fluidParticleIndex = particleIndex;
		m_contact.m_object = 0; 
	}
	
	virtual	btScalar addSingleResult(btManifoldPoint& cp, const btCollisionObjectWrapper* colObj0Wrap, int partId0, int index0,
														  const btCollisionObjectWrapper* colObj1Wrap, int partId1, int index1)
	{
		const btCollisionObject* colObj0 = colObj0Wrap->m_collisionObject;
		const btCollisionObject* colObj1 = colObj1Wrap->m_collisionObject;
		
		m_contact.m_distance = cp.getDistance();
	
		//Assume 0 == A, 1 == B
		if(m_particleObject == colObj0 && colObj1)		
		{
			m_contact.m_object = const_cast<btCollisionObject*>(colObj1);
			m_contact.m_normalOnObject = cp.m_normalWorldOnB;
			m_contact.m_hitPointWorldOnObject = cp.getPositionWorldOnB();
		}
		else if(m_particleObject == colObj1 && colObj0)
		{
			//	this branch is never reached?
		
			m_contact.m_object = const_cast<btCollisionObject*>(colObj0);
			m_contact.m_normalOnObject = -cp.m_normalWorldOnB;
			m_contact.m_hitPointWorldOnObject = cp.getPositionWorldOnA();
		}
		
		//Value returned from btCollisionWorld::ContactResultCallback::addSingleResult() appears to be unused
		const btScalar UNUSED = btScalar(1.0);
		return UNUSED;
	}

	const bool hasHit() const { return (m_contact.m_object != 0); }
};
struct btAabbTestCallback : public btBroadphaseAabbCallback
{
	btAlignedObjectArray<btCollisionObject*> m_results;

	virtual bool process(const btBroadphaseProxy* proxy) 
	{ 
		m_results.push_back( static_cast<btCollisionObject*>(proxy->m_clientObject) ); 
		
		//Return value of btBroadphaseAabbCallback::process() appears to be unused
		const bool UNUSED = true;
		return UNUSED;
	}
};

void btFluidRigidCollisionDetector::detectCollisionsSingleFluid2(const btFluidParametersGlobal& FG, btFluidSph* fluid, btCollisionWorld* world)
{
	//Sample AABB from all fluid particles, detect intersecting 
	//broadphase proxies, and test each btCollisionObject
	//against all fluid grid cells intersecting with the btCollisionObject's AABB.
	btVector3 fluidSystemMin, fluidSystemMax;
	fluid->getCurrentAabb(FG, &fluidSystemMin, &fluidSystemMax);
	
	//m_particleRadius is at simscale; divide by simscale to transform into world scale
	btScalar particleRadius = FG.m_particleRadius / FG.m_simulationScale;
	btSphereShape particleShape(particleRadius);
	
	btCollisionObject particleObject;
	particleObject.setCollisionShape(&particleShape);
					
	btTransform& particleTransform = particleObject.getWorldTransform();
	particleTransform.setRotation( btQuaternion::getIdentity() );
	
	//
	const btFluidSortingGrid& grid = fluid->getGrid();
	
	btAabbTestCallback aabbTest;
	world->getBroadphase()->aabbTest(fluidSystemMin, fluidSystemMax, aabbTest);
	
	btFluidRigidContactGroup resultContacts(fluid);
	for(int i = 0; i < aabbTest.m_results.size(); ++i)
	{
		btCollisionObject* const object = aabbTest.m_results[i];
	
		const btVector3& objectMin = object->getBroadphaseHandle()->m_aabbMin;
		const btVector3& objectMax = object->getBroadphaseHandle()->m_aabbMax;
		
		btAlignedObjectArray<int> gridCellIndicies;
		grid.getGridCellIndiciesInAabb(objectMin, objectMax, &gridCellIndicies);
		
		for(int j = 0; j < gridCellIndicies.size(); ++j)
		{
			btFluidGridIterator FI = grid.getGridCell( gridCellIndicies[j] );
			
			for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
			{
				const btVector3& fluidPos = fluid->getPosition(n);
				
				btVector3 fluidMin( fluidPos.x() - particleRadius, fluidPos.y() - particleRadius, fluidPos.z() - particleRadius );
				btVector3 fluidMax( fluidPos.x() + particleRadius, fluidPos.y() + particleRadius, fluidPos.z() + particleRadius );
						
				if( fluidMin.x() <= objectMax.x() && objectMin.x() <= fluidMax.x()  
				 && fluidMin.y() <= objectMax.y() && objectMin.y() <= fluidMax.y()
				 && fluidMin.z() <= objectMax.z() && objectMin.z() <= fluidMax.z() )
				{
					particleTransform.setOrigin(fluidPos);
					
					btParticleContactResult result(&particleObject, n);
					world->contactPairTest(&particleObject, const_cast<btCollisionObject*>(object), result);
					
					if( result.hasHit() ) resultContacts.m_contacts.push_back(result.m_contact);
				}
			}
		}
	}
	
	if( resultContacts.m_contacts.size() ) m_contactGroups.push_back(resultContacts);
}

void btFluidRigidCollisionDetector::detectCollisionsSingleFluidCcd(const btFluidParametersGlobal& FG, btFluidSph* fluid, btCollisionWorld* world)
{
	//from detectCollisionsSingleFluid()
	//m_particleRadius is at simscale; divide by simscale to transform into world scale
	btSphereShape particleShape( FG.m_particleRadius / FG.m_simulationScale );
	
	btCollisionObject particleObject;
	particleObject.setCollisionShape(&particleShape);
		
	btTransform& particleTransform = particleObject.getWorldTransform();
	particleTransform.setRotation( btQuaternion::getIdentity() );
	//from detectCollisionsSingleFluid()
	
	//Multiply fluid->getVelocity() with velocityScale to get the distance moved in the last frame
	const btScalar velocityScale = FG.m_timeStep / FG.m_simulationScale;
	
	//int numCollided = 0;
	
	btFluidRigidContactGroup resultContacts(fluid);
	for(int i = 0; i < fluid->numParticles(); ++i)
	{
		const btVector3& pos = fluid->getPosition(i);
		btVector3 prev_pos = fluid->getPosition(i) - fluid->getVelocity(i)*velocityScale;
		
		//	investigate:
		//	Calling btCollisionWorld::rayTest() with the ray's origin inside a btCollisionObject
		//	reports no collisions; is this expected behavior when using btCollisionWorld::ClosestRayResultCallback?
		//btCollisionWorld::ClosestRayResultCallback result( btVector3(0.f, -1.0f, 0.f), btVector3(0.f, -2.0f, 0.f) );
		btCollisionWorld::ClosestRayResultCallback result(prev_pos, pos);
		result.m_collisionFilterGroup = btBroadphaseProxy::DefaultFilter;
		result.m_collisionFilterMask = btBroadphaseProxy::AllFilter;
		
		world->rayTest(result.m_rayFromWorld, result.m_rayToWorld, result);
		
		//printf( "result.m_closestHitFraction, result.hasHit(): %f, %d \n", result.m_closestHitFraction, result.hasHit() );
		
		if( result.hasHit() )
		{
			btScalar distanceMoved = result.m_rayFromWorld.distance(result.m_rayToWorld);
			btScalar distanceCollided = result.m_rayFromWorld.distance(result.m_hitPointWorld);
	
			btScalar distance = distanceCollided - distanceMoved;
		
			//++numCollided;
			btFluidRigidContact contact;
			contact.m_fluidParticleIndex = i;
			contact.m_object = const_cast<btCollisionObject*>(result.m_collisionObject);
			contact.m_normalOnObject = result.m_hitNormalWorld;
			contact.m_hitPointWorldOnObject = result.m_hitPointWorld;
			contact.m_distance = distance;
			
			resultContacts.m_contacts.push_back(contact);
		}
		else
		{
			//Since btCollisionWorld::rayTest() does not report a collision if the ray's origin
			//is inside a btCollsionObject, use btCollisionWorld::contactTest() in case
			//the particle is inside.
			
			//from detectCollisionsSingleFluid()
			particleTransform.setOrigin( fluid->getPosition(i) );
					
			btParticleContactResultMulti result(&particleObject, i);
			world->contactTest(&particleObject, result);
					
			for(int n = 0; n < result.m_numCollisions; ++n)
				resultContacts.m_contacts.push_back(result.m_contacts[n]);
			//from detectCollisionsSingleFluid()
		}
	}
	
	if( resultContacts.m_contacts.size() ) m_contactGroups.push_back(resultContacts);
	//printf("numCollided: %d \n", numCollided);
}	

