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
#include "btSphereCapsuleCollisionAlgorithm.h"

#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"

#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btCapsuleShape.h"

#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"

btSphereCapsuleCollisionAlgorithm::btSphereCapsuleCollisionAlgorithm(btPersistentManifold* mf, const btCollisionAlgorithmConstructionInfo& ci,
															const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap,
															bool swapped)
: btActivatingCollisionAlgorithm(ci, col0Wrap, col1Wrap), m_ownManifold(false), m_manifoldPtr(mf), m_swapped(swapped)
{
	if(!m_manifoldPtr)
	{
		const btCollisionObjectWrapper* sphereWrap = (m_swapped) ? col1Wrap : col0Wrap;
		const btCollisionObjectWrapper* capsuleWrap = (m_swapped) ? col0Wrap : col1Wrap;
	
		m_manifoldPtr = m_dispatcher->getNewManifold( sphereWrap->getCollisionObject(), capsuleWrap->getCollisionObject() );
		m_ownManifold = true;
	}
}

void btSphereCapsuleCollisionAlgorithm::processCollision(const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap,
														const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	if(!m_manifoldPtr) return;

	resultOut->setPersistentManifold(m_manifoldPtr);
	
	const btCollisionObjectWrapper* sphereWrap = (m_swapped) ? col1Wrap : col0Wrap;
	const btCollisionObjectWrapper* capsuleWrap = (m_swapped) ? col0Wrap : col1Wrap;

	const btSphereShape* sphereShape = static_cast<const btSphereShape*>( sphereWrap->getCollisionShape() );
	const btCapsuleShape* capsuleShape = static_cast<const btCapsuleShape*>( capsuleWrap->getCollisionShape() );
	
	//Consider the sphere as a point, and the capsule as a constant distance(radius) from a line segment(2 points).
	//Add the sphere radius to the capsule radius.
	{
		const btTransform& capsuleTransform = capsuleWrap->getWorldTransform();
		
		btVector3 capsulePointA(0, 0, 0);
		capsulePointA.m_floats[ capsuleShape->getUpAxis() ] = capsuleShape->getHalfHeight();	//getUpAxis() == 0 for X, 1 for Y, 2 for Z
		
		btVector3 capsulePointB = -capsulePointA;
	
		//Determine sphere position relative to capsule
		btVector3 sphereRelPos = capsuleTransform.invXform( sphereWrap->getWorldTransform().getOrigin() );
	
		//Determine the point on the line segment that is nearest to the sphere center
		//Line equation: capsulePointA + capsuleAToB * t, where (0 <= t <= 1)
		btVector3 capsuleAToB = capsulePointB - capsulePointA;
		
		btScalar t = -(capsulePointA - sphereRelPos).dot(capsuleAToB) / capsuleAToB.length2();	
		t = btMin( btMax(btScalar(0.0), t), btScalar(1.0) );
		
		btVector3 nearestPointOnLine = capsulePointA + capsuleAToB * t;		//Nearest point to sphere on line segment
		btVector3 lineToSphere = sphereRelPos - nearestPointOnLine;
		
		//Calculate distance
		btScalar lineToSphereDistance = lineToSphere.length();
		btScalar distance = lineToSphereDistance - ( capsuleShape->getRadius() + sphereShape->getRadius() );
		
		//Negative distance indicates collision/penetration
		if( distance < btScalar(0.0) ) 
		{
			btVector3 normalOnCapsule = (lineToSphereDistance > SIMD_EPSILON) ? lineToSphere / lineToSphereDistance : btVector3(0, 1, 0);
			btVector3 pointOnCapsule = nearestPointOnLine + normalOnCapsule * capsuleShape->getRadius();
			
			//Convect normal, point into world space
			btVector3 pointOnCapsuleInWorld = capsuleTransform(pointOnCapsule);
			btVector3 normalOnCapsuleInWorld = capsuleTransform.getBasis() * normalOnCapsule;

			resultOut->addContactPoint(normalOnCapsuleInWorld, pointOnCapsuleInWorld, distance);
		}
	}

	if(m_ownManifold) resultOut->refreshContactPoints();
}

btScalar btSphereCapsuleCollisionAlgorithm::calculateTimeOfImpact(btCollisionObject* col0, btCollisionObject* col1,
													const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	(void)resultOut;
	(void)dispatchInfo;
	(void)col0;
	(void)col1;

	//Not implemented
	return btScalar(1.0);
}
