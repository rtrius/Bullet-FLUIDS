/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2009 Erwin Coumans  http://bulletphysics.com

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.

Experimental Buoyancy fluid demo written by John McCutchan
*/
#include "btHfFluidRigidCollisionAlgorithm.h"

#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btBoxShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"

#include "btHfFluidBuoyantConvexShape.h"
#include "btHfFluid.h"


btHfFluidRigidCollisionAlgorithm::btHfFluidRigidCollisionAlgorithm(const btCollisionAlgorithmConstructionInfo& ci, 
																	const btCollisionObjectWrapper* body0Wrap,
																	const btCollisionObjectWrapper* body1Wrap, bool isSwapped)
: btCollisionAlgorithm(ci), m_isSwapped(isSwapped), 
	m_convexTrianglecallback(ci.m_dispatcher1, body0Wrap, body1Wrap, !isSwapped) // we flip the isSwapped because we are hf fluid vs. convex and callback expects convex vs. concave
{
}

void btHfFluidRigidCollisionAlgorithm::processGround (const btCollisionObjectWrapper* hfFluidWrap, const btCollisionObjectWrapper* rigidWrap,
														const btDispatcherInfo& dispatchInfo,btManifoldResult* resultOut)
{
	// to perform the convex shape vs. ground terrain:
	// we pull the convex shape out of the btHfFluidBuoyantConvexShape and replace it temporarily
	const btHfFluidBuoyantConvexShape* tmpRigidShape = static_cast<const btHfFluidBuoyantConvexShape*>( rigidWrap->getCollisionShape() );
	const btConvexShape* convexShape = tmpRigidShape->getConvexShape();
	
	btCollisionObjectWrapper tempRigidWrap( rigidWrap, convexShape, rigidWrap->getCollisionObject(), rigidWrap->getWorldTransform() );
	
	btScalar triangleMargin = m_rigidCollisionObject->getCollisionShape()->getMargin();
	resultOut->setPersistentManifold(m_convexTrianglecallback.m_manifoldPtr);
	
	m_convexTrianglecallback.setTimeStepAndCounters(triangleMargin, dispatchInfo, &tempRigidWrap, hfFluidWrap, resultOut);
	m_hfFluid->foreachGroundTriangle(&m_convexTrianglecallback, m_convexTrianglecallback.getAabbMin(), m_convexTrianglecallback.getAabbMax());
	
	resultOut->refreshContactPoints();
	m_convexTrianglecallback.clearWrapperData();
}

btScalar btHfFluidRigidCollisionAlgorithm::processFluid (const btDispatcherInfo& dispatchInfo, btScalar density, btScalar floatyness)
{
	btRigidBody* rb = btRigidBody::upcast(m_rigidCollisionObject);
	btHfFluidColumnRigidBodyCallback columnCallback (rb, dispatchInfo.m_debugDraw, density, floatyness);
	m_hfFluid->foreachFluidColumn (&columnCallback, m_convexTrianglecallback.getAabbMin(), m_convexTrianglecallback.getAabbMax());
	return columnCallback.getVolume ();
}

void btHfFluidRigidCollisionAlgorithm::applyFluidFriction (btScalar mu, btScalar submerged_percentage)
{
	btRigidBody* rb = btRigidBody::upcast(m_rigidCollisionObject);
	btScalar dt = btScalar(1.0f/60.0f);

//#define OLD_WAY
#ifdef OLD_WAY
	btScalar radius = btScalar(0.0f);
	{
		btVector3 aabbMin, aabbMax;
		btTransform T;
		T.setIdentity();
		rb->getCollisionShape()->getAabb (T, aabbMin, aabbMax);
		radius = (aabbMax-aabbMin).length()*btScalar(0.5);
	}
	btScalar viscosity = btScalar(0.05);
	btVector3 force = btScalar(6.0f) * SIMD_PI * viscosity * radius * -rb->getLinearVelocity();
	
	btVector3 impulse = force * dt;
	rb->applyCentralImpulse (impulse);

	if (true)
	{
		btVector3 av = rb->getAngularVelocity();
		av *= btScalar(0.99);
		rb->setAngularVelocity (av);
	}
#else
	btScalar scaled_mu = mu * submerged_percentage * btScalar(0.4f);
	rb->applyCentralImpulse (dt * scaled_mu * -rb->getLinearVelocity());
	rb->applyTorqueImpulse (dt * scaled_mu * -rb->getAngularVelocity());
#endif
}

void btHfFluidRigidCollisionAlgorithm::processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap,
														const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	const btCollisionObjectWrapper* hfFluidWrap = (m_isSwapped) ? body1Wrap : body0Wrap;
	const btCollisionObjectWrapper* rigidWrap = (m_isSwapped) ? body0Wrap : body1Wrap;	
	
	m_hfFluid = static_cast<btHfFluid*>( const_cast<btCollisionObject*>(hfFluidWrap->getCollisionObject()) );
	m_rigidCollisionObject = const_cast<btCollisionObject*>( rigidWrap->getCollisionObject() );
	
	processGround (hfFluidWrap, rigidWrap, dispatchInfo, resultOut);
	btHfFluidBuoyantConvexShape* buoyantShape = (btHfFluidBuoyantConvexShape*)m_rigidCollisionObject->getCollisionShape();
	btRigidBody* rb = btRigidBody::upcast(m_rigidCollisionObject);
	if (rb)	//	should check if rb->getInvMass() != btScalar(0.0)
	{
		btScalar mass = btScalar(1.0f) / rb->getInvMass ();
		btScalar volume = buoyantShape->getTotalVolume ();
		btScalar density = mass / volume;
		btScalar floatyness = buoyantShape->getFloatyness ();
		btScalar submerged_volume = processFluid (dispatchInfo, density, floatyness);
		if (submerged_volume > btScalar(0.0))
		{
			btScalar submerged_percentage = submerged_volume/buoyantShape->getTotalVolume();
			//printf("%f\n", submerged_percentage);
			btScalar mu = btScalar(6.0f);
			applyFluidFriction (mu, submerged_percentage);
		}
	}
}

btScalar btHfFluidRigidCollisionAlgorithm::calculateTimeOfImpact(btCollisionObject* body0,btCollisionObject* body1,const btDispatcherInfo& dispatchInfo,btManifoldResult* resultOut)
{
	return btScalar(1.0);
}
