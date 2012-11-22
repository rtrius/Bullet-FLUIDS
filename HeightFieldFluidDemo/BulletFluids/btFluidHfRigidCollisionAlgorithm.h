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

#ifndef BT_FLUID_HF_RIGID_COLLISION_ALGORITHM_H
#define BT_FLUID_HF_RIGID_COLLISION_ALGORITHM_H

#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"
#include "BulletCollision/NarrowPhaseCollision/btPersistentManifold.h"
#include "BulletCollision/CollisionShapes/btTriangleCallback.h"
#include "BulletCollision/CollisionDispatch/btConvexConcaveCollisionAlgorithm.h"
#include "LinearMath/btVector3.h"

class btFluidHf;

/// btFluidHfRigidCollisionAlgorithm provides collision detection between btFluidHf and btRigidBody
class btFluidHfRigidCollisionAlgorithm : public btCollisionAlgorithm
{
	btFluidHf*				m_hfFluid;
	btCollisionObject*		m_rigidCollisionObject;

	///for rigid versus fluid (instead of fluid versus rigid), we use this swapped boolean
	bool	m_isSwapped;

	btConvexTriangleCallback m_convexTrianglecallback;

	void processGround (const btCollisionObjectWrapper* hfFluidWrap, const btCollisionObjectWrapper* rigidWrap,
						const btDispatcherInfo& dispatchInfo,btManifoldResult* resultOut);
	void applyFluidFriction (btScalar mu, btScalar submerged_percentage);
	btScalar processFluid (const btDispatcherInfo& dispatchInfo, btScalar density, btScalar floatyness);
	
public:
	btFluidHfRigidCollisionAlgorithm(const btCollisionAlgorithmConstructionInfo& ci,
									const btCollisionObjectWrapper* body0Wrap, 
									const btCollisionObjectWrapper* body1Wrap, bool isSwapped);

	virtual ~btFluidHfRigidCollisionAlgorithm() {}
	
	virtual void processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap,
								  const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);
								  
	virtual btScalar calculateTimeOfImpact(btCollisionObject* body0,btCollisionObject* body1,const btDispatcherInfo& dispatchInfo,btManifoldResult* resultOut);

	virtual	void getAllContactManifolds(btManifoldArray& manifoldArray) 
	{
		if(m_convexTrianglecallback.m_manifoldPtr) manifoldArray.push_back(m_convexTrianglecallback.m_manifoldPtr); 
	}


	struct CreateFunc : public btCollisionAlgorithmCreateFunc
	{
		virtual	btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, 
																const btCollisionObjectWrapper* body0Wrap, 
																const btCollisionObjectWrapper* body1Wrap)
		{
			void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm(sizeof(btFluidHfRigidCollisionAlgorithm));
			return new(mem) btFluidHfRigidCollisionAlgorithm(ci, body0Wrap, body1Wrap, m_swapped);
		}
	};
};

#endif //BT_FLUID_HF_RIGID_COLLISION_ALGORITHM_H
