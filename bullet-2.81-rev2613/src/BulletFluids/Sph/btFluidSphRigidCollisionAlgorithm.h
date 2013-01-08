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
#ifndef BT_FLUID_RIGID_COLLISION_ALGORITHM_H
#define BT_FLUID_RIGID_COLLISION_ALGORITHM_H

#include "BulletCollision/BroadphaseCollision/btCollisionAlgorithm.h"
#include "BulletCollision/BroadphaseCollision/btDispatcher.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"

#include "btFluidSph.h"

///Implements btFluidSph - btRigidBody / btCollisionObject collision detection.
class btFluidSphRigidCollisionAlgorithm : public btCollisionAlgorithm
{
	///If true, the algorithm implements btCollisionObject - btFluidSph collision 
	///instead of btFluidSph - btCollisionObject collision.
	bool m_isSwapped;

public:
	btFluidSphRigidCollisionAlgorithm(btPersistentManifold* mf, const btCollisionAlgorithmConstructionInfo& ci,
									const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap,
									bool isSwapped) : btCollisionAlgorithm(ci), m_isSwapped(isSwapped) {}

	virtual ~btFluidSphRigidCollisionAlgorithm() {}

	virtual void processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap,
									const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
	{
		(void)dispatchInfo;
		(void)resultOut;
		
		const btCollisionObjectWrapper* fluidWrap = (m_isSwapped) ? body1Wrap : body0Wrap;
		const btCollisionObjectWrapper* rigidWrap = (m_isSwapped) ? body0Wrap : body1Wrap;
		
		const btFluidSph* constFluid = static_cast<const btFluidSph*>( fluidWrap->getCollisionObject() );
		btFluidSph* fluid = const_cast<btFluidSph*>(constFluid);
		
		//Use batched contect detection after reporting AABB intersections here
		fluid->internalGetIntersectingRigidAabbs().push_back( rigidWrap->getCollisionObject() );
	}

	virtual btScalar calculateTimeOfImpact(btCollisionObject* body0, btCollisionObject* body1, 
											const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
	{
		(void)resultOut;
		(void)dispatchInfo;
		(void)body0;
		(void)body1;

		//not yet
		return btScalar(1.);
	}

	virtual	void getAllContactManifolds(btManifoldArray& manifoldArray) {}


	struct CreateFunc : public btCollisionAlgorithmCreateFunc
	{
		virtual	btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, 
																const btCollisionObjectWrapper* body0Wrap, 
																const btCollisionObjectWrapper* body1Wrap)
		{
			void* ptr = ci.m_dispatcher1->allocateCollisionAlgorithm( sizeof(btFluidSphRigidCollisionAlgorithm) );
			return new(ptr) btFluidSphRigidCollisionAlgorithm(0, ci, body0Wrap, body1Wrap, m_swapped);
		}
	};
};

#endif