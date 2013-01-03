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
#ifndef BT_SPHERE_CAPSULE_COLLISION_ALGORITHM_H
#define BT_SPHERE_CAPSULE_COLLISION_ALGORITHM_H

#include "BulletCollision/BroadphaseCollision/btDispatcher.h"
#include "BulletCollision/CollisionDispatch/btActivatingCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"

class btPersistentManifold;

///Implements collision between btSphereShape and btCapsuleShape
class btSphereCapsuleCollisionAlgorithm : public btActivatingCollisionAlgorithm
{
	bool m_ownManifold;
	btPersistentManifold* m_manifoldPtr;
	
	bool m_swapped;		///If true, implements capsule-sphere instead of sphere-capsule.
	
public:
	btSphereCapsuleCollisionAlgorithm(btPersistentManifold* mf, const btCollisionAlgorithmConstructionInfo& ci,
									const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap, bool swapped);

	virtual ~btSphereCapsuleCollisionAlgorithm()
	{
		if(m_ownManifold && m_manifoldPtr) m_dispatcher->releaseManifold(m_manifoldPtr);
	}
	
	virtual void processCollision(const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap,
									const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);
									
	virtual btScalar calculateTimeOfImpact(btCollisionObject* col0, btCollisionObject* col1, 
											const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut);
	
	virtual	void getAllContactManifolds(btManifoldArray& manifoldArray)
	{
		if(m_manifoldPtr && m_ownManifold) manifoldArray.push_back(m_manifoldPtr);
	}
	

	struct CreateFunc : public btCollisionAlgorithmCreateFunc
	{
		virtual	btCollisionAlgorithm* CreateCollisionAlgorithm(btCollisionAlgorithmConstructionInfo& ci, 
																const btCollisionObjectWrapper* col0Wrap, 
																const btCollisionObjectWrapper* col1Wrap)
		{
			void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm( sizeof(btSphereCapsuleCollisionAlgorithm) );
			return new(mem) btSphereCapsuleCollisionAlgorithm(ci.m_manifold, ci, col0Wrap, col1Wrap, m_swapped);
		}
	};

};

#endif //BT_SPHERE_CAPSULE_COLLISION_ALGORITHM_H
