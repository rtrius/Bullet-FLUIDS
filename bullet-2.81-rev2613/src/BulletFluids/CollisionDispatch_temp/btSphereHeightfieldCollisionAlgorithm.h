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
#ifndef BT_SPHERE_HEIGHTFIELD_COLLISION_ALGORITHM_H
#define BT_SPHERE_HEIGHTFIELD_COLLISION_ALGORITHM_H

#include "BulletCollision/BroadphaseCollision/btDispatcher.h"
#include "BulletCollision/CollisionDispatch/btActivatingCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btCollisionCreateFunc.h"

class btPersistentManifold;

///Simplified collision algorithm; use btConvexConcaveCollisionAlgorithm for better accuracy.
///@remarks
///Rather than colliding the sphere with each triangle in its AABB, this algorithm simply
///collides the sphere with the triangle directly below its center.
///If there is no heightfield below the sphere center, the collision is ignored.
///@warning
///This algorithm requires that the btHeightfieldTerrainShape upAxis == 1, that is,
///that Y+ is up. It also assumes that the volume below the heightfield is solid.
///Spheres below the heightfield will be pushed up, even if they do not appear to be
///colliding, as long as they are within the heightfield's AABB. Since each sphere is
///collided with only 1 triangle, the radius of spheres should be smaller than the
///edge lengths of triangles.
class btSphereHeightfieldCollisionAlgorithm : public btActivatingCollisionAlgorithm
{
	bool m_ownManifold;
	btPersistentManifold* m_manifoldPtr;
	
	bool m_swapped;		///If true, implements heightfield-sphere instead of sphere-heightfield.
	
public:
	btSphereHeightfieldCollisionAlgorithm(btPersistentManifold* mf, const btCollisionAlgorithmConstructionInfo& ci,
									const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap, bool swapped);

	virtual ~btSphereHeightfieldCollisionAlgorithm()
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
			void* mem = ci.m_dispatcher1->allocateCollisionAlgorithm( sizeof(btSphereHeightfieldCollisionAlgorithm) );
			return new(mem) btSphereHeightfieldCollisionAlgorithm(ci.m_manifold, ci, col0Wrap, col1Wrap, m_swapped);
		}
	};

};

#endif //BT_SPHERE_HEIGHTFIELD_COLLISION_ALGORITHM_H
