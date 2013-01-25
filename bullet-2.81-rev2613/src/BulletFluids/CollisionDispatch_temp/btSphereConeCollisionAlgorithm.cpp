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
#include "btSphereConeCollisionAlgorithm.h"

#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btConeShape.h"

#include "btCollide2dSphereConvex.h"

btSphereConeCollisionAlgorithm::btSphereConeCollisionAlgorithm(btPersistentManifold* mf, const btCollisionAlgorithmConstructionInfo& ci,
															const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap,
															bool swapped)
: btActivatingCollisionAlgorithm(ci, col0Wrap, col1Wrap), m_ownManifold(false), m_manifoldPtr(mf), m_swapped(swapped)
{
	if(!m_manifoldPtr)
	{
		const btCollisionObjectWrapper* sphereWrap = (m_swapped) ? col1Wrap : col0Wrap;
		const btCollisionObjectWrapper* coneWrap = (m_swapped) ? col0Wrap : col1Wrap;
	
		m_manifoldPtr = m_dispatcher->getNewManifold( sphereWrap->getCollisionObject(), coneWrap->getCollisionObject() );
		m_ownManifold = true;
	}
}

void btSphereConeCollisionAlgorithm::processCollision(const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap,
														const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	if(!m_manifoldPtr) return;

	resultOut->setPersistentManifold(m_manifoldPtr);
	
	const btCollisionObjectWrapper* sphereWrap = (m_swapped) ? col1Wrap : col0Wrap;
	const btCollisionObjectWrapper* coneWrap = (m_swapped) ? col0Wrap : col1Wrap;
	const btSphereShape* sphereShape = static_cast<const btSphereShape*>( sphereWrap->getCollisionShape() );
	const btConeShape* coneShape = static_cast<const btConeShape*>( coneWrap->getCollisionShape() );
	
	{
		const btTransform& coneTransform = coneWrap->getWorldTransform();
		
		btScalar coneRadius = coneShape->getRadius();
		btScalar coneHalfHeight = coneShape->getHeight() * btScalar(0.5);
			
		//Determine sphere position relative to cone
		btVector3 sphereRelPos = coneTransform.invXform( sphereWrap->getWorldTransform().getOrigin() );
		
		//Line formed by (coneTop, coneBase) is the cone's axis
		btVector3 coneTop(0, 0, 0);
		coneTop.m_floats[ coneShape->getConeUpIndex() ] = coneHalfHeight;	//getConeUpIndex() == 0 for X, 1 for Y, 2 for Z
		btVector3 coneBase = -coneTop;
		
		btScalar distance;
		btVector3 normalOnCone;
		btVector3 pointOnCone;
		
		//Check if the sphere's center is on the line formed by (coneTop, coneBase)
		btScalar sphereDistance2d = btPointLineDistance3d(coneTop, coneBase, sphereRelPos);
		if(sphereDistance2d < SIMD_EPSILON)
		{
			//The sphere is directly above or below the cone.
			//Handle the collision as sphere v. line segment (coneTop, coneBase).
			
			btVector3 coneCenter = (coneBase + coneTop) * btScalar(0.5);
			
			btVector3 coneCenterToSphere = sphereRelPos - coneCenter;
			btScalar centerToSphereDistance = coneCenterToSphere.length();
			
			distance = centerToSphereDistance - ( coneHalfHeight + sphereShape->getRadius() );
			normalOnCone = ( sphereRelPos.y() > btScalar(0.0) ) ? btVector3(0, 1, 0) : btVector3(0, -1, 0);
			pointOnCone = coneCenter + normalOnCone*coneHalfHeight;
		}
		else
		{
			//Since the sphere is not directly above or below, it is possible to reduce the problem into 2 dimensions
			//by forming a plane with 3 points: the cone top, cone bottom, and sphere center.
			//In this case, the sphere becomes a circle and the cone becomes a triangle.
		
			//Result is a 2D plane with coordinate system (u, v), where 
			//+u is towards the sphere, -u is away from the sphere,
			//+v is towards the cone's pointy end, and -v is towards the cone's base.
			//That is, the u coordinate is the 3D distance to the axis of the cone(point to line distance),
			//and the v coordinate is the point's position along the axis of the cone
			//(simply the X, Y, or Z coordinate since the cone's axis is aligned with one of those axes).
			
			//Convert the sphere's position from 3D to 2D
			btScalar spherePos_u = sphereDistance2d;		//Define the u axis so that +u is towards the sphere
			btScalar spherePos_v = sphereRelPos.m_floats[ coneShape->getConeUpIndex() ];	
			
			//Convert the cone from 3D to 2D; into a triangle composed of the top point and 2 base points
			btScalar coneCenter_u = btScalar(0.0);
			btScalar coneCenter_v = btScalar(0.0);
			
			btScalar coneTop_u = btScalar(0.0);
			btScalar coneTop_v = coneHalfHeight;
			
			btScalar coneBaseA_u = coneRadius;
			btScalar coneBaseA_v = -coneHalfHeight;
			
			btScalar coneBaseB_u = -coneRadius;
			btScalar coneBaseB_v = -coneHalfHeight;
			
			//
			const int NUM_EDGES = 3;
			btScalar trianglePoints_u[NUM_EDGES + 1] = { coneTop_u, coneBaseA_u, coneBaseB_u, coneTop_u };
			btScalar trianglePoints_v[NUM_EDGES + 1] = { coneTop_v, coneBaseA_v, coneBaseB_v, coneTop_v };
			
			//Detect collision between a circle and triangle in 2D
			bool sphereInTriangle;
			btScalar distance2d, contact2d_u, contact2d_v, normal2d_u, normal2d_v;
			
			btCollide2dSphereConvex(spherePos_u, spherePos_v, sphereShape->getRadius(),
									coneCenter_u, coneCenter_v, trianglePoints_u, trianglePoints_v, NUM_EDGES,
									sphereInTriangle, distance2d, contact2d_u, contact2d_v, normal2d_u, normal2d_v);
			
			//Convert the results from 2D to 3D
			btVector3 u_axis = sphereRelPos / sphereDistance2d;
			u_axis.m_floats[ coneShape->getConeUpIndex() ] = btScalar(0.0);
			btVector3 v_axis(0, 0, 0);
			v_axis.m_floats[ coneShape->getConeUpIndex() ] = btScalar(1.0);
			
			distance = distance2d;
			normalOnCone = u_axis * normal2d_u + v_axis * normal2d_v;
			pointOnCone = u_axis * contact2d_u + v_axis * contact2d_v;
		}
		
		//Negative distance indicates collision/penetration
		if( distance < btScalar(0.0) ) 
		{
			btVector3 normalOnConeInWorld = coneTransform.getBasis() * normalOnCone.normalized();
			btVector3 pointOnConeInWorld = coneTransform(pointOnCone);
			
			resultOut->addContactPoint(normalOnConeInWorld, pointOnConeInWorld, distance);
		}
	}
	
	if(m_ownManifold) resultOut->refreshContactPoints();
}

btScalar btSphereConeCollisionAlgorithm::calculateTimeOfImpact(btCollisionObject* col0, btCollisionObject* col1,
													const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	(void)resultOut;
	(void)dispatchInfo;
	(void)col0;
	(void)col1;

	//Not implemented
	return btScalar(1.0);
}
