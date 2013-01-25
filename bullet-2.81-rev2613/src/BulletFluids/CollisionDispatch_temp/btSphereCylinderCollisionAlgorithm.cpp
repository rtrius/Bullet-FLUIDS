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
#include "btSphereCylinderCollisionAlgorithm.h"

#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btCylinderShape.h"

#include "btCollide2dSphereConvex.h"

btSphereCylinderCollisionAlgorithm::btSphereCylinderCollisionAlgorithm(btPersistentManifold* mf, const btCollisionAlgorithmConstructionInfo& ci,
															const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap,
															bool swapped)
: btActivatingCollisionAlgorithm(ci, col0Wrap, col1Wrap), m_ownManifold(false), m_manifoldPtr(mf), m_swapped(swapped)
{
	if(!m_manifoldPtr)
	{
		const btCollisionObjectWrapper* sphereWrap = (m_swapped) ? col1Wrap : col0Wrap;
		const btCollisionObjectWrapper* cylinderWrap = (m_swapped) ? col0Wrap : col1Wrap;
	
		m_manifoldPtr = m_dispatcher->getNewManifold( sphereWrap->getCollisionObject(), cylinderWrap->getCollisionObject() );
		m_ownManifold = true;
	}
}

void btSphereCylinderCollisionAlgorithm::processCollision(const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap,
														const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	if(!m_manifoldPtr) return;

	resultOut->setPersistentManifold(m_manifoldPtr);
	
	const btCollisionObjectWrapper* sphereWrap = (m_swapped) ? col1Wrap : col0Wrap;
	const btCollisionObjectWrapper* cylinderWrap = (m_swapped) ? col0Wrap : col1Wrap;
	const btSphereShape* sphereShape = static_cast<const btSphereShape*>( sphereWrap->getCollisionShape() );
	const btCylinderShape* cylinderShape = static_cast<const btCylinderShape*>( cylinderWrap->getCollisionShape() );
	
	{
		const btTransform& cylinderTransform = cylinderWrap->getWorldTransform();
		
		//cylinderHalfExtents: only 2 scalars are used(1 for radius, 1 for height)
		const btVector3& cylinderHalfExtents = cylinderShape->getHalfExtentsWithoutMargin();	
		btScalar cylinderRadius = cylinderShape->getRadius();
		btScalar cylinderHalfHeight = cylinderHalfExtents.m_floats[ cylinderShape->getUpAxis() ];
	
		//Determine sphere position relative to cylinder
		btVector3 sphereRelPos = cylinderTransform.invXform( sphereWrap->getWorldTransform().getOrigin() );
		
		//Line formed by (cylinderTop, cylinderBase) is the cylinder's axis
		//btCylinderShape::getUpAxis() == 0 for X, 1 for Y, 2 for Z
		btVector3 cylinderTop(0, 0, 0);
		cylinderTop.m_floats[ cylinderShape->getUpAxis() ] = cylinderHalfHeight;	
		btVector3 cylinderBase = -cylinderTop;
		
		btScalar distance;
		btVector3 normalOnCylinder;
		btVector3 pointOnCylinder;
		
		//Check if the sphere's center is on the line formed by (cylinderTop, cylinderBase)
		btScalar sphereDistance2d = btPointLineDistance3d(cylinderTop, cylinderBase, sphereRelPos);
		if(sphereDistance2d < SIMD_EPSILON)
		{
			//The sphere is directly above or below the cylinder.
			//Handle the collision as sphere v. line segment (cylinderTop, cylinderBase).
			
			btVector3 cylinderCenter = (cylinderBase + cylinderTop) * btScalar(0.5);
			
			btVector3 cylinderCenterToSphere = sphereRelPos - cylinderCenter;
			btScalar centerToSphereDistance = cylinderCenterToSphere.length();
			
			distance = centerToSphereDistance - ( cylinderHalfHeight + sphereShape->getRadius() );
			normalOnCylinder = ( sphereRelPos.y() > btScalar(0.0) ) ? btVector3(0, 1, 0) : btVector3(0, -1, 0);
			pointOnCylinder = cylinderCenter + normalOnCylinder*cylinderHalfHeight;
		}
		else
		{
			//Since the sphere is not directly above or below, it is possible to reduce the problem into 2 dimensions
			//by forming a plane with 3 points: the cylinder top, cylinder bottom, and sphere center.
			//In this case, the sphere becomes a circle and the cylinder becomes a rectangle.
		
			//Result is a 2D plane with coordinate system (u, v), where 
			//+u is towards the sphere, -u is away from the sphere,
			//+v is towards the cylinder's top, and -v is towards the cylinder's base.
			//That is, the u coordinate is the 3D distance to the axis of the cylinder(point to line distance),
			//and the v coordinate is the point's position along the axis of the cylinder
			//(simply the X, Y, or Z coordinate since the cylinder's axis is aligned with one of those axes).
			
			//Convert the sphere's position from 3D to 2D
			btScalar spherePos_u = sphereDistance2d;		//Define the u axis so that +u is towards the sphere
			btScalar spherePos_v = sphereRelPos.m_floats[ cylinderShape->getUpAxis() ];	
			
			//Convert the cylinder from 3D to 2D; into a rectangle
			btScalar cylinderCenter_u = btScalar(0.0);
			btScalar cylinderCenter_v = btScalar(0.0);
			
			btScalar cylinderTopRight_u = cylinderRadius;
			btScalar cylinderTopRight_v = cylinderHalfHeight;
			
			btScalar cylinderTopLeft_u = -cylinderRadius;
			btScalar cylinderTopLeft_v = cylinderHalfHeight;
			
			btScalar cylinderBaseRight_u = cylinderRadius;
			btScalar cylinderBaseRight_v = -cylinderHalfHeight;
			
			btScalar cylinderBaseLeft_u = -cylinderRadius;
			btScalar cylinderBaseLeft_v = -cylinderHalfHeight;
			
			//4 edges in clockwise order
			const int NUM_EDGES = 4;
			btScalar rectanglePoints_u[NUM_EDGES + 1] = 
			{ 
				cylinderTopRight_u, 
				cylinderBaseRight_u, 
				cylinderBaseLeft_u, 
				cylinderTopLeft_u, 
				cylinderTopRight_u 
			};
			btScalar rectanglePoints_v[NUM_EDGES + 1] =
			{ 
				cylinderTopRight_v, 
				cylinderBaseRight_v, 
				cylinderBaseLeft_v, 
				cylinderTopLeft_v, 
				cylinderTopRight_v 
			};
			
			//Detect collision between a circle and box in 2D
			bool sphereInRectangle;
			btScalar distance2d, contact2d_u, contact2d_v, normal2d_u, normal2d_v;
			
			btCollide2dSphereConvex(spherePos_u, spherePos_v, sphereShape->getRadius(),
									cylinderCenter_u, cylinderCenter_v, rectanglePoints_u, rectanglePoints_v, NUM_EDGES,
									sphereInRectangle, distance2d, contact2d_u, contact2d_v, normal2d_u, normal2d_v);
			
			//Convert the results from 2D to 3D
			btVector3 u_axis = sphereRelPos / sphereDistance2d;
			u_axis.m_floats[ cylinderShape->getUpAxis() ] = btScalar(0.0);
			btVector3 v_axis(0, 0, 0);
			v_axis.m_floats[ cylinderShape->getUpAxis() ] = btScalar(1.0);
			
			distance = distance2d;
			normalOnCylinder = u_axis * normal2d_u + v_axis * normal2d_v;
			pointOnCylinder = u_axis * contact2d_u + v_axis * contact2d_v;
		}
		
		//Negative distance indicates collision/penetration
		if( distance < btScalar(0.0) ) 
		{
			btVector3 normalOnCylinderInWorld = cylinderTransform.getBasis() * normalOnCylinder.normalized();
			btVector3 pointOnCylinderInWorld = cylinderTransform(pointOnCylinder);
			
			resultOut->addContactPoint(normalOnCylinderInWorld, pointOnCylinderInWorld, distance);
		}
	}
	
	if(m_ownManifold) resultOut->refreshContactPoints();
}

btScalar btSphereCylinderCollisionAlgorithm::calculateTimeOfImpact(btCollisionObject* col0, btCollisionObject* col1,
													const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	(void)resultOut;
	(void)dispatchInfo;
	(void)col0;
	(void)col1;

	//Not implemented
	return btScalar(1.0);
}
