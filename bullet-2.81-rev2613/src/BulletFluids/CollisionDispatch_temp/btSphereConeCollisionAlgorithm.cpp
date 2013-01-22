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

#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btConeShape.h"

#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"

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

//Given a line segment(or line) and a point in 2D,
//return a normal on the line segment, nearest point on line segment, distance from point to line segment, and
//normal on the line.
inline void pointLineSegmentCollision2d(btScalar lineA_u, btScalar lineA_v, btScalar lineB_u, btScalar lineB_v, 
										btScalar point_u, btScalar point_v, 
										btScalar& out_normalOnSegment_u, btScalar& out_normalOnSegment_v,
										btScalar& out_pointOnSegment_u, btScalar& out_pointOnSegment_v, btScalar& out_segmentDistance,
										btScalar& out_normalOnLine_u, btScalar& out_normalOnLine_v)
{
	btScalar lineAToB_u = lineB_u - lineA_u;
	btScalar lineAToB_v = lineB_v - lineA_v;
	
	btScalar pointToLineA_u = lineA_u - point_u;
	btScalar pointToLineA_v = lineA_v - point_v;
	
	btScalar lineAToB_length2 = lineAToB_u*lineAToB_u + lineAToB_v*lineAToB_v;
	
	btScalar t_line = -(pointToLineA_u * lineAToB_u + pointToLineA_v * lineAToB_v) / lineAToB_length2;
	
	//Line segment: return normal, nearest point, and distance
	{
		btScalar t_lineSegment = btMin( btMax(btScalar(0.0), t_line), btScalar(1.0) );
		
		btScalar nearestPointOnSegment_u = lineA_u + lineAToB_u * t_lineSegment; 
		btScalar nearestPointOnSegment_v = lineA_v + lineAToB_v * t_lineSegment; 
		
		btScalar segmentToPoint_u = point_u - nearestPointOnSegment_u;
		btScalar segmentToPoint_v = point_v - nearestPointOnSegment_v;
		
		btScalar segmentToPoint_length = btSqrt(segmentToPoint_u*segmentToPoint_u + segmentToPoint_v*segmentToPoint_v);
		
		out_normalOnSegment_u = segmentToPoint_u /= segmentToPoint_length;
		out_normalOnSegment_v = segmentToPoint_v /= segmentToPoint_length;
		
		out_segmentDistance = segmentToPoint_length;
		out_pointOnSegment_u = nearestPointOnSegment_u;
		out_pointOnSegment_v = nearestPointOnSegment_v;
	}
	
	//Line: return normal on line
	{
		btScalar nearestPointOnLine_u = lineA_u + lineAToB_u * t_line; 
		btScalar nearestPointOnLine_v = lineA_v + lineAToB_v * t_line; 
		
		btScalar lineToPoint_u = point_u - nearestPointOnLine_u;
		btScalar lineToPoint_v = point_v - nearestPointOnLine_v;
		
		btScalar lineToPoint_length = btSqrt(lineToPoint_u*lineToPoint_u + lineToPoint_v*lineToPoint_v);
		
		out_normalOnLine_u = lineToPoint_u / lineToPoint_length;
		out_normalOnLine_v = lineToPoint_v / lineToPoint_length;
	}
}

inline btScalar pointLineDistance3d(const btVector3& lineA, const btVector3& lineB, const btVector3& point)
{
	//Line equation: lineA + lineAToB * t
	btVector3 lineAToB = (lineB - lineA).normalized();
	
	btVector3 pointToLineA = lineA - point;

	btVector3 pointToLineNearest = pointToLineA - pointToLineA.dot(lineAToB)*lineAToB;
	
	return pointToLineNearest.length();
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
	
	//Initially support only Y axis cones
	btAssert(coneShape->getConeUpIndex() == 1);
	
	{
		const btTransform& coneTransform = coneWrap->getWorldTransform();
		
		//Determine sphere position relative to cone
		btVector3 sphereRelPos = coneTransform.invXform( sphereWrap->getWorldTransform().getOrigin() );
		
		//(Assuming a cone that points upwards in the Y direction - the 'circle/disk' end is at Y-, while the pointy end is at Y+),
		btVector3 coneTop(0, 0, 0);
		coneTop.m_floats[ coneShape->getConeUpIndex() ] = coneShape->getHeight()*btScalar(0.5);	//getConeUpIndex() == 0 for X, 1 for Y, 2 for Z
		btVector3 coneBase = -coneTop;
		
		btScalar distance;
		btVector3 normalOnCone;
		btVector3 pointOnCone;
		
		//Check if the sphere's center is on the line formed by (coneTop, coneBase)
		btScalar sphereDistance2d = pointLineDistance3d(coneTop, coneBase, sphereRelPos);
		if(sphereDistance2d < SIMD_EPSILON)
		{
			//The sphere is directly above or below the cone.
			//Handle the collision as sphere v. line segment (coneTop, coneBase).
			
			btVector3 coneCenter = (coneBase + coneTop) * btScalar(0.5);
			
			btVector3 coneCenterToSphere = sphereRelPos - coneCenter;
			btScalar centerToSphereDistance = coneCenterToSphere.length();
			
			distance = centerToSphereDistance - ( coneShape->getHeight()*btScalar(0.5) + sphereShape->getRadius() );
			normalOnCone = ( sphereRelPos.y() > btScalar(0.0) ) ? btVector3(0, 1, 0) : btVector3(0, -1, 0);
			pointOnCone = coneCenter + normalOnCone*( distance * btScalar(0.5) );
		}
		else
		{
			//Since the sphere is not directly above or below, it is possible to reduce the problem into 2 dimensions
			//by forming a plane with 3 points: the cone top, cone bottom, and sphere center.
			//In this case, the sphere becomes a circle and the cone becomes a triangle.
		
			//Create a normalized 2D vector from the cone center to the sphere center
			//Since the cone center is at (0, y, 0), the vector is simply the normalized 2d position of the sphere
			btScalar normal_x = sphereRelPos.x() / sphereDistance2d;
			btScalar normal_z = sphereRelPos.z() / sphereDistance2d;
			
			//Result is a 2D plane with coordinate system (u, v)
			//where the conversion between 2D(u, v) and 3D(x, y, z) is that 
			//u == btVector3(normal_x, 0, normal_z) (that is, u == the distance in the xz plane), and 
			//v == btVector3(0, v, 0);
			
			//Convert from 2D to 3D: btVector3(normal_x*u, v, normal_z*u);		//input(u, v)
			//Convert from 3D to 2D: btVector2(distance_from_origin_xz, y);		//input(x, y, z)
			
			//Use that vector to create a triangle, composed of the top point and 2 base points
			btScalar coneRadius = coneShape->getRadius();
			
			btScalar spherePos_u = sphereDistance2d;	//sign is included in normal_x/z -- define u coordinate so that +u is towards the sphere
			btScalar spherePos_v = sphereRelPos.y(); 
			
			btScalar coneTop_u = btScalar(0.0);
			btScalar coneTop_v = coneShape->getHeight() * btScalar(0.5);
			btScalar coneCenter_u = btScalar(0.0);
			btScalar coneCenter_v = btScalar(0.0);
			
			btScalar coneBaseA_u = coneRadius;
			btScalar coneBaseA_v = -coneShape->getHeight() * btScalar(0.5);
			
			btScalar coneBaseB_u = -coneRadius;
			btScalar coneBaseB_v = -coneShape->getHeight() * btScalar(0.5);
			
			//Check if the sphere is inside the triangle by comparing the position of the sphere
			//to the triangle's edges(which side of line)
			btScalar trianglePoints_u[4] = { coneTop_u, coneBaseA_u, coneBaseB_u, coneTop_u };	//3 edges
			btScalar trianglePoints_v[4] = { coneTop_v, coneBaseA_v, coneBaseB_v, coneTop_v };
			
			
			btScalar nearestLineSegmentDistance(BT_LARGE_FLOAT);
			btScalar nearestLineSegmentContact_u, nearestLineSegmentContact_v;
			btScalar nearestLineSegmentNormal_u, nearestLineSegmentNormal_v;
			
			bool sphereInTriangle = true;
			for(int i = 0; i < 3; ++i)
			{
				btScalar lineA_u = trianglePoints_u[i];
				btScalar lineA_v = trianglePoints_v[i];
				
				btScalar lineB_u = trianglePoints_u[i+1];
				btScalar lineB_v = trianglePoints_v[i+1];
				
				//Calculate normal by projecting the cone center onto the line
				btScalar discard[5];	//Unused
				btScalar triangleNormal_u, triangleNormal_v;
				pointLineSegmentCollision2d(lineA_u, lineA_v, lineB_u, lineB_v, coneCenter_u, coneCenter_v, 
											discard[0], discard[1], discard[2], discard[3], discard[4],
											triangleNormal_u, triangleNormal_v);
				triangleNormal_u = -triangleNormal_u;	//Since the cone center is used, the normal will point inwards, towards the triangle
				triangleNormal_v = -triangleNormal_v;	//Invert it to make it point outwards
				
				//Use normal to test which side of line sphere is on
				btScalar normalOnSegment_u, normalOnSegment_v, pointOnSegment_u, pointOnSegment_v, distanceToSegment; 
				btScalar normalOnLine_u, normalOnLine_v; 
				pointLineSegmentCollision2d(lineA_u, lineA_v, lineB_u, lineB_v, spherePos_u, spherePos_v, 
											normalOnSegment_u, normalOnSegment_v, pointOnSegment_u, pointOnSegment_v, distanceToSegment,
											normalOnLine_u, normalOnLine_v);
				
				btScalar dotProduct = triangleNormal_u*normalOnLine_u + triangleNormal_v*normalOnLine_v;
				
				if(dotProduct > 0.0) 	//Sphere is not on the same side of the line as the cone center
				{
					sphereInTriangle = false;
				}
				
				if(distanceToSegment < nearestLineSegmentDistance)
				{
					nearestLineSegmentDistance = distanceToSegment;
					nearestLineSegmentContact_u = pointOnSegment_u;
					nearestLineSegmentContact_v = pointOnSegment_v;
					nearestLineSegmentNormal_u = normalOnSegment_u;
					nearestLineSegmentNormal_v = normalOnSegment_v;
				}
			}
			
			btScalar distance2d;
			btScalar contact2d_u, contact2d_v;
			btScalar normal2d_u, normal2d_v;
			
			contact2d_u = nearestLineSegmentContact_u;
			contact2d_v = nearestLineSegmentContact_v;
			
			if(!sphereInTriangle)
			{
				distance2d = nearestLineSegmentDistance - sphereShape->getRadius();
				
				//Since the sphere is outside the triangle, the line normal already points outwards
				normal2d_u = nearestLineSegmentNormal_u;
				normal2d_v = nearestLineSegmentNormal_v;
			}
			else
			{
				distance2d = -nearestLineSegmentDistance - sphereShape->getRadius();
			
				//The sphere is inside the triangle, so the line normal will point into the triangle
				//Invert the line segment normal to make it point away from the triangle center
				normal2d_u = -nearestLineSegmentNormal_u;
				normal2d_v = -nearestLineSegmentNormal_v;
			}
			
			//Convert the results from 2D to 3D
			btVector3 u_axis(normal_x, 0.0, normal_z);	//	sphereRelPos * u_axis_pt / distance2d;
			btVector3 v_axis(0.0, 1.0, 0.0);
			
			normalOnCone = u_axis * normal2d_u + v_axis * normal2d_v;
			
			distance = distance2d;
			normalOnCone.setValue(normal2d_u*normal_x, normal2d_v, normal2d_u*normal_z);
			pointOnCone.setValue(contact2d_u*normal_x, contact2d_v, contact2d_u*normal_z);
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
