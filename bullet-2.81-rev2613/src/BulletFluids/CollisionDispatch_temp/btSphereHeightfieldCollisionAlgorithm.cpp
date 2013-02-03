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
#include "btSphereHeightfieldCollisionAlgorithm.h"

#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionDispatch/btCollisionObjectWrapper.h"

#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btHeightfieldTerrainShape.h"
#include "BulletCollision/CollisionShapes/btTriangleBuffer.h"
#include "BulletCollision/CollisionShapes/btTriangleShape.h"

#include "btCollide2dSphereConvex.h"

btSphereHeightfieldCollisionAlgorithm::btSphereHeightfieldCollisionAlgorithm(btPersistentManifold* mf, const btCollisionAlgorithmConstructionInfo& ci,
															const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap,
															bool swapped)
: btActivatingCollisionAlgorithm(ci, col0Wrap, col1Wrap), m_ownManifold(false), m_manifoldPtr(mf), m_swapped(swapped)
{
	if(!m_manifoldPtr)
	{
		const btCollisionObjectWrapper* sphereWrap = (m_swapped) ? col1Wrap : col0Wrap;
		const btCollisionObjectWrapper* heightfieldWrap = (m_swapped) ? col0Wrap : col1Wrap;
	
		m_manifoldPtr = m_dispatcher->getNewManifold( sphereWrap->getCollisionObject(), heightfieldWrap->getCollisionObject() );
		m_ownManifold = true;
	}
}

//Adds functions for simplified sphere collisions to btHeightfieldTerrainShape. Do not add any data members to this class.
class btHeightfieldTerrainShapeSimple : public btHeightfieldTerrainShape
{
public:

	int getUpAxis() const { return m_upAxis; }

	//Returns true if the point is in the heightfield's XZ AABB
	bool getTriangleInHeightfieldSpace(const btVector3& point, btVector3& out_triangle0, btVector3& out_triangle1, btVector3& out_triangle2) const
	{
		btAssert(m_upAxis == 1);	//Only Y axis heightfields are supported
		
		//Scale down the point so that it is in local (non-scaled) coordinates
		btVector3 localPoint = point * btVector3(1.f/m_localScaling[0], 1.f/m_localScaling[1], 1.f/m_localScaling[2]);
		
		//Account for local origin
		localPoint += m_localOrigin;
		
		//Quantize, with truncation
		int x = static_cast<int>( localPoint.x() );
		int z = static_cast<int>( localPoint.z() );
		
		//We want to find the lower (x,z) indicies of a quad, so the cutoff is at e.g. m_heightStickWidth-1
		bool quantizedPointOnGrid = (0 <= x && x < m_heightStickWidth-1 && 0 <= z && z < m_heightStickLength-1);
		if(!quantizedPointOnGrid) return false;
		
		//Get triangles
		btVector3 triangleA[3];
		btVector3 triangleB[3];
		
		//Shared vertices of the 2 triangles in the quad/rectangle
		btScalar seperatingLineA_x;
		btScalar seperatingLineA_z;
		btScalar seperatingLineB_x;
		btScalar seperatingLineB_z;
		
		btScalar normalTowardsA_x;		//Normal on seperating line pointing towards a vertex of triangle A
		btScalar normalTowardsA_z;
		
		btScalar discard[5];	//Unused; required by function parameters of btPointLineSegmentCollision2d()
		
		//This section should match btHeightfieldTerrainShape::processAllTriangles()
		{
			if( m_flipQuadEdges || (m_useDiamondSubdivision && !((z+x) & 1))|| (m_useZigzagSubdivision  && !(z & 1)) )
			{
				//First triangle
				getVertex(x, z, triangleA[0]);
				getVertex(x+1, z, triangleA[1]);
				getVertex(x+1, z+1, triangleA[2]);
				
				//Second triangle
				getVertex(x, z, triangleB[0]);
				getVertex(x+1, z+1, triangleB[1]);
				getVertex(x, z+1, triangleB[2]);	

				//Get the normal on the line seperating the triangles, pointing towards triangleA
				btScalar uniqueToTriangleA_x = triangleA[1].x();
				btScalar uniqueToTriangleA_z = triangleA[1].z();
				
				seperatingLineA_x = triangleA[0].x();
				seperatingLineA_z = triangleA[0].z();
				seperatingLineB_x = triangleA[2].x();
				seperatingLineB_z = triangleA[2].z();
				
				btPointLineSegmentCollision2d(seperatingLineA_x, seperatingLineA_z, seperatingLineB_x, seperatingLineB_z,
												uniqueToTriangleA_x, uniqueToTriangleA_z, 
												discard[0], discard[1], discard[2], discard[3], discard[4],
												normalTowardsA_x, normalTowardsA_z);
			}
			else
			{
				//First triangle
				getVertex(x, z, triangleA[0]);
				getVertex(x, z+1, triangleA[1]);
				getVertex(x+1, z, triangleA[2]);
				
				//Second triangle
				getVertex(x+1, z, triangleB[0]);
				getVertex(x, z+1, triangleB[1]);
				getVertex(x+1, z+1, triangleB[2]);
				
				//Get the normal on the line seperating the triangles, pointing towards triangleA
				btScalar uniqueToTriangleA_x = triangleA[0].x();
				btScalar uniqueToTriangleA_z = triangleA[0].z();
				
				seperatingLineA_x = triangleA[1].x();
				seperatingLineA_z = triangleA[1].z();
				seperatingLineB_x = triangleA[2].x();
				seperatingLineB_z = triangleA[2].z();
				
				btPointLineSegmentCollision2d(seperatingLineA_x, seperatingLineA_z, seperatingLineB_x, seperatingLineB_z,
												uniqueToTriangleA_x, uniqueToTriangleA_z, 
												discard[0], discard[1], discard[2], discard[3], discard[4],
												normalTowardsA_x, normalTowardsA_z);
			}
		}
		
		//Get the normal on the line, pointing towards the sphere
		btScalar normalTowardsPoint_x;
		btScalar normalTowardsPoint_z;
		
		btPointLineSegmentCollision2d(seperatingLineA_x, seperatingLineA_z, seperatingLineB_x, seperatingLineB_z,
										point.x(), point.z(), 
										discard[0], discard[1], discard[2], discard[3], discard[4],
										normalTowardsPoint_x, normalTowardsPoint_z);
		
		btScalar dot = normalTowardsA_x*normalTowardsPoint_x + normalTowardsA_z*normalTowardsPoint_z;
		
		//Determine which triangle the point is in
		//If the dot product is negative, the point is not on the same side of the line as triangleA
		//Since the point is in the quad, it must be in triangleB
		bool inTriangleA = (dot > 0.0);
		
		if(inTriangleA)
		{
			out_triangle0 = triangleA[0];
			out_triangle1 = triangleA[1];
			out_triangle2 = triangleA[2];
		}
		else
		{
			out_triangle0 = triangleB[0];
			out_triangle1 = triangleB[1];
			out_triangle2 = triangleB[2];
		}
		
		return true;
	}
};

void btSphereHeightfieldCollisionAlgorithm::processCollision(const btCollisionObjectWrapper* col0Wrap, const btCollisionObjectWrapper* col1Wrap,
														const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	if(!m_manifoldPtr) return;

	resultOut->setPersistentManifold(m_manifoldPtr);
	
	const btCollisionObjectWrapper* sphereWrap = (m_swapped) ? col1Wrap : col0Wrap;
	const btCollisionObjectWrapper* heightfieldWrap = (m_swapped) ? col0Wrap : col1Wrap;
	
	const btSphereShape* sphereShape = static_cast<const btSphereShape*>( sphereWrap->getCollisionShape() );
	const btHeightfieldTerrainShape* hfShape = static_cast<const btHeightfieldTerrainShape*>( heightfieldWrap->getCollisionShape() );
	
	const btHeightfieldTerrainShapeSimple* hfSimpleShape = static_cast<const btHeightfieldTerrainShapeSimple*>(hfShape);
	
	btAssert( sizeof(btHeightfieldTerrainShape) == sizeof(btHeightfieldTerrainShapeSimple) );
	btAssert( hfSimpleShape->getUpAxis() == 1 );	//Only y-axis == up heightfields are supported
	
	{
		//Determine sphere position relative to heightfield
		const btVector3& spherePosition = sphereWrap->getWorldTransform().getOrigin();
		const btTransform& heightfieldTransform = heightfieldWrap->getWorldTransform();
		btVector3 sphereRelPos = heightfieldTransform.invXform(spherePosition);
		
		btVector3 triangle[3];
		bool pointInAabbXZ = hfSimpleShape->getTriangleInHeightfieldSpace(sphereRelPos, triangle[0], triangle[1], triangle[2]);
		if(!pointInAabbXZ) return;
		
		btTriangleShape triShape(triangle[0], triangle[1], triangle[2]);
		btVector3 triangleNormal;
		triShape.calcNormal(triangleNormal);
		
		//Since the entire volume below the heightfield is considered as solid, the normal must point up(Y+)
		if( triangleNormal.y() < btScalar(0.0) ) triangleNormal = -triangleNormal;	
		
		//Point-plane distance
		btVector3 trianglePointToSphere = sphereRelPos - triangle[0];
		btScalar distance = trianglePointToSphere.dot(triangleNormal) - sphereShape->getRadius();
	
		if( distance < btScalar(0.0) )
		{
			btVector3 normalOnHeightfieldInWorld = heightfieldTransform.getBasis() * triangleNormal;
			btVector3 pointOnHeightfieldInWorld = spherePosition + normalOnHeightfieldInWorld*distance;
			
			resultOut->addContactPoint(normalOnHeightfieldInWorld, pointOnHeightfieldInWorld, distance);
		}
	}
	
	if(m_ownManifold) resultOut->refreshContactPoints();
}

btScalar btSphereHeightfieldCollisionAlgorithm::calculateTimeOfImpact(btCollisionObject* col0, btCollisionObject* col1,
													const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	(void)resultOut;
	(void)dispatchInfo;
	(void)col0;
	(void)col1;

	//Not implemented
	return btScalar(1.0);
}
