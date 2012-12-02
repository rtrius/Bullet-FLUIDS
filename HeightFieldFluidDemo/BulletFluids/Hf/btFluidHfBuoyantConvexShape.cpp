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
//This is an altered source version based on the HeightFieldFluidDemo included with Bullet Physics 2.80(bullet-2.80-rev2531).

#include "btFluidHfBuoyantConvexShape.h"

#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/NarrowPhaseCollision/btVoronoiSimplexSolver.h"


btFluidHfBuoyantConvexShape::btFluidHfBuoyantConvexShape (btConvexShape* convexShape)
{
	m_convexShape = convexShape;
	m_shapeType = HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE;
	m_radius = btScalar(0.f);
	m_totalVolume = btScalar(0.0f);
	m_floatyness = btScalar(1.5f);
}


//must be above the machine epsilon
#define REL_ERROR2 btScalar(1.0e-6)
static bool intersect(btVoronoiSimplexSolver* simplexSolver, 
					  const btTransform& transformA, 
					  const btTransform& transformB,
					  btConvexShape* a,
					  btConvexShape* b)
{
	
	btScalar squaredDistance = SIMD_INFINITY;
	btTransform	localTransA = transformA;
	btTransform localTransB = transformB;
	btVector3 positionOffset = (localTransA.getOrigin() + localTransB.getOrigin()) * btScalar(0.5);
	localTransA.getOrigin() -= positionOffset;
	localTransB.getOrigin() -= positionOffset;
	btScalar delta = btScalar(0.);
	btVector3 v = btVector3(1.0f, 0.0f, 0.0f);
	simplexSolver->reset ();
	do
	{
		btVector3 seperatingAxisInA = (-v)* transformA.getBasis();
		btVector3 seperatingAxisInB = v* transformB.getBasis();

		btVector3 pInA = a->localGetSupportVertexNonVirtual(seperatingAxisInA);
		btVector3 qInB = b->localGetSupportVertexNonVirtual(seperatingAxisInB);

		btVector3  pWorld = localTransA(pInA);	
		btVector3  qWorld = localTransB(qInB);

		btVector3 w	= pWorld - qWorld;
		delta = v.dot(w);

		// potential exit, they don't overlap
		if ((delta > btScalar(0.0)))
		{
			return false;
		}

		if (simplexSolver->inSimplex (w))
		{
			return false;
		}

		simplexSolver->addVertex (w, pWorld, qWorld);

		if (!simplexSolver->closest(v))
		{
			return false;
		}

		btScalar previousSquaredDistance = squaredDistance;
		squaredDistance = v.length2();

		if (previousSquaredDistance - squaredDistance <= SIMD_EPSILON * previousSquaredDistance) 
		{ 
			return false;
		}
	} 
	while (!simplexSolver->fullSimplex() && squaredDistance > REL_ERROR2 * simplexSolver->maxVertex());

    return true;
}
void btFluidHfBuoyantConvexShape::generateShape(btScalar radius, btScalar gap)
{
	btTransform identity = btTransform::getIdentity();
	btVector3 aabbMin, aabbMax;
	getAabb(identity, aabbMin, aabbMax);

	btVoronoiSimplexSolver simplexSolver;
	btSphereShape sphereShape(radius);
	
	m_voxelPositions.resize(0);
	
	btTransform voxelTransform = btTransform::getIdentity();
	
	const btScalar spacing = btScalar(2.0) * radius + gap;
	const int MAX_VOXEL_DIMENSION = 32;
	for(int i = 0; i < MAX_VOXEL_DIMENSION; i++)
	{
		for(int j = 0; j < MAX_VOXEL_DIMENSION; j++)
		{
			for(int k = 0; k < MAX_VOXEL_DIMENSION; k++)
			{
				btVector3 point( aabbMin.x() + i * spacing, aabbMin.y() + j * spacing, aabbMin.z() + k * spacing );
				if( TestPointAgainstAabb2(aabbMin, aabbMax, point) )
				{
					voxelTransform.setOrigin(point);

					if( intersect(&simplexSolver, identity, voxelTransform, m_convexShape, &sphereShape) ) m_voxelPositions.push_back(point);
				}
			}
		}
	}
	
	m_volumePerVoxel = btScalar(4.0/3.0)*SIMD_PI*radius*radius*radius;
	m_totalVolume = getNumVoxels() * m_volumePerVoxel;
	m_radius = radius;
}

