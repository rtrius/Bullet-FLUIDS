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
	m_collisionRadius = btScalar(0.0);
	m_totalVolume = btScalar(0.0);
	
	m_useFluidDensity = false;
	m_buoyancyScale = btScalar(1.0);
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
void btFluidHfBuoyantConvexShape::generateShape(btScalar collisionRadius, btScalar volumeEstimationRadius)
{
	m_voxelPositions.resize(0);
	
	btVoronoiSimplexSolver simplexSolver;
	btTransform voxelTransform = btTransform::getIdentity();
	
	btVector3 aabbMin, aabbMax;
	getAabb( btTransform::getIdentity(), aabbMin, aabbMax );	//AABB of the convex shape
	
	//Generate voxels for collision
	{
		btSphereShape collisionShape(collisionRadius);
	
		const btScalar collisionDiameter = btScalar(2.0) * collisionRadius;
		btVector3 numVoxels = (aabbMax - aabbMin) / collisionDiameter;
		for(int i = 0; i < ceil( numVoxels.x() ); i++)
		{
			for(int j = 0; j < ceil( numVoxels.y() ); j++)
			{
				for(int k = 0; k < ceil( numVoxels.z() ); k++)
				{
					btVector3 voxelPosition = aabbMin + btVector3(i * collisionDiameter, j * collisionDiameter, k * collisionDiameter);
					voxelTransform.setOrigin(voxelPosition);

					if( intersect(&simplexSolver, btTransform::getIdentity(), voxelTransform, m_convexShape, &collisionShape) ) 
						m_voxelPositions.push_back(voxelPosition);
				}
			}
		}
		
		const bool CENTER_VOXELS_AABB_ON_ORIGIN = true;
		if(CENTER_VOXELS_AABB_ON_ORIGIN)
		{
			btVector3 diameter(collisionDiameter, collisionDiameter, collisionDiameter);
		
			btVector3 voxelAabbMin(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT); 
			btVector3 voxelAabbMax(-BT_LARGE_FLOAT, -BT_LARGE_FLOAT, -BT_LARGE_FLOAT);
			for(int i = 0; i < m_voxelPositions.size(); ++i)
			{
				voxelAabbMin.setMin(m_voxelPositions[i] - diameter);
				voxelAabbMax.setMax(m_voxelPositions[i] + diameter);
			}
		
			btVector3 offset = (voxelAabbMax - voxelAabbMin)*btScalar(0.5) - voxelAabbMax;
		
			for(int i = 0; i < m_voxelPositions.size(); ++i) m_voxelPositions[i] += offset;
		}
	}
	
	//Estimate volume with smaller spheres
	btScalar estimatedVolume;
	{
		btSphereShape volumeEstimationShape(volumeEstimationRadius);	
		
		int numCollidingVoxels = 0;
		const btScalar estimationDiameter = btScalar(2.0) * volumeEstimationRadius;
		btVector3 numEstimationVoxels = (aabbMax - aabbMin) / estimationDiameter;
		
		for(int i = 0; i < ceil( numEstimationVoxels.x() ); i++)
		{
			for(int j = 0; j < ceil( numEstimationVoxels.y() ); j++)
			{
				for(int k = 0; k < ceil( numEstimationVoxels.z() ); k++)
				{
					btVector3 voxelPosition = aabbMin + btVector3(i * estimationDiameter, j * estimationDiameter, k * estimationDiameter);
					voxelTransform.setOrigin(voxelPosition);

					if( intersect(&simplexSolver, btTransform::getIdentity(), voxelTransform, m_convexShape, &volumeEstimationShape) ) 
						++numCollidingVoxels;
				}
			}
		}
		

		
		//Although the voxels are spherical, it is better to use the volume of a cube
		//for volume estimation. Since convex shapes are completely solid and the voxels
		//are generated by moving along a cubic lattice, using the volume of a sphere
		//would result in gaps. Comparing the volume of a cube with edge length 2(8 m^3) 
		//and the volume of a sphere with diameter 2(~4 m^3), the estimated volume would
		//be off by about 1/2.
		btScalar volumePerEstimationVoxel = btPow( estimationDiameter, btScalar(3.0) );
		//btScalar volumePerEstimationVoxel =  btScalar(4.0/3.0) * SIMD_PI * btPow( volumeEstimationRadius, btScalar(3.0) );
		
		estimatedVolume = static_cast<btScalar>(numCollidingVoxels) * volumePerEstimationVoxel;
	}
	
	m_volumePerVoxel = estimatedVolume / static_cast<btScalar>( getNumVoxels() );
	m_totalVolume = estimatedVolume;
	
	m_collisionRadius = collisionRadius;
}

