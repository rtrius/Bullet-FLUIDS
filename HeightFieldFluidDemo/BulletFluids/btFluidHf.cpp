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
#include "btFluidHf.h"

#include "BulletDynamics/Dynamics/btRigidBody.h"
#include "LinearMath/btQuickProf.h"

#include "btFluidHfCollisionShape.h"
#include "btFluidHfBuoyantConvexShape.h"


btFluidHf::btFluidHf (btScalar gridCellWidth, int numNodesX, int numNodesZ)
{
	m_columns.setGridDimensions(gridCellWidth, numNodesX, numNodesZ);
	
	m_aabbMin = btVector3(0.0, 0.0, 0.0);
	m_aabbMax = btVector3( m_columns.m_numNodesX * m_columns.m_gridCellWidth, 
							0.0, 
							m_columns.m_numNodesZ * m_columns.m_gridCellWidth );
	
	//btCollisionObject
	{
		m_internalType = CO_HF_FLUID;
		m_collisionShape = new btFluidHfCollisionShape(this);
	}
}
btFluidHf::~btFluidHf ()
{
	delete m_collisionShape;
}

void btFluidHf::stepSimulation(btScalar dt)
{
	BT_PROFILE("btFluidHf::stepSimulation()");

	m_solver.stepSimulation(dt, m_hfParameters, m_columns);
	
	//Update AABB; x and z is constant
	{
		btScalar minY = m_columns.m_ground[0];
		for(int i = 1; i < m_columns.m_ground.size(); ++i) minY = btMin(minY, m_columns.m_ground[i]);
		m_aabbMin.setY(minY);
		
		btScalar maxY = m_columns.m_combinedHeight[0];
		for(int i = 1; i < m_columns.m_combinedHeight.size(); ++i) maxY = btMax(maxY, m_columns.m_combinedHeight[i]);
		m_aabbMax.setY(maxY);
	}
}

void btFluidHf::prep ()
{
	for(int i = 0; i < m_columns.m_numNodesZ*m_columns.m_numNodesX; i++) m_columns.m_combinedHeight[i] = m_columns.m_fluidDepth[i] + m_columns.m_ground[i];
	btFluidHfSolverDefault::computeFlagsAndFillRatio(m_hfParameters, m_columns);
}

void btFluidHf::setFluidHeight (int index, btScalar height)
{
	m_columns.m_fluidDepth[index] = height;
	m_columns.m_combinedHeight[index] = m_columns.m_ground[index] + m_columns.m_fluidDepth[index];
	m_columns.m_active[index] = true;
}

void btFluidHf::addFluidHeight (int x, int y, btScalar height)
{
	int index = arrayIndex (x,y);
	m_columns.m_fluidDepth[index] += height;
	m_columns.m_combinedHeight[index] = m_columns.m_ground[index] + m_columns.m_fluidDepth[index];
	m_columns.m_active[index] = true;
}

void btFluidHf::getAabbForColumn (int i, int j, btVector3& aabbMin, btVector3& aabbMax)
{
	const btVector3& origin = getWorldTransform().getOrigin();
	int index = arrayIndex(i, j);

	aabbMin = btVector3( getCellPosX(i), m_columns.m_ground[index], getCellPosZ(j) ) + origin;
	aabbMax = btVector3( getCellPosX(i+1), m_columns.m_combinedHeight[index], getCellPosZ(j+1) ) + origin;
}

void btFluidHf::forEachFluidColumn (btFluidHfColumnCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax)
{
	int startNodeX, startNodeZ, endNodeX, endNodeZ;

	startNodeX = (int)( aabbMin.x() * m_columns.m_gridCellWidthInv );
	startNodeZ = (int)( aabbMin.z() * m_columns.m_gridCellWidthInv );

	endNodeX = (int)( aabbMax.x() * m_columns.m_gridCellWidthInv ) + 1;
	endNodeZ = (int)( aabbMax.z() * m_columns.m_gridCellWidthInv ) + 1;
	
	startNodeX = btMax (1, startNodeX);
	startNodeZ = btMax (1, startNodeZ);
	endNodeX = btMin (m_columns.m_numNodesX-2, endNodeX);
	endNodeZ = btMin (m_columns.m_numNodesZ-2, endNodeZ);

	for (int i = startNodeX; i < endNodeX; i++)
	{
		for (int j = startNodeZ; j < endNodeZ; j++)
		{
			if ( !m_columns.m_active[arrayIndex(i, j)] ) continue;

			if ( !callback->processColumn(this, i, j) ) return;
		}
	}
}

void btFluidHf::forEachTriangle(btTriangleCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax, 
								const btAlignedObjectArray<btScalar>& heightArray) const
{
	//If aabbMax is not clamped, it may produce a negative endNodeX if aabbMax.x() == BT_LARGE_FLOAT
	//static_cast<unsigned int>(BT_LARGE_FLOAT * m_columns.m_gridCellWidthInv) == 0
	//static_cast<int>(BT_LARGE_FLOAT * m_columns.m_gridCellWidthInv) == lowest(negative) int
	btVector3 clampedMin = aabbMin;
	btVector3 clampedMax = aabbMax;
	clampedMin.setMax(m_aabbMin);		//Not 'reversed'; clamp aabbMin/aabbMax to the heightfield's m_aabbMin/m_aabbMax
	clampedMax.setMin(m_aabbMax);

	int startNodeX, startNodeZ, endNodeX, endNodeZ;
	
	startNodeX = (int)( clampedMin.x() * m_columns.m_gridCellWidthInv );
	startNodeZ = (int)( clampedMin.z() * m_columns.m_gridCellWidthInv );

	endNodeX = (int)( clampedMax.x() / m_columns.m_gridCellWidth) + 1;
	endNodeZ = (int)( clampedMax.z() * m_columns.m_gridCellWidthInv ) + 1;
	
	startNodeX = btMax (0, startNodeX);
	startNodeZ = btMax (0, startNodeZ);
	endNodeX = btMin (m_columns.m_numNodesX-1, endNodeX);
	endNodeZ = btMin (m_columns.m_numNodesZ-1, endNodeZ);
	
	for (int j = startNodeZ; j < endNodeZ; j++)
	{
		for (int i = startNodeX; i < endNodeX; i++)
		{
			btVector3 sw = btVector3( getCellPosX(i), heightArray[arrayIndex(i, j)], getCellPosZ(j) );
			btVector3 se = btVector3( getCellPosX(i+1), heightArray[arrayIndex(i+1, j)], getCellPosZ(j) );
			btVector3 nw = btVector3( getCellPosX(i), heightArray[arrayIndex(i, j+1)], getCellPosZ(j+1) );
			btVector3 ne = btVector3( getCellPosX(i+1), heightArray[arrayIndex(i+1, j+1)], getCellPosZ(j+1) );
		
			btVector3 verts[3];
			// triangle 1
			verts[0] = sw;
			verts[1] = nw;
			verts[2] = se;
			callback->processTriangle(verts,i,j);
			
			// triangle 2
			verts[0] = se;
			verts[1] = nw;
			verts[2] = ne;
			callback->processTriangle(verts,i,j);
		}
	}
}

static btScalar rangeOverlap (btScalar lo1, btScalar hi1, btScalar lo2, btScalar hi2, btScalar& loOut, btScalar& hiOut)
{
	if (!(lo1 <= hi2 && lo2 <= hi1))
		return btScalar(0.0f);

	if (lo1 >= lo2 && lo1 <= hi2 &&
		hi1 >= lo2 && hi1 <= hi2)
	{
		hiOut = hi1;
		loOut = lo1;
		return hi1 - lo1;
	} else if (lo2 >= lo1 && lo2 <= hi1 &&
			   hi2 >= lo1 && hi2 <= hi1)
	{
		hiOut = hi2;
		loOut = lo2;
		return hi2 - lo2;
	} else if (hi1 >= lo2 && lo1 <= lo2) {
		hiOut = hi1;
		loOut = lo2;
		return hi1 - lo2;
	} else {
		hiOut = hi2;
		loOut = lo1;
		return hi2 - lo1;
	}
}

btFluidHfColumnRigidBodyCallback::btFluidHfColumnRigidBodyCallback (btRigidBody* rigidBody, btIDebugDraw* debugDraw, btScalar density, btScalar floatyness)
{
	m_rigidBody = rigidBody;
	m_buoyantShape = (btFluidHfBuoyantConvexShape*)rigidBody->getCollisionShape();
	m_debugDraw = debugDraw;
	m_rigidBody->getAabb (m_aabbMin, m_aabbMax);
	m_volume = btScalar(0.0f);
	m_density = density;
	m_floatyness = floatyness;
	
	m_numVoxels = m_buoyantShape->getNumVoxels();
	m_voxelPositionsXformed = (btVector3*)btAlignedAlloc(sizeof(btVector3)*m_numVoxels, 16);
	m_voxelSubmerged = (bool*)btAlignedAlloc(sizeof(bool)*m_numVoxels, 16);
	
	for (int i = 0; i < m_numVoxels; i++)
	{
		btVector3 p = m_buoyantShape->getVoxelPositionsArray()[i];
		p = m_rigidBody->getWorldTransform().getBasis() * p;
		p += m_rigidBody->getWorldTransform().getOrigin();
		m_voxelPositionsXformed[i] = p;
		m_voxelSubmerged[i] = false;
	}
}

btFluidHfColumnRigidBodyCallback::~btFluidHfColumnRigidBodyCallback ()
{
	if (m_voxelPositionsXformed)
		btAlignedFree (m_voxelPositionsXformed);
	if (m_voxelSubmerged)
		btAlignedFree (m_voxelSubmerged);
}

static bool sphereVsAABB (const btVector3& aabbMin, const btVector3& aabbMax, const btVector3& sphereCenter, const btScalar sphereRadius)
{
  btScalar totalDistance = 0;

  // Accumulate the distance of the sphere's center on each axis
  for(int i = 0; i < 3; ++i) {

    // If the sphere's center is outside the aabb, we've got distance on this axis
    if(sphereCenter[i] < aabbMin[i]) {
      btScalar borderDistance = aabbMin[i] - sphereCenter[i];
      totalDistance += borderDistance * borderDistance;

    } else if(sphereCenter[i] > aabbMax[i]) {
      btScalar borderDistance = sphereCenter[i] - aabbMax[i];
      totalDistance += borderDistance * borderDistance;

    }
    // Otherwise the sphere's center is within the box on this axis, so the
    // distance will be 0 and we do not need to accumulate anything at all

  }

  // If the distance to the box is lower than the sphere's radius, both are overlapping
  return (totalDistance <= (sphereRadius * sphereRadius));
}

bool btFluidHfColumnRigidBodyCallback::processColumn (btFluidHf* fluid, int w, int l)
{
	btVector3 columnAabbMin,columnAabbMax;

	fluid->getAabbForColumn (w, l, columnAabbMin, columnAabbMax);

	bool applyBuoyancyImpulse = true;
	bool applyFluidVelocityImpulse = true;
	bool applyFluidDisplace = true;

	btScalar dt = btScalar(1.0f/60.0f);

	btScalar columnVolume = btScalar(0.0f);

	for (int i = 0; i < m_buoyantShape->getNumVoxels(); i++)
	{
		if (m_voxelSubmerged[i])
			continue;

		if (sphereVsAABB(columnAabbMin, columnAabbMax, m_voxelPositionsXformed[i], m_buoyantShape->getVoxelRadius()))
		{
			m_voxelSubmerged[i] = true;
			btScalar voxelVolume = m_buoyantShape->getVolumePerVoxel();
			columnVolume += voxelVolume;

			btVector3 application_point = m_voxelPositionsXformed[i];
			btVector3 relative_position = application_point - m_rigidBody->getCenterOfMassPosition();

			if (applyBuoyancyImpulse)
			{
				btScalar massDisplacedWater = voxelVolume * m_density * m_floatyness;
				btScalar force = massDisplacedWater * -fluid->getGravity();
				btScalar impulseMag = force * dt;
				btVector3 impulseVec = btVector3(btScalar(0.0f), btScalar(1.0f), btScalar(0.0f)) * impulseMag;
//#define CENTER_IMPULSE_ONLY 1
#ifdef CENTER_IMPULSE_ONLY
				m_rigidBody->applyCentralImpulse (impulseVec);
#else
				m_rigidBody->applyImpulse (impulseVec, relative_position);
#endif
			}
		}
	}

	if (columnVolume > btScalar(0.0))
	{
		m_volume += columnVolume;

		if (applyFluidDisplace)
		{
			fluid->addDisplaced (w, l, columnVolume);
		}

		if (applyFluidVelocityImpulse)
		{
			int index = fluid->arrayIndex (w,l);
			btScalar u = fluid->getVelocityXArray()[index];
			btScalar v = fluid->getVelocityZArray()[index];
			btVector3 vd = btVector3(u, btScalar(0.0f), v);
			btVector3 impulse = vd * dt * fluid->getHorizontalVelocityScale();
			m_rigidBody->applyCentralImpulse (impulse);
		}
	}

	return true;
}

