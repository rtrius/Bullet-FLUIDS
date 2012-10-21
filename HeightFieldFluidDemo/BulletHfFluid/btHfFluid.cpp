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

#include <stdio.h>
#include "btHfFluid.h"
#include "btHfFluidCollisionShape.h"
#include "btHfFluidBuoyantConvexShape.h"
#include "BulletCollision/NarrowPhaseCollision/btPersistentManifold.h"
#include "BulletCollision/CollisionDispatch/btConvexConcaveCollisionAlgorithm.h"
#include "BulletCollision/CollisionDispatch/btCollisionWorld.h"
#include "BulletCollision/CollisionDispatch/btManifoldResult.h"
#include "BulletCollision/CollisionShapes/btTriangleShape.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"
#include "../../OpenGL/GLDebugDrawer.h"
								
btHfFluid::btHfFluid (btScalar gridCellWidth, int numNodesWidth, int numNodesLength)
{
	m_internalType = CO_HF_FLUID;
	m_rIndex = 0;
	setGridDimensions (gridCellWidth, numNodesWidth, numNodesLength);

	btScalar maxHeight = 20.0;
	m_aabbMin = btVector3 (0.0, 0.0, 0.0);
	m_aabbMax = btVector3 (m_gridWidth, maxHeight, m_gridLength);

	setCollisionShape( new btHfFluidCollisionShape(this) );

	m_globalVelocityU = btScalar(0.0f);
	m_globalVelocityV = btScalar(0.0f);
	m_gravity = btScalar(-10.0f);

	m_volumeDisplacementScale = btScalar(0.5f);
	m_horizontalVelocityScale = btScalar(0.5f);

	m_epsEta = btScalar(0.01f);
	m_epsHeight = btScalar(0.001f);
}

btHfFluid::~btHfFluid ()
{
	btCollisionShape* collisionShape = getCollisionShape();
	delete collisionShape;
}

void btHfFluid::stepSimulation(btScalar dt)
{
	transferDisplaced (dt);
	
	advectEta (dt);
	advectVelocityU (dt);
	advectVelocityV (dt);
	updateHeight (dt);
	computeFlagsAndFillRatio ();
	updateVelocity (dt);
	
	setReflectBoundaryLeft ();
	setReflectBoundaryRight ();
	setReflectBoundaryTop ();
	setReflectBoundaryBottom ();
	
	{
		static btScalar total_volume = btScalar(0.0f);
		btScalar new_total_volume = btScalar(0.0f);
		for (int i = 0; i < m_numNodesWidth*m_numNodesLength; i++)
		{
			new_total_volume += m_eta[i] * m_gridCellWidth * m_gridCellWidth;
		}
		printf("volume = %f volume delta = %f\n", new_total_volume, new_total_volume - total_volume);
		total_volume = new_total_volume;
	}
}

void btHfFluid::prep ()
{
	for(int i = 0; i < m_numNodesLength*m_numNodesWidth; i++) m_height[i] = m_eta[i] + m_ground[i];
	computeFlagsAndFillRatio ();
}


int btHfFluid::arrayIndex (int i, int j) const
{
	btAssert (i >= 0);
	btAssert (i < m_numNodesWidth);
	btAssert (j >= 0);
	btAssert (j < m_numNodesLength);
	int index = i + (j * m_numNodesWidth);
	return index;
}


void btHfFluid::setFluidHeight (int index, btScalar height)
{
	m_eta[index] = height;
	m_height[index] = m_ground[index] + m_eta[index];
	m_flags[index] = true;
}

void btHfFluid::addFluidHeight (int x, int y, btScalar height)
{
	int index = arrayIndex (x,y);
	m_eta[index] += height;
	m_height[index] = m_ground[index] + m_eta[index];
	m_flags[index] = true;
}

void btHfFluid::getAabbForColumn (int i, int j, btVector3& aabbMin, btVector3& aabbMax)
{
	btVector3 com = getWorldTransform().getOrigin();
	int sw = arrayIndex (i, j);
	int se = arrayIndex (i+1, j);
	int nw = arrayIndex (i, j+1);
	int ne = arrayIndex (i+1, j+1);

	btScalar h = m_height[sw];
	btScalar g = m_ground[sw];

	aabbMin = btVector3(widthPos (i), g, lengthPos (j));
	aabbMax = btVector3(widthPos (i+1), h, lengthPos (j+1));
	aabbMin += com;
	aabbMax += com;
}

void btHfFluid::foreachGroundTriangle(btTriangleCallback* callback,const btVector3& aabbMin,const btVector3& aabbMax)
{
	btVector3 verts[3];

	btScalar minX, minZ, maxX, maxZ;
	int startNodeX, startNodeZ, endNodeX, endNodeZ;

	minX = aabbMin.getX();
	minZ = aabbMin.getZ();
	maxX = aabbMax.getX();
	maxZ = aabbMax.getZ();

	startNodeX = (int)(minX * m_gridCellWidthInv);
	startNodeZ = (int)(minZ * m_gridCellWidthInv);

	endNodeX = (int)(maxX * m_gridCellWidthInv);
	endNodeZ = (int)(maxZ * m_gridCellWidthInv);

	endNodeX++;
	endNodeZ++;

	startNodeX = btMax (1, startNodeX);
	startNodeZ = btMax (1, startNodeZ);
	endNodeX = btMin (m_numNodesWidth-1, endNodeX);
	endNodeZ = btMin (m_numNodesLength-1, endNodeZ);

#ifdef __BRUTE__
	for (int i = 1; i < m_numNodesWidth-1; i++)
	{
		for (int j = 1; j < m_numNodesLength-1; j++)
		{
			// triangle 1
			verts[0] = btVector3(widthPos(i), m_ground[arrayIndex(i,j)], lengthPos(j));
			verts[1] = btVector3(widthPos(i), m_ground[arrayIndex(i,j+1)], lengthPos(j+1));
			verts[2] = btVector3(widthPos(i+1), m_ground[arrayIndex(i+1,j)], lengthPos(j));
			callback->processTriangle(verts,i,j);
			// triangle 2
			verts[0] = btVector3(widthPos(i+1), m_ground[arrayIndex(i+1,j)], lengthPos(j));
			verts[1] = btVector3(widthPos(i), m_ground[arrayIndex(i,j+1)], lengthPos(j+1));
			verts[2] = btVector3(widthPos(i+1), m_ground[arrayIndex(i+1,j+1)], lengthPos(j+1));
			callback->processTriangle(verts,i,j);
		}
	}
#else

	for (int i = startNodeX; i < endNodeX; i++)
	{
		for (int j = startNodeZ; j < endNodeZ; j++)
		{
			// triangle 1
			verts[0] = btVector3(widthPos(i), m_ground[arrayIndex(i,j)], lengthPos(j));
			verts[1] = btVector3(widthPos(i), m_ground[arrayIndex(i,j+1)], lengthPos(j+1));
			verts[2] = btVector3(widthPos(i+1), m_ground[arrayIndex(i+1,j)], lengthPos(j));
			callback->processTriangle(verts,i,j);
			// triangle 2
			verts[0] = btVector3(widthPos(i+1), m_ground[arrayIndex(i+1,j)], lengthPos(j));
			verts[1] = btVector3(widthPos(i), m_ground[arrayIndex(i,j+1)], lengthPos(j+1));
			verts[2] = btVector3(widthPos(i+1), m_ground[arrayIndex(i+1,j+1)], lengthPos(j+1));
			callback->processTriangle(verts,i,j);
		}
	}
#endif
}

void btHfFluid::foreachFluidColumn (btHfFluidColumnCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax)
{
	btScalar minX, minZ, maxX, maxZ;
	int startNodeX, startNodeZ, endNodeX, endNodeZ;

	minX = aabbMin.getX();
	minZ = aabbMin.getZ();
	maxX = aabbMax.getX();
	maxZ = aabbMax.getZ();

	startNodeX = (int)(minX * m_gridCellWidthInv);
	startNodeZ = (int)(minZ * m_gridCellWidthInv);

	endNodeX = (int)(maxX * m_gridCellWidthInv);
	endNodeZ = (int)(maxZ * m_gridCellWidthInv);

	endNodeX++;
	endNodeZ++;

	startNodeX = btMax (1, startNodeX);
	startNodeZ = btMax (1, startNodeZ);
	endNodeX = btMin (m_numNodesWidth-2, endNodeX);
	endNodeZ = btMin (m_numNodesLength-2, endNodeZ);

	bool r;
	for (int i = startNodeX; i < endNodeX; i++)
	{
		for (int j = startNodeZ; j < endNodeZ; j++)
		{
			if ( !m_flags[arrayIndex(i, j)] )
				continue;

			r = callback->processColumn (this, i, j);
			if (!r)
				return;
		}
	}
}

void btHfFluid::foreachSurfaceTriangle (btTriangleCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax) const
{
	btVector3 verts[3];

	btScalar minX, minZ, maxX, maxZ;
	int startNodeX, startNodeZ, endNodeX, endNodeZ;

	minX = aabbMin.getX();
	minZ = aabbMin.getZ();
	maxX = aabbMax.getX();
	maxZ = aabbMax.getZ();

	startNodeX = (int)(minX * m_gridCellWidthInv);
	startNodeZ = (int)(minZ * m_gridCellWidthInv);

	endNodeX = (int)(maxX * m_gridCellWidthInv);
	endNodeZ = (int)(maxZ * m_gridCellWidthInv);

	endNodeX++;
	endNodeZ++;

	startNodeX = btMax (1, startNodeX);
	startNodeZ = btMax (1, startNodeZ);
	endNodeX = m_numNodesWidth-1;
	endNodeZ = m_numNodesLength-1;

	for (int i = startNodeX; i < endNodeX; i++)
	{
		for (int j = startNodeZ; j < endNodeZ; j++)
		{
			if (!m_flags[arrayIndex(i,j)])
				continue;
			// triangle 1
			verts[0] = btVector3(widthPos(i), m_height[arrayIndex(i,j)], lengthPos(j));
			verts[1] = btVector3(widthPos(i), m_height[arrayIndex(i,j+1)], lengthPos(j+1));
			verts[2] = btVector3(widthPos(i+1), m_height[arrayIndex(i+1,j)], lengthPos(j));
			callback->processTriangle(verts,i,j);
			// triangle 2
			verts[0] = btVector3(widthPos(i+1), m_height[arrayIndex(i+1,j)], lengthPos(j));
			verts[1] = btVector3(widthPos(i), m_height[arrayIndex(i,j+1)], lengthPos(j+1));
			verts[2] = btVector3(widthPos(i+1), m_height[arrayIndex(i+1,j+1)], lengthPos(j+1));
			callback->processTriangle(verts,i,j);
		}
	}
}

void btHfFluid::setGridDimensions (btScalar gridCellWidth,
								   int numNodesWidth, int numNodesLength)
{
	m_gridWidth = gridCellWidth * numNodesWidth;
	m_gridLength = gridCellWidth * numNodesLength;
	m_gridCellWidth = gridCellWidth;
	m_numNodesWidth = numNodesWidth;
	m_numNodesLength = numNodesLength;
	m_gridCellWidthInv = btScalar(1.0) / gridCellWidth;

	allocateArrays ();
}

btScalar btHfFluid::bilinearInterpolate (const btAlignedObjectArray<btScalar>& array, btScalar iPos, btScalar jPos)
{
	int i = (int)iPos;
	int j = (int)jPos;

	btScalar iParam1 = iPos - btScalar(i);
	btScalar iParam0 = btScalar(1.0) - iParam1;

	btScalar jParam1 = jPos - btScalar(j);
	btScalar jParam0 = btScalar(1.0) - jParam1;

	btScalar SW = array[arrayIndex(i, j)];
	btScalar SE = array[arrayIndex(i+1, j)];
	btScalar NW = array[arrayIndex(i, j+1)];
	btScalar NE = array[arrayIndex(i+1, j+1)];

	btScalar a = jParam0 * SW + jParam1 * NW;
	btScalar b = jParam0 * SE + jParam1 * NE;
	return iParam0 * a + iParam1 * b;
}

btScalar btHfFluid::advect (const btAlignedObjectArray<btScalar>& array, btScalar i, btScalar j, btScalar di, btScalar dj,btScalar dt)
{
	// trace particle backwards in time
	btScalar srcI = i - di * dt * m_gridCellWidthInv;
	btScalar srcJ = j - dj * dt * m_gridCellWidthInv;

	// srcI and srcJ are indices into the array,
	// we need to clamp them to within the domain
	srcI = btMax (srcI, btScalar(0.0));
	srcI = btMin (srcI, btScalar(m_numNodesWidth-1));
	srcJ = btMax (srcJ, btScalar(0.0));
	srcJ = btMin (srcJ, btScalar(m_numNodesLength-1));

	return bilinearInterpolate (array, srcI, srcJ);
}

void btHfFluid::advectEta (btScalar dt)
{
	for (int i = 1; i < m_numNodesWidth-1; i++)
	{
		for (int j = 1; j < m_numNodesLength-1; j++)
		{
			int index = arrayIndex (i, j);

			if (!m_flags[index])
				continue;

			btScalar u = m_globalVelocityU;
			btScalar v = m_globalVelocityV;

			u += (m_u[index]+m_u[index+1]) * btScalar(0.5);
			v += (m_v[index]+m_v[index+m_numNodesWidth]) * btScalar(0.5);

			m_temp[index] = advect (m_eta, btScalar(i), btScalar(j), u, v, dt);
		}
	}
	for (int i = 1; i < m_numNodesWidth-1; i++)
	{
		for (int j = 1; j < m_numNodesLength-1; j++)
		{
			int index = arrayIndex (i, j);
			m_eta[index] = m_temp[index];
		}
	}
}

void btHfFluid::updateHeight (btScalar dt)
{
	for (int j = 1; j < m_numNodesLength-1; j++)
	{
		for (int i = 1; i < m_numNodesWidth-1; i++)
		{
			int index = arrayIndex (i, j);
			if (!m_flags[index])
			{
				m_height[index] = m_ground[index] + m_eta[index];
				continue;
			}
			btScalar deta = -btScalar(0.5f) * m_eta[index] * dt * m_gridCellWidthInv * ( (m_u[index+1] - m_u[index]) + (m_v[index+m_numNodesWidth] - m_v[index]));
			m_eta[index] += deta;
			m_height[index] = m_ground[index] + btMax(m_eta[index],btScalar(0.0f));
		}
	}
}

void btHfFluid::advectVelocityU (btScalar dt)
{
	for (int i = 1; i < m_numNodesWidth-1; i++)
	{
		for (int j = 1; j < m_numNodesLength-1; j++)
		{
			int index = arrayIndex (i, j);
			if (!m_flags[index])
			{
				continue;
			}
			btScalar u = m_globalVelocityU;
			btScalar v = m_globalVelocityV;

			u += m_u[index];
			v += (m_v[index]+m_v[index+1]+m_v[index+m_numNodesWidth]+m_v[index+m_numNodesWidth+1]) * btScalar(0.25);

			m_temp[index] = advect (m_u, btScalar(i), btScalar(j), u, v, dt);
		}
	}

	for (int i = 1; i < m_numNodesWidth-1; i++)
	{
		for (int j = 1; j < m_numNodesLength-1; j++)
		{
			int index = arrayIndex (i, j);
			m_u[index] = m_temp[index];
		}
	}
}

void btHfFluid::advectVelocityV (btScalar dt)
{
	for (int i = 1; i < m_numNodesWidth-1; i++)
	{
		for (int j = 1; j < m_numNodesLength-1; j++)
		{
			int index = arrayIndex (i, j);
			if (!m_flags[index])
			{
				continue;
			}
			btScalar u = m_globalVelocityU;
			btScalar v = m_globalVelocityV;

			u += (m_u[index]+m_u[index+1]+m_u[index+m_numNodesWidth]+m_u[index+m_numNodesWidth+1]) * btScalar(0.25);
			v += m_v[index];
			

			m_temp[index] = advect (m_v, btScalar(i), btScalar(j), u, v, dt);
		}
	}
	for (int i = 1; i < m_numNodesWidth-1; i++)
	{
		for (int j = 1; j < m_numNodesLength-1; j++)
		{
			int index = arrayIndex (i, j);
			m_v[index] = m_temp[index];
		}
	}
}


void btHfFluid::transferDisplaced (btScalar dt)
{
	for (int i = 2; i < m_numNodesWidth - 2; i++)
	{
		for (int j = 2; j < m_numNodesLength - 2; j++)
		{
			btScalar deltaR = m_r[m_rIndex][arrayIndex(i,j)] - m_r[(m_rIndex+1)%2][arrayIndex(i,j)];
			deltaR /= m_gridCellWidth * m_gridCellWidth;	//deltaR is in volume, but we want to change the height
			deltaR *= m_volumeDisplacementScale;
			
			btScalar qdeltaR = deltaR / btScalar(4.0f);
			m_eta[arrayIndex(i-1,j-1)] += qdeltaR;
			m_eta[arrayIndex(i-1,j+1)] += qdeltaR;
			m_eta[arrayIndex(i+1,j-1)] += qdeltaR;
			m_eta[arrayIndex(i+1,j+1)] += qdeltaR;
			m_eta[arrayIndex(i,j)] -= deltaR;
			
			// OPTIMIZATION: zero out next frames r value
			m_r[(m_rIndex+1)%2][arrayIndex(i,j)] = btScalar(0.0);
		}
	}
	
	m_rIndex = (m_rIndex + 1) % 2; // flip frame
}

void btHfFluid::updateVelocity (btScalar dt)
{
	for (int j = 1; j < m_numNodesLength-1; j++)
	{
		for (int i = 2; i < m_numNodesWidth-1; i++)
		{
			int index = arrayIndex (i, j);
			if (!m_flags[index])
			{
				continue;
			}
			m_u[index] += m_gravity * dt * m_gridCellWidthInv * (m_height[index]-m_height[index-1]);
		}
	}

	for (int j = 2; j < m_numNodesLength-1; j++)
	{
		for (int i = 1; i < m_numNodesWidth-1; i++)
		{
			int index = arrayIndex (i, j);
			if (!m_flags[index])
			{
				continue;
			}
			m_v[index] += m_gravity * dt * m_gridCellWidthInv * (m_height[index]-m_height[index-m_numNodesWidth]);
		}
	}
}

void btHfFluid::setReflectBoundaryLeft ()
{
	for (int j = 0; j < m_numNodesLength; j++)
	{
		int indexL = arrayIndex (0, j);

		m_height[indexL] = m_height[indexL+1];
		m_u[indexL+1] = btScalar(0.0);
		m_v[indexL] = btScalar(0.0);
	}
}

void btHfFluid::setReflectBoundaryRight ()
{
	for (int j = 0; j < m_numNodesLength; j++)
	{
		int indexR = arrayIndex (m_numNodesWidth-1, j);

		m_height[indexR] = m_height[indexR-1];
		m_u[indexR-1] = btScalar(0.0);
		m_v[indexR] = btScalar(0.0);
	}
}

void btHfFluid::setReflectBoundaryBottom ()
{
	for (int i = 0; i < m_numNodesWidth; i++)
	{
		int indexT = arrayIndex (i, 0);
		
		m_height[indexT] = m_height[indexT+m_numNodesWidth];
		m_v[indexT+m_numNodesWidth] = btScalar(0.0);
		m_u[indexT] = btScalar(0.0);
	}
}

void btHfFluid::setReflectBoundaryTop ()
{
	for (int i = 0; i < m_numNodesWidth; i++)
	{
		int indexB = arrayIndex (i, m_numNodesLength-1);

		m_height[indexB] = m_height[indexB-m_numNodesWidth];
		m_v[indexB-m_numNodesWidth] = btScalar(0.0);
		m_u[indexB] = btScalar(0.0);
	}
}

/*
void btHfFluid::setAbsorbBoundaryLeft (btScalar dt)
{
	for (int j = 0; j < m_numNodesLength; j++)
	{
		int indexL = arrayIndex (0, j);

		btScalar c = btSqrt(m_eta[indexL+1]*m_gravity);
		m_height[indexL] = ((m_gridCellWidthInv * m_height[(m_heightIndex+1)%2][indexL+1])+(dt*c*m_height[indexL+1]))/(m_gridCellWidthInv + dt * c);
		m_u[indexL+1] = btScalar(0.0);
		m_v[indexL+1] = btScalar(0.0);
	}
}
*/


void btHfFluid::computeFlagsAndFillRatio ()
{
	for (int i = 1; i < m_numNodesWidth-1; i++)
	{
		for (int j = 1; j < m_numNodesLength-1; j++)
		{
			btScalar h = m_height[arrayIndex(i,j)];
			btScalar hMin = computeHmin(i,j);
			btScalar hMax = computeHmax(i,j);
			btScalar etaMax = computeEtaMax(i,j);
			if (h <= hMin && etaMax < m_epsEta)
			{
				m_flags[arrayIndex(i,j)] = false;
				m_fillRatio[arrayIndex(i,j)] = btScalar(0.0f);
			} 
			else if (h > hMax)
			{
				m_flags[arrayIndex(i,j)] = true;
				m_fillRatio[arrayIndex(i,j)] = btScalar(1.0f);
			} 
			else
			{
				m_flags[arrayIndex(i,j)] = true;
				m_fillRatio[arrayIndex(i,j)] = (h - hMin)/(hMax - hMin);
			}
			
		}
	}
}

btScalar btHfFluid::computeHmin (int i, int j)
{
	btAssert (i > 0);
	btAssert (i < m_numNodesWidth-1);
	btAssert (j > 0);
	btAssert (j < m_numNodesLength-1);

	btScalar h1 = m_height[arrayIndex(i-1,j-1)];
	btScalar h2 = m_height[arrayIndex(i-1,j+1)];
	btScalar h3 = m_height[arrayIndex(i+1,j-1)];
	btScalar h4 = m_height[arrayIndex(i+1,j+1)];
	btScalar h = m_height[arrayIndex(i,j)];
	btScalar minh = btMin(h1, btMin(h2, btMin(h3,h4)));

	return (minh + h) * btScalar(0.5f);
}

btScalar btHfFluid::computeHmax (int i, int j)
{
	btAssert (i > 0);
	btAssert (i < m_numNodesWidth-1);
	btAssert (j > 0);
	btAssert (j < m_numNodesLength-1);

	btScalar h1 = m_height[arrayIndex(i-1,j-1)];
	btScalar h2 = m_height[arrayIndex(i-1,j+1)];
	btScalar h3 = m_height[arrayIndex(i+1,j-1)];
	btScalar h4 = m_height[arrayIndex(i+1,j+1)];
	btScalar h = m_height[arrayIndex(i,j)];
	btScalar maxh = btMax(h1, btMax(h2, btMax(h3,h4)));

	return (maxh + h) * btScalar(0.5f) + m_epsHeight;
}

btScalar btHfFluid::computeEtaMax (int i, int j)
{
	btAssert (i > 0);
	btAssert (i < m_numNodesWidth-1);
	btAssert (j > 0);
	btAssert (j < m_numNodesLength-1);

	btScalar eta1 = m_eta[arrayIndex(i-1,j-1)];
	btScalar eta2 = m_eta[arrayIndex(i-1,j+1)];
	btScalar eta3 = m_eta[arrayIndex(i+1,j-1)];
	btScalar eta4 = m_eta[arrayIndex(i+1,j+1)];
	btScalar eta = m_eta[arrayIndex(i,j)];
	btScalar maxeta = btMax(eta1, btMax(eta2, btMax(eta3,eta4)));

	return (maxeta + eta) * btScalar(0.5f);
}

void btHfFluid::allocateArrays ()
{
	int numNodes = m_numNodesWidth * m_numNodesLength;
	m_temp.resize(numNodes);
	m_height.resize(numNodes);
	m_ground.resize(numNodes);
	m_eta.resize(numNodes);
	m_u.resize(numNodes);
	m_v.resize(numNodes);
	m_r[0].resize(numNodes);
	m_r[1].resize(numNodes);
	m_fillRatio.resize(numNodes);
	m_flags.resize(numNodes);

	for (int i = 0; i < numNodes; i++)
	{
		m_temp[i] = btScalar(0.0);
		m_height[i] = btScalar(0.0);
		m_eta[i] = btScalar(0.0);
		m_u[i] = btScalar(0.0);
		m_v[i] = btScalar(0.0);
		m_r[0][i] = btScalar(0.0);
		m_r[1][i] = btScalar(0.0);
		m_ground[i] = btScalar(0.0);
		m_fillRatio[i] = btScalar(0.0);
		m_flags[i] = false;
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

btHfFluidColumnRigidBodyCallback::btHfFluidColumnRigidBodyCallback (btRigidBody* rigidBody, btIDebugDraw* debugDraw, btScalar density, btScalar floatyness)
{
	m_rigidBody = rigidBody;
	m_buoyantShape = (btHfFluidBuoyantConvexShape*)rigidBody->getCollisionShape();
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

btHfFluidColumnRigidBodyCallback::~btHfFluidColumnRigidBodyCallback ()
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

bool btHfFluidColumnRigidBodyCallback::processColumn (btHfFluid* fluid, int w, int l)
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
			btScalar u = fluid->getVelocityUArray()[index];
			btScalar v = fluid->getVelocityVArray()[index];
			btVector3 vd = btVector3(u, btScalar(0.0f), v);
			btVector3 impulse = vd * dt * fluid->getHorizontalVelocityScale();
			m_rigidBody->applyCentralImpulse (impulse);
		}
	}

	return true;
}

