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
#ifndef BT_FLUID_HF_H
#define BT_FLUID_HF_H

#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/CollisionShapes/btTriangleCallback.h"

#include "btFluidColumns.h"
#include "btFluidHfSolver.h"
#include "btFluidHfSolverExperimental.h"

class btPersistentManifold;
class btManifoldResult;


// Fix flags and fill ratio
// -> Figure out the constants used in flags and fill ratio code
// Fix volume removal
// add buoyant convex vs. convex / concave
// add buoyant concave support (try bunny model)


class btFluidHf : public btCollisionObject
{
protected:
	btFluidColumns m_columns;
	
	btFluidHfParameters m_hfParameters;

#define USE_NEW_FLUID_HF_SOLVER
#ifndef USE_NEW_FLUID_HF_SOLVER
	btFluidHfSolverDefault m_solver;
#else
	btFluidHfSolverExperimental m_solver;
#endif

	btVector3 m_aabbMin;
	btVector3 m_aabbMax;
	
public:
	btFluidHf(btScalar gridCellWidth, int numNodesX, int numNodesZ);
	~btFluidHf();

	///Prep does some initial setup of the height field fluid. Call this at initialization time.
	void prep();
	
	void stepSimulation(btScalar dt);
	
	int arrayIndex(int i, int j) const { return m_columns.getIndex(i,j); }

	int getNumNodesX() const { return m_columns.m_numNodesX; }
	int getNumNodesZ() const { return m_columns.m_numNodesZ; }

	btScalar getGridCellWidth() const { return m_columns.m_gridCellWidth; }
	btScalar getCellPosX(int i) const { return getGridCellWidth() * i; }
	btScalar getCellPosZ(int j) const { return getGridCellWidth() * j; }
	
	const btAlignedObjectArray<btScalar>& getCombinedHeightArray() const { return m_columns.m_combinedHeight; }
	const btAlignedObjectArray<btScalar>& getGroundArray() const { return m_columns.m_ground; }
	const btAlignedObjectArray<btScalar>& getFluidArray() const { return m_columns.m_fluidDepth; }
	const btAlignedObjectArray<btScalar>& getVelocityXArray() const { return m_columns.m_vel_x; }
	const btAlignedObjectArray<btScalar>& getVelocityZArray() const { return m_columns.m_vel_z; }
	const btAlignedObjectArray<bool>& getFlagsArray() const { return m_columns.m_active; }
	
	btAlignedObjectArray<btScalar>& getCombinedHeightArray() { return m_columns.m_combinedHeight; }
	btAlignedObjectArray<btScalar>& getGroundArray() { return m_columns.m_ground; }
	btAlignedObjectArray<btScalar>& getFluidArray() { return m_columns.m_fluidDepth; }
	btAlignedObjectArray<btScalar>& getVelocityXArray() { return m_columns.m_vel_x; }
	btAlignedObjectArray<btScalar>& getVelocityZArray() { return m_columns.m_vel_z; }
	btAlignedObjectArray<bool>& getFlagsArray() { return m_columns.m_active; }

	void setFluidHeight (int index, btScalar height);

	void addFluidHeight (int x, int y, btScalar height);
	void addDisplaced (int i, int j, btScalar r) { m_columns.m_displaced[m_columns.m_displacedIndex][arrayIndex(i,j)] += r; }

	void getAabbForColumn (int x, int y, btVector3& aabbMin, btVector3& aabbMax);
	
	class btFluidHfColumnCallback 
	{
	public:
		virtual ~btFluidHfColumnCallback () {}

		///btFluidHf::forEachFluidColumn() will continue calling processColumn() if this returns true
		virtual bool processColumn (btFluidHf* fluid, int w, int l) { return true; }
	};

	void forEachFluidColumn(btFluidHfColumnCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax);
	void forEachGroundTriangle(btTriangleCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax) const
	{
		forEachTriangle(callback, aabbMin, aabbMax, m_columns.m_ground);
	}
	void forEachSurfaceTriangle(btTriangleCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax) const
	{
		forEachTriangle(callback, aabbMin, aabbMax, m_columns.m_combinedHeight);
	}

	void setGlobalVelocity(btScalar globalVelocityX, btScalar globalVelocityZ)
	{
		m_hfParameters.m_globalVelocityX = globalVelocityX;
		m_hfParameters.m_globalVelocityZ = globalVelocityZ;
	}
	void getGlobalVelocity(btScalar& globalVelocityX, btScalar& globalVelocityZ) const
	{
		globalVelocityX = m_hfParameters.m_globalVelocityX;
		globalVelocityZ = m_hfParameters.m_globalVelocityZ;
	}
	
	void setGravity (btScalar gravity) { m_hfParameters.m_gravity = gravity; }
	btScalar getGravity () const { return m_hfParameters.m_gravity; }
	
	void setVolumeDisplacementScale (btScalar volumeDisplacementScale) { m_hfParameters.m_volumeDisplacementScale = volumeDisplacementScale; }
	btScalar getVolumeDisplacementScale () const { return m_hfParameters.m_volumeDisplacementScale; }
	
	void setHorizontalVelocityScale (btScalar horizontalVelocityScale) { m_hfParameters.m_horizontalVelocityScale = horizontalVelocityScale; }
	btScalar getHorizontalVelocityScale () const { return m_hfParameters.m_horizontalVelocityScale; }
	
	virtual void getAabb(btVector3& aabbMin, btVector3& aabbMax) const
	{
		aabbMin = m_aabbMin;
		aabbMax = m_aabbMax;
	}
	
	//btCollisionObject
	virtual void setCollisionShape(btCollisionShape *collisionShape) { btAssert(0); }
	
	static const btFluidHf*	upcast(const btCollisionObject* colObj)
	{
		return (colObj->getInternalType() == CO_HF_FLUID) ? (const btFluidHf*)colObj : 0;
	}
	static btFluidHf* upcast(btCollisionObject* colObj)
	{
		return (colObj->getInternalType() == CO_HF_FLUID) ? (btFluidHf*)colObj : 0;
	}
	
protected:
	void forEachTriangle(btTriangleCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax, 
						const btAlignedObjectArray<btScalar>& heightArray) const;
};



#endif

