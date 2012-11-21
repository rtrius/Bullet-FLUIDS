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
#ifndef BT_HF_FLUID_H
#define BT_HF_FLUID_H

#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletCollision/CollisionShapes/btTriangleCallback.h"

#include "btFluidColumns.h"
#include "btFluidHfSolver.h"

class btPersistentManifold;
class btManifoldResult;

// FIX AABB calculation for whole btHfFluid shape
// Fix flags and fill ratio
// -> Figure out the constants used in flags and fill ratio code
// Fix volume removal
// add buoyant convex vs. convex / concave
// add buoyant concave support (try bunny model)


class btHfFluid : public btCollisionObject
{
protected:
	btFluidColumns m_columns;
	
	btFluidHfParameters m_hfParameters;
	
	btFluidHfSolverDefault m_solver;
	
	btVector3 m_aabbMin;
	btVector3 m_aabbMax;
	
public:
	btHfFluid (btScalar gridCellWidth, int numNodesWidth, int numNodesLength);
	~btHfFluid ();

	///Prep does some initial setup of the height field fluid. Call this at initialization time.
	void prep ();
	
	void stepSimulation(btScalar dt);
	
	int arrayIndex (int i, int j) const { return m_columns.getIndex(i,j); }

	int getNumNodesWidth () const { return m_columns.m_numNodesWidth; }
	int getNumNodesLength () const { return m_columns.m_numNodesLength; }

	btScalar getGridCellWidth () const { return m_columns.m_gridCellWidth; }
	btScalar widthPos (int i) const { return m_columns.m_gridCellWidth * i; }
	btScalar lengthPos (int j) const { return m_columns.m_gridCellWidth * j; }
	btScalar getTotalWidth() const { return m_columns.m_gridCellWidth * static_cast<btScalar>(m_columns.m_numNodesWidth); }
	btScalar getTotalLength() const { return m_columns.m_gridCellWidth * static_cast<btScalar>(m_columns.m_numNodesLength); }
	
	const btAlignedObjectArray<btScalar>& getHeightArray() const { return m_columns.m_height; }
	const btAlignedObjectArray<btScalar>& getGroundArray() const { return m_columns.m_ground; }
	const btAlignedObjectArray<btScalar>& getEtaArray() const { return m_columns.m_eta; }
	const btAlignedObjectArray<btScalar>& getVelocityUArray() const { return m_columns.m_u; }
	const btAlignedObjectArray<btScalar>& getVelocityVArray() const { return m_columns.m_v; }
	const btAlignedObjectArray<bool>& getFlagsArray() const { return m_columns.m_flags; }
	
	btAlignedObjectArray<btScalar>& getHeightArray() { return m_columns.m_height; }
	btAlignedObjectArray<btScalar>& getGroundArray() { return m_columns.m_ground; }
	btAlignedObjectArray<btScalar>& getEtaArray() { return m_columns.m_eta; }
	btAlignedObjectArray<btScalar>& getVelocityUArray() { return m_columns.m_u; }
	btAlignedObjectArray<btScalar>& getVelocityVArray() { return m_columns.m_v; }
	btAlignedObjectArray<bool>& getFlagsArray() { return m_columns.m_flags; }

	void setFluidHeight (int index, btScalar height);

	void addFluidHeight (int x, int y, btScalar height);
	void addDisplaced (int i, int j, btScalar r) { m_columns.m_r[m_columns.m_rIndex][arrayIndex(i,j)] += r; }

	void getAabbForColumn (int x, int y, btVector3& aabbMin, btVector3& aabbMax);
	
	class btHfFluidColumnCallback 
	{
	public:
		virtual ~btHfFluidColumnCallback () {}

		///btHfFluid::forEachFluidColumn() will continue calling processColumn() if this returns true
		virtual bool processColumn (btHfFluid* fluid, int w, int l) { return true; }
	};

	void forEachFluidColumn(btHfFluidColumnCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax);
	void forEachGroundTriangle(btTriangleCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax) const
	{
		forEachTriangle(callback, aabbMin, aabbMax, m_columns.m_ground);
	}
	void forEachSurfaceTriangle(btTriangleCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax) const
	{
		forEachTriangle(callback, aabbMin, aabbMax, m_columns.m_height);
	}

	///You can enforce a global velocity at the surface of the fluid; default: 0.0 and 0.0
	void setGlobalVelocity (btScalar globalVelocityU, btScalar globalVelocityV)
	{
		m_hfParameters.m_globalVelocityU = globalVelocityU;
		m_hfParameters.m_globalVelocityV = globalVelocityV;
	}
	void getGlobalVelocity (btScalar& globalVelocityU, btScalar& globalVelocityV) const
	{
		globalVelocityU = m_hfParameters.m_globalVelocityU;
		globalVelocityV = m_hfParameters.m_globalVelocityV;
	}
	
	///Force of gravity, should match physics world; default: -9.8
	void setGravity (btScalar gravity) { m_hfParameters.m_gravity = gravity; }
	btScalar getGravity () const { return m_hfParameters.m_gravity; }

	///Percentage of fluid displaced into adjacent cells when a body is submerged; range [0.0, 1.0]; default: 0.5
	void setVolumeDisplacementScale (btScalar volumeDisplacementScale) { m_hfParameters.m_volumeDisplacementScale = volumeDisplacementScale; }
	btScalar getVolumeDisplacementScale () const { return m_hfParameters.m_volumeDisplacementScale; }
	
	///Influence of the fluid's horizontal velocity on submerged bodies; range [0.0, 1.0]; default: 0.5
	void setHorizontalVelocityScale (btScalar horizontalVelocityScale) { m_hfParameters.m_horizontalVelocityScale = horizontalVelocityScale; }
	btScalar getHorizontalVelocityScale () const { return m_hfParameters.m_horizontalVelocityScale; }
	
	virtual void getAabb(btVector3& aabbMin, btVector3& aabbMax) const
	{
		aabbMin = m_aabbMin;
		aabbMax = m_aabbMax;
	}
	
	//btCollisionObject
	virtual void setCollisionShape(btCollisionShape *collisionShape) { btAssert(0); }
	
	static const btHfFluid*	upcast(const btCollisionObject* colObj)
	{
		return (colObj->getInternalType() == CO_HF_FLUID) ? (const btHfFluid*)colObj : 0;
	}
	static btHfFluid* upcast(btCollisionObject* colObj)
	{
		return (colObj->getInternalType() == CO_HF_FLUID) ? (btHfFluid*)colObj : 0;
	}
	
protected:
	void forEachTriangle(btTriangleCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax, 
						const btAlignedObjectArray<btScalar>& heightArray) const;
};

class btRigidBody;
class btIDebugDraw;
class btHfFluidBuoyantConvexShape;

class btHfFluidColumnRigidBodyCallback : public btHfFluid::btHfFluidColumnCallback
{
protected:
	btRigidBody* m_rigidBody;
	btHfFluidBuoyantConvexShape* m_buoyantShape;
	btIDebugDraw* m_debugDraw;
	
	int m_numVoxels;
	btVector3* m_voxelPositionsXformed;
	bool* m_voxelSubmerged;
	
	btVector3 m_aabbMin;
	btVector3 m_aabbMax;
	btScalar m_volume;
	btScalar m_density;
	btScalar m_floatyness;
public:
	btHfFluidColumnRigidBodyCallback (btRigidBody* rigidBody, btIDebugDraw* debugDraw, btScalar density, btScalar floatyness);
	~btHfFluidColumnRigidBodyCallback ();
	bool processColumn (btHfFluid* fluid, int w, int l);
	btScalar getVolume () const { return m_volume; }
};

#endif

