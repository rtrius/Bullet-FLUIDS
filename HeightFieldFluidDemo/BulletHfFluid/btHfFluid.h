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
	int m_numNodesWidth;
	int m_numNodesLength;

	btScalar m_gridCellWidth;
	btScalar m_gridWidth;
	btScalar m_gridLength;

	btScalar m_gridCellWidthInv;
	
	btVector3 m_aabbMin;
	btVector3 m_aabbMax;

	int m_rIndex;
	
	btAlignedObjectArray<btScalar> m_temp;
	btAlignedObjectArray<btScalar> m_height;	//Combined height; fluid + ground
	btAlignedObjectArray<btScalar> m_ground;	//Heightfield
	btAlignedObjectArray<btScalar> m_eta; 		//Depth of fluid; height - ground
	btAlignedObjectArray<btScalar> m_u;			//x velocity
	btAlignedObjectArray<btScalar> m_v;			//z velocity
	btAlignedObjectArray<btScalar> m_r[2];
	btAlignedObjectArray<btScalar> m_fillRatio;
	btAlignedObjectArray<bool> m_flags;

	// tweakables
	btScalar m_globalVelocityU;
	btScalar m_globalVelocityV;
	btScalar m_gravity;
	btScalar m_volumeDisplacementScale;
	btScalar m_horizontalVelocityScale;

	btScalar m_epsHeight;
	btScalar m_epsEta;
	
public:
	btHfFluid (btScalar gridCellWidth, int numNodesWidth, int numNodesLength);
	~btHfFluid ();

	void stepSimulation(btScalar dt);

	///Prep does some initial setup of the height field fluid. Call this at initialization time.
	void prep ();
		
	static const btHfFluid*	upcast(const btCollisionObject* colObj)
	{
		return (colObj->getInternalType() == CO_HF_FLUID) ? (const btHfFluid*)colObj : 0;
	}
	static btHfFluid* upcast(btCollisionObject* colObj)
	{
		return (colObj->getInternalType() == CO_HF_FLUID) ? (btHfFluid*)colObj : 0;
	}
	virtual void getAabb(btVector3& aabbMin,btVector3& aabbMax) const
	{
		aabbMin = m_aabbMin;
		aabbMax = m_aabbMax;
	}

	int getNumNodesWidth () const { return m_numNodesWidth; }
	int getNumNodesLength () const { return m_numNodesLength; }

	btScalar getGridCellWidth () const { return m_gridCellWidth; }
	btScalar widthPos (int i) const { return m_gridCellWidth * i; }
	btScalar lengthPos (int j) const { return m_gridCellWidth * j; }

	int arrayIndex (int i, int j) const;

	const btScalar* getHeightArray () const { return &m_height[0]; }
	const btScalar* getGroundArray () const { return &m_ground[0]; }
	const btScalar* getEtaArray () const { return &m_eta[0]; }
	const btScalar* getVelocityUArray () const { return &m_u[0]; }
	const btScalar* getVelocityVArray () const { return &m_v[0]; }
	const bool* getFlagsArray () const { return &m_flags[0]; }

	btScalar* getHeightArray () {return &m_height[0]; }
	btScalar* getGroundArray () { return &m_ground[0]; }
	btScalar* getEtaArray () { return &m_eta[0]; }
	btScalar* getVelocityUArray () { return &m_u[0]; }
	btScalar* getVelocityVArray () { return &m_v[0]; }
	bool* getFlagsArray() { return &m_flags[0]; }

	void setFluidHeight (int index, btScalar height);

	void addFluidHeight (int x, int y, btScalar height);
	void addDisplaced (int i, int j, btScalar r) { m_r[m_rIndex][arrayIndex(i,j)] += r; }

	void getAabbForColumn (int x, int y, btVector3& aabbMin, btVector3& aabbMax);


	void foreachGroundTriangle(btTriangleCallback* callback,const btVector3& aabbMin,const btVector3& aabbMax);
	
	class btHfFluidColumnCallback 
	{
	public:
		virtual ~btHfFluidColumnCallback () {}

		virtual bool processColumn (btHfFluid* fluid, int w, int l)
		{
			return true; // keep going
		}
	};

	void foreachFluidColumn (btHfFluidColumnCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax);

	void foreachSurfaceTriangle (btTriangleCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax) const;

	///You can enforce a global velocity at the surface of the fluid; default: 0.0 and 0.0
	void setGlobaVelocity (btScalar globalVelocityU, btScalar globalVelocityV)
	{
		m_globalVelocityU = globalVelocityU;
		m_globalVelocityV = globalVelocityV;
	}
	void getGlobalVelocity (btScalar& globalVelocityU, btScalar& globalVelocityV) const
	{
		globalVelocityU = m_globalVelocityU;
		globalVelocityV = m_globalVelocityV;
	}
	
	///Force of gravity, should match physics world; default: -10.0
	void setGravity (btScalar gravity) { m_gravity = gravity; }
	btScalar getGravity () const { return m_gravity; }

	///Percentage of fluid displaced into adjacent cells when a body is submerged; range [0.0, 1.0]; default: 0.5
	void setVolumeDisplacementScale (btScalar volumeDisplacementScale) { m_volumeDisplacementScale = volumeDisplacementScale; }
	btScalar getVolumeDisplacementScale () const { return m_volumeDisplacementScale; }
	
	///Influence of the fluid's horizontal velocity on submerged bodies; range [0.0, 1.0]; default: 0.5
	void setHorizontalVelocityScale (btScalar horizontalVelocityScale) { m_horizontalVelocityScale = horizontalVelocityScale; }
	btScalar getHorizontalVelocityScale () const { return m_horizontalVelocityScale; }
	
protected:
	void setGridDimensions (btScalar gridCellWidth,
							int numNodesWidth, int numNodesLength);

	btScalar bilinearInterpolate (const btAlignedObjectArray<btScalar>& array, btScalar i, btScalar j);

	btScalar advect (const btAlignedObjectArray<btScalar>& array, btScalar i, btScalar j, btScalar di, btScalar dj, btScalar dt);

	void advectEta (btScalar dt);
	void updateHeight (btScalar dt);

	void advectVelocityU (btScalar dt);
	void advectVelocityV (btScalar dt);
	void updateVelocity (btScalar dt);

	void transferDisplaced (btScalar dt);

	void setReflectBoundaryLeft ();
	void setReflectBoundaryRight ();
	void setReflectBoundaryTop ();
	void setReflectBoundaryBottom ();

	void computeFlagsAndFillRatio ();
	btScalar computeHmin (int i, int j);
	btScalar computeHmax (int i, int j);
	btScalar computeEtaMax (int i, int j);

	void allocateArrays ();
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

