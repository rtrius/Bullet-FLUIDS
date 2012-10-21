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
#ifndef __BT_HFFLUID_BUOYANT_CONVEX_SHAPE_H
#define __BT_HFFLUID_BUOYANT_CONVEX_SHAPE_H

#include "LinearMath/btVector3.h"
#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionShapes/btConvexShape.h"


///Wrapper for btConvexShape that allows it to interact with a btHfFluid by adding a voxelized representation.
class btHfFluidBuoyantConvexShape : public btCollisionShape
{
protected:
	btScalar m_floatyness;
	btScalar m_radius;
	btScalar m_totalVolume;
	btScalar m_volumePerVoxel;
	int m_numVoxels;
	btVector3* m_voxelPositions;
	btConvexShape* m_convexShape;
	
public:
	btHfFluidBuoyantConvexShape (btConvexShape* convexShape);
	
	virtual ~btHfFluidBuoyantConvexShape () { if(m_voxelPositions) btAlignedFree (m_voxelPositions); }
	
	void generateShape (btScalar radius, btScalar gap);

	btConvexShape* getConvexShape () { return m_convexShape; }
	const btConvexShape* getConvexShape() const { return m_convexShape; }

	virtual void getAabb(const btTransform& t,btVector3& aabbMin,btVector3& aabbMax) const 
	{ 
		return m_convexShape->getAabb (t, aabbMin, aabbMax);
	}
	virtual void setMargin(btScalar margin) { m_convexShape->setMargin (margin); }
	virtual btScalar getMargin() const { return m_convexShape->getMargin(); }
	virtual void setLocalScaling(const btVector3& scaling) { m_convexShape->setLocalScaling (scaling); }
	virtual const btVector3& getLocalScaling() const { return m_convexShape->getLocalScaling(); }
	virtual void calculateLocalInertia(btScalar mass, btVector3& inertia) const
	{ 
		m_convexShape->calculateLocalInertia (mass, inertia); 
	}
	virtual const char*	getName() const { return "HF_FLUID_BUOYANT_CONVEX_SHAPE"; }

	btScalar getVoxelRadius () const { return m_radius; }
	btScalar getTotalVolume () const { return m_totalVolume; }
	btScalar getVolumePerVoxel () const { return m_volumePerVoxel; }
	btScalar getFloatyness () const { return m_floatyness; }
	void setFloatyness (btScalar floatyness) { m_floatyness = floatyness; }
	int getNumVoxels () const { return m_numVoxels; }
	const btVector3* getVoxelPositionsArray() { return m_voxelPositions; }
};

#endif
