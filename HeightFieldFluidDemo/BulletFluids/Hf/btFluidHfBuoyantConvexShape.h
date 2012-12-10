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

#ifndef BT_FLUID_HF_BUOYANT_CONVEX_SHAPE_H
#define BT_FLUID_HF_BUOYANT_CONVEX_SHAPE_H

#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btVector3.h"
#include "BulletCollision/CollisionShapes/btCollisionShape.h"
#include "BulletCollision/CollisionShapes/btConvexShape.h"


///Wrapper for btConvexShape that allows it to interact with a btFluidHf by adding a voxelized representation.
class btFluidHfBuoyantConvexShape : public btCollisionShape
{
protected:
	btScalar m_floatyness;
	btScalar m_radius;
	btScalar m_totalVolume;
	btScalar m_volumePerVoxel;
	
	btAlignedObjectArray<btVector3> m_voxelPositions;
	
	btConvexShape* m_convexShape;
	
public:
	btFluidHfBuoyantConvexShape(btConvexShape* convexShape);
	virtual ~btFluidHfBuoyantConvexShape() {}
	
	///This function must be called before the shape may interact with a btFluidHf; it voxelizes the shape.
	///@param radius Radius of a sphere/voxel; it should not be much larger than the width of a fluid column.
	///@param gap Spacing between each voxel.
	void generateShape(btScalar radius, btScalar gap);

	btConvexShape* getConvexShape() { return m_convexShape; }
	const btConvexShape* getConvexShape() const { return m_convexShape; }

	virtual void getAabb(const btTransform& t,btVector3& aabbMin,btVector3& aabbMax) const 
	{ 
		return m_convexShape->getAabb(t, aabbMin, aabbMax);
	}
	virtual void setMargin(btScalar margin) { m_convexShape->setMargin (margin); }
	virtual btScalar getMargin() const { return m_convexShape->getMargin(); }
	virtual void setLocalScaling(const btVector3& scaling) { m_convexShape->setLocalScaling (scaling); }
	virtual const btVector3& getLocalScaling() const { return m_convexShape->getLocalScaling(); }
	virtual void calculateLocalInertia(btScalar mass, btVector3& inertia) const
	{ 
		m_convexShape->calculateLocalInertia (mass, inertia); 
	}
	virtual const char*	getName() const { return "FLUID_HF_BUOYANT_CONVEX_SHAPE"; }

	btScalar getVoxelRadius() const { return m_radius; }
	btScalar getTotalVolume() const { return m_totalVolume; }
	btScalar getVolumePerVoxel() const { return m_volumePerVoxel; }
	int getNumVoxels() const { return m_voxelPositions.size(); }
	const btVector3& getVoxelPosition(int index) { return m_voxelPositions[index]; }
	
	///Scales the buoyancy force.
	btScalar getFloatyness() const { return m_floatyness; }
	void setFloatyness(btScalar floatyness) { m_floatyness = floatyness; }
};

#endif
