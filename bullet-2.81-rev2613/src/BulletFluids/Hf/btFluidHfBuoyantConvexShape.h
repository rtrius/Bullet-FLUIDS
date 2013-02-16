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


///Wrapper for btConvexShape that allows a btRigidBody to interact with a btFluidHf by adding a voxelized representation.
///@remarks
///Make sure to call generateShape() after creating this class.
class btFluidHfBuoyantConvexShape : public btCollisionShape
{
protected:
	bool m_useFluidDensity;
	btScalar m_buoyancyScale;
	
	btScalar m_collisionRadius;
	btScalar m_totalVolume;
	btScalar m_volumePerVoxel;
	
	btAlignedObjectArray<btVector3> m_voxelPositions;
	
	btConvexShape* m_convexShape;
	
public:
	btFluidHfBuoyantConvexShape(btConvexShape* convexShape);
	virtual ~btFluidHfBuoyantConvexShape() {}
	
	///This function must be called before the shape may interact with a btFluidHf; it voxelizes the shape.
	///@param collisionRadius Radius of a sphere/voxel when colliding; it should not be much larger than the width of a fluid column.
	///@param volumeEstimationRadius Radius of a voxel used for estimating the volume of the shape; 
	///this should be equal to or smaller than collisionRadius, and is not used outside of generateShape().
	void generateShape(btScalar collisionRadius, btScalar volumeEstimationRadius);

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

	btScalar getVoxelRadius() const { return m_collisionRadius; }
	btScalar getTotalVolume() const { return m_totalVolume; }
	btScalar getVolumePerVoxel() const { return m_volumePerVoxel; }
	int getNumVoxels() const { return m_voxelPositions.size(); }
	const btVector3& getVoxelPosition(int index) { return m_voxelPositions[index]; }
	
	
	///If true, the rigid body is treated as having the same density as the fluid.
	///This makes it easier to control buoyancy; when enabled a setBuoyancyScale() of 1/2 means that the
	///rigid is twice as dense compared to the fluid, 1/4 is 4 times as dense, and 2 is half as dense. 
	///False by default.
	void useFluidDensity(bool useFluidDensity) { m_useFluidDensity = useFluidDensity; }
	bool isUsingFluidDensity() const { return m_useFluidDensity; }
	
	///Scales the magnitude of the buoyancy force; default 1.0.
	btScalar getBuoyancyScale() const { return m_buoyancyScale; }
	void setBuoyancyScale(btScalar scale) { m_buoyancyScale = scale; }	
};

#endif
