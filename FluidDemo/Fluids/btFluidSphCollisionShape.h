/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. 
   If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_FLUID_SPH_COLLISION_SHAPE_H
#define BT_FLUID_SPH_COLLISION_SHAPE_H

#include "BulletCollision/CollisionShapes/btConcaveShape.h"

///Internal class needed for SPH particle collision detection. Each btFluidSph has a unique btFluidSphCollisionShape.
class btFluidSphCollisionShape : public btConcaveShape
{
public:
	btFluidSph*	m_owner;

	btFluidSphCollisionShape(btFluidSph* owner)
	{
		m_shapeType = CUSTOM_CONCAVE_SHAPE_TYPE; 
		m_owner = owner; 
	}

	void processAllTriangles(btTriangleCallback* callback, const btVector3& aabbMin, const btVector3& aabbMax) const {}
	
	virtual void getAabb(const btTransform& t, btVector3& aabbMin, btVector3& aabbMax) const
	{
		// t is usually identity, except when colliding against btCompoundShape. See Issue 512
		btVector3 mins, maxs;
		m_owner->getAabb(mins, maxs);
		
		const btVector3	corners[] = {	t * btVector3(mins.x(),mins.y(),mins.z()),
										t * btVector3(maxs.x(),mins.y(),mins.z()),
										t * btVector3(mins.x(),maxs.y(),mins.z()),
										t * btVector3(maxs.x(),maxs.y(),mins.z()),
										t * btVector3(mins.x(),mins.y(),maxs.z()),
										t * btVector3(maxs.x(),mins.y(),maxs.z()),
										t * btVector3(mins.x(),maxs.y(),maxs.z()),
										t * btVector3(maxs.x(),maxs.y(),maxs.z()) };
		aabbMin = aabbMax = corners[0];
		for(int i = 1; i < 8;++i)
		{
			aabbMin.setMin(corners[i]);
			aabbMax.setMax(corners[i]);
		}
	}
	
	virtual void setLocalScaling(const btVector3& scaling) {}
	virtual const btVector3& getLocalScaling() const
	{
		static const btVector3 dummy(1,1,1);
		return dummy;
	}
	virtual void calculateLocalInertia(btScalar mass, btVector3& inertia) const { btAssert(0); }
	virtual const char* getName() const { return "FluidSph"; }
};

#endif