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

#ifndef BT_FLUID_HF_RIGID_DYNAMICS_WORLD_H
#define BT_FLUID_HF_RIGID_DYNAMICS_WORLD_H

#include "btFluidRigidDynamicsWorld.h"

#define DRAWMODE_NORMAL 0
#define DRAWMODE_VELOCITY 1
#define DRAWMODE_MAX 2

#define BODY_DRAWMODE_NORMAL 0
#define BODY_DRAWMODE_VOXEL 1
#define BODY_DRAWMODE_MAX 2

class btFluidHf;
class btFluidHfBuoyantConvexShape;

///World for experimental heightfield fluid simulation
class btFluidHfRigidDynamicsWorld : public btFluidRigidDynamicsWorld
{
	btAlignedObjectArray<btFluidHf*> m_hfFluids;
	int m_drawMode;
	int m_bodyDrawMode;
	
protected:
	virtual void internalSingleStepSimulation(btScalar timeStep);

	void drawFluidHfGround (btIDebugDraw* debugDraw, btFluidHf* fluid);
	void drawFluidHfVelocity (btIDebugDraw* debugDraw, btFluidHf* fluid);
	void drawFluidHfBuoyantConvexShape (btIDebugDraw* debugDrawer, btCollisionObject* object, btFluidHfBuoyantConvexShape* buoyantShape, int voxelDraw);
	void drawFluidHfNormal (btIDebugDraw* debugDraw, btFluidHf* fluid);
	
public:
	btFluidHfRigidDynamicsWorld(btDispatcher* dispatcher, btBroadphaseInterface* pairCache, 
								btConstraintSolver* constraintSolver, btCollisionConfiguration* collisionConfiguration, 
								btFluidSphSolver* sphSolver = 0);
	
	void addFluidHf(btFluidHf* fluid);
	void removeFluidHf(btFluidHf* fluid);
	
	void setDrawMode (int drawMode) { m_drawMode = drawMode; }
	void setBodyDrawMode (int bodyDrawMode) { m_bodyDrawMode = bodyDrawMode; }

	btAlignedObjectArray<btFluidHf*>& getFluidHfArray() { return m_hfFluids; }
	const btAlignedObjectArray<btFluidHf*>& getFluidHfArray() const { return m_hfFluids; }
	
	virtual void debugDrawWorld();
};

#endif //BT_FLUID_HF_RIGID_DYNAMICS_WORLD_H
