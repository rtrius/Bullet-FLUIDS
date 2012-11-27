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
#ifndef HF_FLUID_DEMO_H
#define HF_FLUID_DEMO_H

#ifdef _WINDOWS
#include "Win32DemoApplication.h"
#define PlatformDemoApplication Win32DemoApplication
#else
#include "GlutDemoApplication.h"
#define PlatformDemoApplication GlutDemoApplication
#endif

#include "LinearMath/btAlignedObjectArray.h"
#include "BulletFluids/Hf/btFluidHf.h"

class btBroadphaseInterface;
class btCollisionShape;
class btCollisionDispatcher;
class btConstraintSolver;
class btDefaultCollisionConfiguration;

class btFluidHfRigidDynamicsWorld;

class FluidHfDemo_GL_ShapeDrawer;

///FluidHfDemo demonstrates buoyancy / Heightfield fluids
class FluidHfDemo : public PlatformDemoApplication
{
	FluidHfDemo_GL_ShapeDrawer* m_hfFluidShapeDrawer;

public:

	//keep the collision shapes, for deletion/cleanup
	btAlignedObjectArray<btCollisionShape*>	m_collisionShapes;

	btBroadphaseInterface* m_broadphase;
	btCollisionDispatcher* m_dispatcher;
	btConstraintSolver*	m_solver;
	
	btDefaultCollisionConfiguration* m_collisionConfiguration;

public:
	FluidHfDemo ();
	virtual ~FluidHfDemo() { exitPhysics(); }

	void initPhysics();
	void exitPhysics();
	
	virtual void clientMoveAndDisplay();
	virtual void displayCallback();
	void renderme();
	
	//
	void clientResetScene();
	void keyboardCallback(unsigned char key, int x, int y);

	virtual void setShootBoxShape();
	
	virtual const btFluidHfRigidDynamicsWorld* getFluidHfDynamicsWorld() const
	{
		///just make it a btFluidHfRigidDynamicsWorld please or we will add type checking
		return (btFluidHfRigidDynamicsWorld*) m_dynamicsWorld;
	}
	virtual btFluidHfRigidDynamicsWorld* getFluidHfDynamicsWorld()
	{
		///just make it a btFluidHfRigidDynamicsWorld please or we will add type checking
		return (btFluidHfRigidDynamicsWorld*) m_dynamicsWorld;
	}
	
	static DemoApplication* Create()
	{
		FluidHfDemo* demo = new FluidHfDemo;
		demo->myinit();
		demo->initPhysics();
		return demo;
	}
};

#endif //HF_FLUID_DEMO_H





