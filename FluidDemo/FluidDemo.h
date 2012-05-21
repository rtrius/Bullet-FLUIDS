/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2006 Erwin Coumans  http://continuousphysics.com/Bullet/

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef FLUID_DEMO_H
#define FLUID_DEMO_H

#ifdef _WINDOWS
#include "Win32DemoApplication.h"
#define PlatformDemoApplication Win32DemoApplication
#else
#include "GlutDemoApplication.h"
#define PlatformDemoApplication GlutDemoApplication
#endif

#include "LinearMath/btAlignedObjectArray.h"

//
#include "Fluids/fluid_rendering.h"
#include "Fluids/fluid_system.h"
#include "Fluids/SPHInterface.h"


class btCollisionShape;
class btBroadphaseInterface;
class btCollisionDispatcher;
class btConstraintSolver;
class btDefaultCollisionConfiguration;


///FluidDemo demonstrates Bullet-SPH interactions
class FluidDemo : public PlatformDemoApplication
{
	btAlignedObjectArray<btCollisionShape*>	m_collisionShapes;	//Keep the collision shapes, for deletion/cleanup
	btBroadphaseInterface* m_broadphase;
	btCollisionDispatcher* m_dispatcher;
	btConstraintSolver*	m_solver;
	btDefaultCollisionConfiguration* m_collisionConfiguration;

	FluidsWrapper m_fluids;
	bool m_useMarchingCubes;		//If false, fluid is rendered as spheres
	
public:
	FluidDemo() : m_useMarchingCubes(false) { m_fluids.initialize(); }
	virtual ~FluidDemo() { exitPhysics(); }
	
	void initPhysics();
	void exitPhysics();

	virtual void clientResetScene();
	
	virtual void clientMoveAndDisplay();
	virtual void displayCallback();
	
	virtual void keyboardCallback(unsigned char key, int x, int y);
	
	static DemoApplication* Create()
	{
		FluidDemo* demo = new FluidDemo;
		demo->myinit();
		demo->initPhysics();
		return demo;
	}
};

#endif //FLUID_DEMO_H


