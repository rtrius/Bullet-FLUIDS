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
#ifndef FLUID_SPH_HF_DEMO_H
#define FLUID_SPH_HF_DEMO_H

#ifdef _WINDOWS
#include "Win32DemoApplication.h"
#define PlatformDemoApplication Win32DemoApplication
#else
#include "GlutDemoApplication.h"
#define PlatformDemoApplication GlutDemoApplication
#endif

#include "LinearMath/btAlignedObjectArray.h"

//
#include "BulletFluids/Sph/btFluidSphSolver.h"
#include "BulletFluids/btFluidHfRigidCollisionConfiguration.h"
#include "BulletFluids/btFluidHfRigidDynamicsWorld.h"

#include "FluidSphHfInteraction.h"

class btCollisionShape;
class btBroadphaseInterface;
class btCollisionDispatcher;
class btConstraintSolver;
class btDefaultCollisionConfiguration;

class FluidSphHfDemo_GL_ShapeDrawer;

///Highly experimental demonstration of SPH and Heightfield fluid interaction
class FluidSphHfDemo : public PlatformDemoApplication
{
	//Bullet
	btAlignedObjectArray<btCollisionShape*>	m_collisionShapes;	//Keep the collision shapes, for deletion/cleanup
	btBroadphaseInterface* m_broadphase;
	btCollisionDispatcher* m_dispatcher;
	btConstraintSolver*	m_solver;
	btDefaultCollisionConfiguration* m_collisionConfiguration;

	//Fluid system
	btFluidHfRigidDynamicsWorld* m_fluidWorld;
	btFluidSph* m_fluidSph;
	btFluidHf* m_fluidHf;
	
	btFluidSphSolver* m_fluidSphSolver;
	
	FluidSphHfInteractor m_sphHfInteractor;
	
	//Rendering
	FluidSphHfDemo_GL_ShapeDrawer* m_hfFluidShapeDrawer;
	
public:
	FluidSphHfDemo();
	virtual ~FluidSphHfDemo();
	
	void initPhysics();		//Initialize Bullet/fluid system here
	void exitPhysics();		//Deactivate Bullet/fluid system here
	
	//
	virtual void clientMoveAndDisplay();	//Simulation is updated/stepped here
	virtual void displayCallback();			//Rendering occurs here
	
	//
	virtual void keyboardCallback(unsigned char key, int x, int y);
	virtual void specialKeyboard(int key, int x, int y);
	virtual void setShootBoxShape();
	virtual void clientResetScene();
	
	static DemoApplication* Create()
	{
		FluidSphHfDemo* demo = new FluidSphHfDemo;
		demo->myinit();
		return demo;
	}
};

#endif //FLUID_SPH_HF_DEMO_H


