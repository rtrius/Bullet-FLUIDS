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
#include "Fluids/FluidSph.h"
#include "Fluids/SPHInterface.h"
#include "Fluids/FluidSolver.h"

#include "demos.h"

class btCollisionShape;
class btBroadphaseInterface;
class btCollisionDispatcher;
class btConstraintSolver;
class btDefaultCollisionConfiguration;


const int MIN_FLUID_PARTICLES = 64;


enum FluidRenderMode
{
	FRM_Points = 0,
	FRM_Spheres = 1,
	FRM_MarchingCubes = 2
};

#define ENABLE_OPENCL_FLUID_SOLVER
#ifdef ENABLE_OPENCL_FLUID_SOLVER
	#include "Fluids/OpenCL_support/FluidSolverOpenCL.h"
#endif

///FluidDemo demonstrates Bullet-SPH interactions
class FluidDemo : public PlatformDemoApplication
{
	btAlignedObjectArray<btCollisionShape*>	m_collisionShapes;	//Keep the collision shapes, for deletion/cleanup
	btBroadphaseInterface* m_broadphase;
	btCollisionDispatcher* m_dispatcher;
	btConstraintSolver*	m_solver;
	btDefaultCollisionConfiguration* m_collisionConfiguration;

	FluidWorld *m_fluidWorld;
	FluidSph m_fluid;
	FluidRenderMode m_fluidRenderMode;
	
	btAlignedObjectArray<FluidSystemDemo*> m_demos;
	int m_currentDemoIndex;
	int m_maxFluidParticles;
	
	bool m_useFluidSolverOpenCL;
	FluidSolver *m_fluidSolverCPU;
	FluidSolver *m_fluidSolverGPU;
	
public:
	FluidDemo() : m_fluidRenderMode(FRM_Points), m_maxFluidParticles(MIN_FLUID_PARTICLES), m_useFluidSolverOpenCL(false)
	{	
		m_fluidWorld = 0;
		m_fluidSolverCPU = 0;
		m_fluidSolverGPU = 0;
		
		//m_fluidSolverCPU = new FluidSolverGridNeighbor();
		m_fluidSolverCPU = new FluidSolverReducedGridNeighbor();
		
#ifdef ENABLE_OPENCL_FLUID_SOLVER
		m_fluidSolverGPU = new FluidSolverOpenCL();
#endif

		m_fluidWorld = new FluidWorld(m_fluidSolverCPU);
		m_fluidWorld->addFluid(&m_fluid);
	
		initDemos();
	}
	virtual ~FluidDemo() 
	{
		if(m_fluidWorld) delete m_fluidWorld;
		if(m_fluidSolverCPU)delete m_fluidSolverCPU;
		if(m_fluidSolverGPU)delete m_fluidSolverGPU;
	
		exitDemos();
		exitPhysics(); 
	}
	
	void initPhysics();
	void exitPhysics();

	void initDemos();
	void exitDemos();
	
	void startDemo(int index)
	{
		m_demos[index]->addToWorld(m_dynamicsWorld);
		m_demos[index]->reset(*m_fluidWorld, &m_fluid, m_maxFluidParticles);
	}
	void stopDemo(int index) { m_demos[index]->removeFromWorld(m_dynamicsWorld); }
	void resetCurrentDemo()
	{
		printf("m_maxFluidParticles: %d\n", m_maxFluidParticles);
		m_demos[m_currentDemoIndex]->reset(*m_fluidWorld, &m_fluid, m_maxFluidParticles);
	}
	
	void prevDemo();
	void nextDemo();
	
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


