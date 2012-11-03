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

//Contains '#include <GL/glew.h>', which must be included before <GL/gl.h>
#include "FluidRendering/ScreenSpaceFluidRendererGL.h"	

#ifdef _WINDOWS
#include "Win32DemoApplication.h"
#define PlatformDemoApplication Win32DemoApplication
#else
#include "GlutDemoApplication.h"
#define PlatformDemoApplication GlutDemoApplication
#endif

#include "LinearMath/btAlignedObjectArray.h"

//
#include "FluidRendering/MarchingCubes.h"
#include "Fluids/btFluidSph.h"
#include "Fluids/btFluidSolver.h"
#include "Fluids/btFluidSolverMultiphase.h"
#include "Fluids/btFluidRigidCollisionConfiguration.h"


#include "demos.h"

class btCollisionShape;
class btBroadphaseInterface;
class btCollisionDispatcher;
class btConstraintSolver;
class btDefaultCollisionConfiguration;


const int MIN_FLUID_PARTICLES = 64;
const int MAX_FLUID_PARTICLES = 65536;


enum FluidRenderMode
{
	FRM_Points = 0,
	FRM_MediumSpheres,
	FRM_LargeSpheres,
	FRM_ScreenSpace,
	FRM_MarchingCubes
};


#define ENABLE_OPENCL_FLUID_SOLVER
#ifdef ENABLE_OPENCL_FLUID_SOLVER
	#include "Fluids/OpenCL_support/btExperimentsOpenCL/btOpenCLUtils.h"
	#include "Fluids/OpenCL_support/btFluidSolverOpenCL.h"
	class OpenCLConfig
	{
		cl_platform_id m_platformId;
		
	public:	
		cl_device_id m_device;
		
		cl_context m_context;
		cl_command_queue m_commandQueue;
		
		OpenCLConfig()
		{
			cl_int error;
			
			m_context = btOpenCLUtils::createContextFromType(CL_DEVICE_TYPE_GPU, &error, 0, 0, -1, -1, &m_platformId);
			oclCHECKERROR(error, CL_SUCCESS);
			
			if(m_context && m_platformId) btOpenCLUtils::printPlatformInfo(m_platformId);
			
			if( btOpenCLUtils::getNumDevices(m_context) > 0 )
			{
				m_device = btOpenCLUtils::getDevice(m_context, 0);
				if(m_device) btOpenCLUtils::printDeviceInfo(m_device);
				
				m_commandQueue = clCreateCommandQueue(m_context, m_device, 0, &error);
				oclCHECKERROR(error, CL_SUCCESS);
			}
		}
		~OpenCLConfig()
		{
			clReleaseCommandQueue(m_commandQueue);
			clReleaseContext(m_context);
		}
	};
#endif

///FluidDemo demonstrates Bullet-SPH interactions
class FluidDemo : public PlatformDemoApplication
{
	//Bullet
	btAlignedObjectArray<btCollisionShape*>	m_collisionShapes;	//Keep the collision shapes, for deletion/cleanup
	btBroadphaseInterface* m_broadphase;
	btCollisionDispatcher* m_dispatcher;
	btConstraintSolver*	m_solver;
	btDefaultCollisionConfiguration* m_collisionConfiguration;

	//Fluid system
	btFluidRigidDynamicsWorld* m_fluidWorld;
	btAlignedObjectArray<btFluidSph*> m_fluids;
			
	bool m_useFluidSolverOpenCL;
	btFluidSolver* m_fluidSolverCPU;
	btFluidSolver* m_fluidSolverGPU;
	
	//Rendering
	FluidRenderMode m_fluidRenderMode;
	ScreenSpaceFluidRendererGL* m_screenSpaceRenderer;
	
	//Demos
	btAlignedObjectArray<FluidSystemDemo*> m_demos;
	int m_currentDemoIndex;
	int m_maxFluidParticles;
	
public:
	FluidDemo();
	virtual ~FluidDemo();
	
	virtual void clientMoveAndDisplay();	//Simulation is updated/stepped here
	virtual void displayCallback();			//Rendering occurs here
	
	void renderFluids();
	
	//
	void initPhysics();		//Initialize Bullet and Fluid System
	void exitPhysics();		//Deactivate Bullet and Fluid System

	//
	void initDemos();
	void exitDemos();
	
	void startDemo(int index)
	{
		m_demos[index]->addToWorld(m_fluidWorld);
		m_demos[index]->reset(*m_fluidWorld, &m_fluids, m_maxFluidParticles, false);
	}
	void stopDemo(int index) { m_demos[index]->removeFromWorld(m_fluidWorld); }
	void resetCurrentDemo()
	{
		printf("m_maxFluidParticles: %d\n", m_maxFluidParticles);
		m_demos[m_currentDemoIndex]->reset(*m_fluidWorld, &m_fluids, m_maxFluidParticles, true);
	}
	
	void prevDemo();
	void nextDemo();
	
	//
	virtual void keyboardCallback(unsigned char key, int x, int y);
	virtual void setShootBoxShape();
	virtual void myinit();
	virtual void reshape(int w, int h);
	virtual void clientResetScene();
	
	static DemoApplication* Create()
	{
		FluidDemo* demo = new FluidDemo;
		demo->myinit();
		return demo;
	}
};

#endif //FLUID_DEMO_H


