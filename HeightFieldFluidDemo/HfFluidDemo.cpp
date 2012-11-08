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
#include "HfFluidDemo.h"

#include <stdio.h> //printf debugging

#include "btBulletDynamicsCommon.h"
#include "LinearMath/btQuickprof.h"

#include "BulletHfFluid/btHfFluidRigidDynamicsWorld.h"
#include "BulletHfFluid/btHfFluid.h"
#include "BulletHfFluid/btHfFluidRigidCollisionConfiguration.h"
#include "BulletHfFluid/btHfFluidBuoyantConvexShape.h"

#include "hfDemos.h"
#include "HfFluidDemo_GL_ShapeDrawer.h"


unsigned int current_draw_mode = DRAWMODE_NORMAL;
unsigned int current_body_draw_mode = 0;
unsigned	current_demo = 0;

HfFluidDemo::HfFluidDemo()
{
	overrideGLShapeDrawer( new HfFluidDemo_GL_ShapeDrawer() );
	setTexturing(true);
	setShadows(true);
}
void HfFluidDemo::initPhysics()
{
	m_collisionConfiguration = new btHfFluidRigidCollisionConfiguration();
	m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);
	m_broadphase = new btDbvtBroadphase();
	m_solver = new btSequentialImpulseConstraintSolver();
	
	m_dynamicsWorld = new btHfFluidRigidDynamicsWorld(m_dispatcher, m_broadphase, m_solver, m_collisionConfiguration);
	m_dynamicsWorld->setGravity( btVector3(0, -10, 0) );
	
	//
	const btScalar CUBE_HALF_EXTENTS = btScalar(1.5);
	btCollisionShape* groundShape = new btBoxShape( btVector3(200, CUBE_HALF_EXTENTS, 200) );
	{
		m_collisionShapes.push_back(groundShape);
	
		btTransform transform( btQuaternion::getIdentity(), btVector3(0, -12, 0) );
		localCreateRigidBody(0.f, transform, groundShape);
	}
	
	//
	clientResetScene();
}
void HfFluidDemo::exitPhysics()
{
	//cleanup in the reverse order of creation/initialization

	//remove the rigidbodies from the dynamics world and delete them
	for(int i=m_dynamicsWorld->getNumCollisionObjects()-1; i>=0 ;i--)
	{
		btCollisionObject* obj = m_dynamicsWorld->getCollisionObjectArray()[i];
		btRigidBody* body = btRigidBody::upcast(obj);
		if (body && body->getMotionState())
		{
			delete body->getMotionState();
		}
		m_dynamicsWorld->removeCollisionObject( obj );
		delete obj;
	}

	//delete collision shapes
	for (int j=0;j<m_collisionShapes.size();j++)
	{
		btCollisionShape* shape = m_collisionShapes[j];
		m_collisionShapes[j] = 0;
		delete shape;
	}

	delete m_dynamicsWorld;

	delete m_solver;
	delete m_broadphase;
	delete m_dispatcher;
	delete m_collisionConfiguration;
}

void HfFluidDemo::clientMoveAndDisplay()
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
	
	if(m_dynamicsWorld)
	{	
		if(demo_run_functions[current_demo]) demo_run_functions[current_demo]( getHfFluidDynamicsWorld() );
		
		const btScalar DELTA_TIME = btScalar(1.0 / 60.0);	
		m_dynamicsWorld->stepSimulation(DELTA_TIME, 0);
	}
	
	renderme();  //render the graphics objects, with center of mass shift

	updateCamera();

	glFlush();
	glutSwapBuffers();
}
void HfFluidDemo::displayCallback(void) 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

	renderme();

	glFlush();
	glutSwapBuffers();
}
void HfFluidDemo::renderme()
{
	m_dynamicsWorld->debugDrawWorld();
	
	DemoApplication::renderme();	
}



void HfFluidDemo::clientResetScene()
{
	DemoApplication::clientResetScene();
	
	//Clean up	 
	for(int i = m_dynamicsWorld->getNumCollisionObjects() - 1; i > 0; i--)
	{
		btCollisionObject* obj = m_dynamicsWorld->getCollisionObjectArray()[i];
		
		btRigidBody* body = btRigidBody::upcast(obj);
		if( body && body->getMotionState() ) delete body->getMotionState();
		
		while(m_dynamicsWorld->getNumConstraints())
		{
			btTypedConstraint* pc = m_dynamicsWorld->getConstraint(0);
			m_dynamicsWorld->removeConstraint(pc);
			delete pc;
		}
		
		btHfFluid* hfFluid = btHfFluid::upcast(obj);
		if(hfFluid) getHfFluidDynamicsWorld()->removeHfFluid(hfFluid);
		else m_dynamicsWorld->removeCollisionObject(obj);
		
		delete obj;
	}
	
	//Init
	printf("current_demo = %d\n", current_demo);
	m_azi = g_azi_array[current_demo];
	m_ele = g_ele_array[current_demo];
	m_cameraDistance = g_cameraDistance_array[current_demo];
	updateCamera();
	demo_init_functions[current_demo]( getHfFluidDynamicsWorld(), m_collisionShapes );
}



void HfFluidDemo::keyboardCallback(unsigned char key, int x, int y)
{
	switch(key)
	{
		case ']':
			current_demo = (current_demo+1)%NUM_DEMOS;
			clientResetScene();
			break;
		case '[':
			current_demo = (current_demo-1)%NUM_DEMOS;
			clientResetScene();
			break;
		case '.':
			current_draw_mode = (current_draw_mode+1) % DRAWMODE_MAX;
			getHfFluidDynamicsWorld()->setDrawMode (current_draw_mode);
			break;
		case 'v':
			current_body_draw_mode = (current_body_draw_mode+1) % BODY_DRAWMODE_MAX;
			getHfFluidDynamicsWorld()->setBodyDrawMode (current_body_draw_mode);
			break;
		default:
			DemoApplication::keyboardCallback(key,x,y);
			break;
	}
}

void HfFluidDemo::setShootBoxShape ()
{
	if (!m_shootBoxShape)
	{
		m_shootBoxShape = new btBoxShape(btVector3(0.3f,1.f,0.2f));
		btHfFluidBuoyantConvexShape* buoyantShape = new btHfFluidBuoyantConvexShape((btConvexShape*)m_shootBoxShape);
		buoyantShape->generateShape (btScalar(0.25f), btScalar(0.05f));
		m_shootBoxShape = buoyantShape;
	}
}

