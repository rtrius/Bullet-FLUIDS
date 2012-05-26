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

#include <stdio.h> 	//printf debugging

#include "FluidDemo.h"
#include "GlutStuff.h"

#include "btBulletDynamicsCommon.h"




void FluidDemo::initPhysics()
{
	setTexturing(true);
	setShadows(true);

	setCameraDistance(50.0);

	//
	m_collisionConfiguration = new btDefaultCollisionConfiguration();
	m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);
	m_broadphase = new btDbvtBroadphase();
	m_solver = new btSequentialImpulseConstraintSolver;
	m_dynamicsWorld = new btDiscreteDynamicsWorld(m_dispatcher, m_broadphase, m_solver, m_collisionConfiguration);
	
	m_dynamicsWorld->setGravity( btVector3(0.0, -9.8, 0.0) );

	//Create a very large static box as the ground
	//We can also use DemoApplication::localCreateRigidBody, but for clarity it is provided here:
	{
		btCollisionShape* groundShape = new btBoxShape( btVector3(50.0, 50.0, 50.0) );
		m_collisionShapes.push_back(groundShape);

		btScalar mass(0.f);

		//Rigid bodies are dynamic if and only if mass is non zero, otherwise static
		btVector3 localInertia(0,0,0);
		if(mass != 0.f) groundShape->calculateLocalInertia(mass,localInertia);

		btTransform groundTransform;
		groundTransform.setIdentity();
		groundTransform.setOrigin( btVector3(0,-50,0) );
		
		//Using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
		btDefaultMotionState* myMotionState = new btDefaultMotionState(groundTransform);
		btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, groundShape, localInertia);
		btRigidBody* body = new btRigidBody(rbInfo);

		m_dynamicsWorld->addRigidBody(body);
	}
}
void FluidDemo::exitPhysics()
{
	//Cleanup in the reverse order of creation/initialization

	//Remove the rigidbodies from the dynamics world and delete them
	for(int i = m_dynamicsWorld->getNumCollisionObjects() - 1; i >= 0; i--)
	{
		btCollisionObject* obj = m_dynamicsWorld->getCollisionObjectArray()[i];
		btRigidBody* body = btRigidBody::upcast(obj);
		if( body && body->getMotionState() ) delete body->getMotionState();
		
		m_dynamicsWorld->removeCollisionObject(obj);
		delete obj;
	}

	//Delete collision shapes
	for(int j = 0; j < m_collisionShapes.size(); j++)
	{
		btCollisionShape* shape = m_collisionShapes[j];
		delete shape;
	}
	m_collisionShapes.clear();

	//
	delete m_dynamicsWorld;
	delete m_solver;
	delete m_broadphase;
	delete m_dispatcher;
	delete m_collisionConfiguration;
}

void FluidDemo::initDemos()
{
	m_demos.push_back( new Demo_DamBreak() );
	m_demos.push_back( new Demo_Drop() );
	m_demos.push_back( new Demo_EmitterAndAbsorber() );
	m_demos.push_back( new Demo_DynamicBox() );
	
	for(int i = 0; i < m_demos.size(); ++i) m_demos[i]->initialize(&m_collisionShapes);
	
	m_currentDemoIndex = 0;
	m_maxFluidParticles = 1024;
	startDemo(m_currentDemoIndex);
}
void FluidDemo::exitDemos()
{
	stopDemo(m_currentDemoIndex);

	for(int i = 0; i < m_demos.size(); ++i) 
	{
		if(m_demos[i])
		{
			m_demos[i]->deactivate();
			delete m_demos[i];
		}
	}
	m_demos.clear();
}

void FluidDemo::prevDemo()
{
	if(m_currentDemoIndex - 1 >= 0)
	{
		stopDemo(m_currentDemoIndex);
		startDemo(--m_currentDemoIndex);
	}
}
void FluidDemo::nextDemo()
{
	if(m_currentDemoIndex + 1 <= m_demos.size() - 1)
	{
		stopDemo(m_currentDemoIndex);
		startDemo(++m_currentDemoIndex);
	}
}

void FluidDemo::clientResetScene()
{
	exitPhysics();
	initPhysics();
}

void FluidDemo::clientMoveAndDisplay()
{
	//Simple dynamics world doesn't handle fixed-time-stepping
	float ms = getDeltaTimeMicroseconds();
	float secondsElapsed = ms * 0.000001f;

	if(m_dynamicsWorld)
	{
		m_dynamicsWorld->stepSimulation(secondsElapsed);
		BulletFluidsInterface::stepSimulation(&m_fluids, m_dynamicsWorld, secondsElapsed);
		if( m_demos.size() ) m_demos[m_currentDemoIndex]->update(&m_fluids);
	}	
	
	displayCallback();
}
void FluidDemo::displayCallback(void) 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
	
	renderme();

	
	static bool areSpheresGenerated = false;
	static GLuint glSmallSphereList;
	static GLuint glLargeSphereList;
	if(!areSpheresGenerated)
	{
		glSmallSphereList = generateSphereList(0.1f);
		glLargeSphereList = generateSphereList(0.75f);
	}
	
	switch(m_fluidRenderMode)
	{
		case FRM_Points:
		{
			BT_PROFILE("Draw fluids - points");
			
			for(int i = 0; i < m_fluids.numParticles(); ++i)
			{
				const Fluid &F = const_cast<const FluidSystem&>(m_fluids).getFluid(i);
				drawSphere(glSmallSphereList, F.pos, F.vel.length() * 2.0f);
			}
		}
			break;
			
		case FRM_Spheres:
		{
			BT_PROFILE("Draw fluids - spheres");
			
			for(int i = 0; i < m_fluids.numParticles(); ++i)
			{
				const Fluid &F = const_cast<const FluidSystem&>(m_fluids).getFluid(i);
				drawSphere(glLargeSphereList, F.pos, F.vel.length() * 2.0f);
			}
		}
			break;
		
		case FRM_MarchingCubes:
		{
			BT_PROFILE("Draw fluids - marching cubes");
			
			const int CELLS_PER_EDGE = 32;
			static MarchingCubes *marchingCubes = 0;
			if(!marchingCubes) 
			{
				marchingCubes = new MarchingCubes;
				marchingCubes->initialize(CELLS_PER_EDGE);
			}
			
			marchingCubes->generateMesh(m_fluids);
			
			const std::vector<float> &vertices = marchingCubes->getTriangleVertices();
			if( vertices.size() )
			{
				glEnable(GL_BLEND);
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 

				glEnableClientState(GL_VERTEX_ARRAY);
					glColor4f(0.9f, 0.9f, 0.9f, 0.6f);
					glVertexPointer( 3, GL_FLOAT, 0, &vertices.front() );
					glDrawArrays( GL_TRIANGLES, 0, vertices.size()/3 );	
				glDisableClientState(GL_VERTEX_ARRAY);
			}
		}
			break;
	}
	
	if(m_dynamicsWorld) m_dynamicsWorld->debugDrawWorld();		//Optional but useful: debug drawing to detect problems

	glFlush();
	swapBuffers();
}

void FluidDemo::keyboardCallback(unsigned char key, int x, int y)
{
	switch(key)
	{
		case 'q':
		{
			int currentRenderMode = static_cast<int>(m_fluidRenderMode);
			m_fluidRenderMode = static_cast<FluidRenderMode>(currentRenderMode + 1);
		
			if(m_fluidRenderMode > FRM_MarchingCubes)m_fluidRenderMode = FRM_Points;
			return;
		}
		
		case 'e':
			m_fluids.toggleOpenCL();
			return;
			
		case ' ':
			resetCurrentDemo();
			return;
		
		case '[':
			prevDemo();
			return;
			
		case ']':
			nextDemo();
			return;
			
		case 'n':
			if(m_maxFluidParticles/2 >= MIN_FLUID_PARTICLES) m_maxFluidParticles /= 2;
			resetCurrentDemo();
			return;
		case 'm':
			if(m_maxFluidParticles*2 <= 32768) m_maxFluidParticles *= 2;
			resetCurrentDemo();
			return;
	}
	
	PlatformDemoApplication::keyboardCallback(key, x, y);
}


