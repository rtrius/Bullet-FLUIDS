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
		
		//
		const bool SPAWN_WALLS = false;
		if(SPAWN_WALLS)
		{
			btCollisionShape* wallShape = new btBoxShape( btVector3(5.0, 15.0, 50.0) );
			m_collisionShapes.push_back(wallShape);
			
			btTransform transformA( btQuaternion::getIdentity(), btVector3(25.0, 10.0, 0.0) );
			btDefaultMotionState* motionStateA = new btDefaultMotionState(transformA);
			btRigidBody::btRigidBodyConstructionInfo rbInfoA(mass, motionStateA, wallShape, localInertia);
			btRigidBody* bodyA = new btRigidBody(rbInfoA);
			m_dynamicsWorld->addRigidBody(bodyA);
			
			btTransform transformB( btQuaternion::getIdentity(), btVector3(-25.0, 10.0, 0.0) );
			btDefaultMotionState* motionStateB = new btDefaultMotionState(transformB);
			btRigidBody::btRigidBodyConstructionInfo rbInfoB(mass, motionStateB, wallShape, localInertia);
			btRigidBody* bodyB = new btRigidBody(rbInfoB);
			m_dynamicsWorld->addRigidBody(bodyB);
			
			btTransform transformC( btQuaternion(3.141*0.5, 0.0, 0.0), btVector3(0.0, 10.0, 25.0) );
			btDefaultMotionState* motionStateC = new btDefaultMotionState(transformC);
			btRigidBody::btRigidBodyConstructionInfo rbInfoC(mass, motionStateC, wallShape, localInertia);
			btRigidBody* bodyC = new btRigidBody(rbInfoC);
			m_dynamicsWorld->addRigidBody(bodyC);
			
			btTransform transformD( btQuaternion(3.141*0.5, 0.0, 0.0), btVector3(55.0, 10.0, -25.0) );
			btDefaultMotionState* motionStateD = new btDefaultMotionState(transformD);
			btRigidBody::btRigidBodyConstructionInfo rbInfoD(mass, motionStateD, wallShape, localInertia);
			btRigidBody* bodyD = new btRigidBody(rbInfoD);
			m_dynamicsWorld->addRigidBody(bodyD);
			
			btTransform transformE( btQuaternion(3.141*0.5, 0.0, 0.0), btVector3(-55.0, 10.0, -25.0) );
			btDefaultMotionState* motionStateE = new btDefaultMotionState(transformE);
			btRigidBody::btRigidBodyConstructionInfo rbInfoE(mass, motionStateE, wallShape, localInertia);
			btRigidBody* bodyE = new btRigidBody(rbInfoE);
			m_dynamicsWorld->addRigidBody(bodyE);
		}
	}


	//Create a large dynamic box
	const bool SPAWN_LARGE_DYNAMIC_BOX = false;
	if(SPAWN_LARGE_DYNAMIC_BOX)
	{
		btCollisionShape* colShape = new btBoxShape( btVector3(10.0, 10.0, 10.0) );
		m_collisionShapes.push_back(colShape);

		btScalar mass(1.f);

		//Rigid bodies are dynamic if and only if mass is non zero, otherwise static
		btVector3 localInertia(0,0,0);
		if(mass != 0.f) colShape->calculateLocalInertia(mass, localInertia);
			
		float start_x = -5.0 - 1/2;
		float start_y = -5.0;
		float start_z = -3.0 - 1/2;
		
		btTransform startTransform;
		startTransform.setIdentity();
		startTransform.setOrigin( btVector3(start_x, start_y + 20.0, start_z) );

		//Using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
		btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
		btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, myMotionState, colShape, localInertia);
		btRigidBody* body = new btRigidBody(rbInfo);

		m_dynamicsWorld->addRigidBody(body);
	}

	//	next: add a bucket composed of constrained btBoxShapes
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
		m_fluids.stepSimulation(m_dynamicsWorld, secondsElapsed);
	}	
	
	displayCallback();
}
void FluidDemo::displayCallback(void) 
{
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
	
	renderme();

	if(m_useMarchingCubes) 
	{
		BT_PROFILE("Draw fluids - marching cubes");
		
		const int CELLS_PER_EDGE = 32;
		static MarchingCubes *marchingCubes = 0;
		if(!marchingCubes) 
		{
			marchingCubes = new MarchingCubes;
			marchingCubes->initialize(CELLS_PER_EDGE);
		}
		
		marchingCubes->generateMesh(m_fluids.m_fluidSystem);
		
		const std::vector<float> &vertices = marchingCubes->getTriangleVertices();
	
		glEnable(GL_BLEND);
		glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 

		glEnableClientState(GL_VERTEX_ARRAY);
			glColor4f(0.9f, 0.9f, 0.9f, 0.6f);
			glVertexPointer( 3, GL_FLOAT, 0, &vertices.front() );
			glDrawArrays( GL_TRIANGLES, 0, vertices.size()/3 );		
		glDisableClientState(GL_VERTEX_ARRAY);
	}
	else 
	{
		BT_PROFILE("Draw fluids - spheres");
		
		for(int i = 0; i < m_fluids.m_fluidSystem.numParticles(); ++i)
		{
			const Fluid &F = const_cast<const FluidSystem&>(m_fluids.m_fluidSystem).getFluid(i);
			drawSphere(F.pos, F.vel.length() * 2.0f);
		}
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
			m_useMarchingCubes = !m_useMarchingCubes;
			return;
			
		case 'e':
			m_fluids.toggleOpenCL();
			return;
	}
	
	PlatformDemoApplication::keyboardCallback(key, x, y);
}





