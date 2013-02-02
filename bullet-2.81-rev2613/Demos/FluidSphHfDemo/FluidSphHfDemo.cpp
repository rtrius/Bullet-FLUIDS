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
#include "FluidSphHfDemo.h"

#include <stdio.h> 	//printf debugging

#include "GlutStuff.h"

#include "btBulletDynamicsCommon.h"
#include "LinearMath/btRandom.h"		//GEN_rand(), GEN_RAND_MAX

#include "BulletFluids/Sph/btFluidSph.h"
#include "BulletFluids/Hf/btFluidHf.h"


#include "FluidSphHfDemo_GL_ShapeDrawer.h"

FluidSphHfDemo::FluidSphHfDemo()
{
	m_hfFluidShapeDrawer = new FluidSphHfDemo_GL_ShapeDrawer();
	overrideGLShapeDrawer(m_hfFluidShapeDrawer);

	setTexturing(true);
	setShadows(true);
	setCameraDistance(50.0);
	
	initPhysics();
}
FluidSphHfDemo::~FluidSphHfDemo() 
{
	exitPhysics(); 
}

void FluidSphHfDemo::initPhysics()
{
	//btFluidRigidCollisionConfiguration adds fluid collision algorithms
	m_collisionConfiguration = new btFluidHfRigidCollisionConfiguration();

	//'standard' Bullet configuration
	m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);
	m_broadphase = new btDbvtBroadphase();
	m_solver = new btSequentialImpulseConstraintSolver();
	
	//btFluidSphSolver determines how the particles move and interact with each other
	m_fluidSphSolver = new btFluidSphSolverDefault();						//Standard optimized CPU solver
	
	m_dynamicsWorld = new btFluidHfRigidDynamicsWorld(m_dispatcher, m_broadphase, m_solver, m_collisionConfiguration, m_fluidSphSolver);
	m_fluidWorld = static_cast<btFluidHfRigidDynamicsWorld*>(m_dynamicsWorld);
	
	m_fluidWorld->setGravity( btVector3(0.0, -9.8, 0.0) );

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

		m_fluidWorld->addRigidBody(body);
	}
	
	//Create btFluidSph(s), which contain groups of particles
	{
		const int MAX_PARTICLES = 8192;
		btFluidSph* fluidSph = new btFluidSph(m_fluidWorld->getGlobalParameters(), MAX_PARTICLES);
		
		{
			btFluidSphParametersLocal FL = fluidSph->getLocalParameters();
			
			const btScalar AABB_EXTENT(25.0);
			FL.m_aabbBoundaryMin = btVector3(-AABB_EXTENT, -AABB_EXTENT, -AABB_EXTENT);
			FL.m_aabbBoundaryMax = btVector3(AABB_EXTENT, AABB_EXTENT, AABB_EXTENT);
			FL.m_enableAabbBoundary = 1;
			
			FL.m_stiffness = btScalar(0.5);		//Use lower stiffness to prevent particles from pushing each other away
			
			fluidSph->setLocalParameters(FL);
		}
		
		const bool ENABLE_CCD = true;
		if(ENABLE_CCD) fluidSph->setCcdMotionThreshold( fluidSph->getLocalParameters().m_particleRadius );
		
		m_fluidWorld->addFluidSph(fluidSph);
		m_fluidSph = fluidSph;
	}
	
	
	{
		btFluidHf* fluidHf = new btFluidHf( btScalar(1.0), 50, 50 );
	
		btTransform transform( btQuaternion::getIdentity(), btVector3(-25.0, 1.0, -25.0) );
		fluidHf->setWorldTransform(transform);
		
		const bool START_WITH_FLUID = false;
		if(START_WITH_FLUID)
		{
			for(int z = 0; z < fluidHf->getNumNodesZ(); z++)
				for(int x = 0; x < fluidHf->getNumNodesX(); x++) 
				{
					int index = x + fluidHf->getNumNodesX() * z;
					//if(z < fluidHf->getNumNodesZ()/2)
					//	fluidHf->setFluidHeight( index, btScalar(8.0) );
				}
		}
		
		fluidHf->prep();
		
		m_fluidWorld->addFluidHf(fluidHf);
		m_fluidHf = fluidHf;
	}
	
}
void FluidSphHfDemo::exitPhysics()
{
	//Cleanup in the reverse order of creation/initialization

	//Remove the rigidbodies from the dynamics world and delete them
	for(int i = m_fluidWorld->getNumCollisionObjects() - 1; i >= 0; i--)
	{
		btCollisionObject* obj = m_fluidWorld->getCollisionObjectArray()[i];
		btRigidBody* body = btRigidBody::upcast(obj);
		if( body && body->getMotionState() ) delete body->getMotionState();
		
		m_fluidWorld->removeCollisionObject(obj);	//Removes btCollisionObject, btRigidBody, and btFluidSph
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
	delete m_fluidWorld;
	delete m_solver;
	delete m_broadphase;
	delete m_dispatcher;
	delete m_collisionConfiguration;
	delete m_fluidSphSolver;
}

void FluidSphHfDemo::clientMoveAndDisplay()
{
	btScalar secondsElapsed = getDeltaTimeMicroseconds() * btScalar(0.000001);
	
	const btFluidSphParametersGlobal& FG = m_fluidWorld->getGlobalParameters();
	if(m_fluidWorld)
	{
		const bool USE_SYNCRONIZED_TIME_STEP = false;	//Default: Rigid bodies/HfFluid == 16ms, Sph particles == 3ms time step
		if(USE_SYNCRONIZED_TIME_STEP)
		{
			btScalar timeStep = m_fluidWorld->getGlobalParameters().m_timeStep;
			m_fluidWorld->stepSimulation(timeStep, 0, timeStep);
		}
		else m_fluidWorld->stepSimulation( btScalar(1.0 / 60.0), 0 );
		//else m_fluidWorld->stepSimulation( secondsElapsed, 10, btScalar(1.0 / 60.0) );
		
		
		if(m_fluidSph && m_fluidHf)
		{
			btFluidSph* fluidSph = m_fluidSph;
			const btFluidSphParametersLocal& FL = fluidSph->getLocalParameters();
			
			btFluidHf* fluidHf = m_fluidHf;
			
			for(int n = 0; n < fluidSph->numParticles(); ++n)
			{
				const btVector3& transformedParticlePos = fluidHf->getWorldTransform().invXform( fluidSph->getPosition(n) );
				
				int columnIndex = fluidHf->arrayIndex( transformedParticlePos.x(), transformedParticlePos.z() );
				if( 0 <= columnIndex && columnIndex < fluidHf->getNumNodes() )	//AABB test
				{
					const btScalar MARGIN(0.05);
					if( transformedParticlePos.y() - FL.m_particleRadius - MARGIN < fluidHf->getCombinedHeight(columnIndex) )
					{
						const btScalar PARTICLE_HEIGHT_CONTRIBUTION(3.0);
					
						//fluidSph->applyForce(n, btVector3(0, 0.001, 0) );
					
						fluidSph->markParticleForRemoval(n);
						fluidHf->addFluidHeight(columnIndex, PARTICLE_HEIGHT_CONTRIBUTION);
					}
				}
				
			}
		}
	}
	
	{
		static int counter = 0;
		if(++counter > 100)
		{
			counter = 0;
			
			for(int i = 0; i < m_fluidWorld->getNumFluids(); ++i)
				printf( "m_fluidWorld->getFluid(%d)->numParticles(): %d \n", i, m_fluidWorld->getFluid(i)->numParticles() );
		}
	}
		
	displayCallback();
}

void FluidSphHfDemo::displayCallback(void) 
{
	//BT_PROFILE() does not work correctly in this function;
	//timings are captured only when the camera is moving.
	//BT_PROFILE("FluidSphHfDemo::displayCallback()");

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

	renderme();
	
	if(m_fluidWorld) m_fluidWorld->debugDrawWorld();		//Optional but useful: debug drawing to detect problems

	glFlush();
	swapBuffers();
}


inline int emitParticle(btFluidSph* fluid, const btVector3& position, const btVector3& velocity)
{
	int index = fluid->addParticle(position);
	if( index != fluid->numParticles() ) fluid->setVelocity(index, velocity);
	else
	{
		index = ( fluid->numParticles() - 1 ) * GEN_rand() / GEN_RAND_MAX;		//Random index
	
		fluid->setPosition(index, position);
		fluid->setVelocity(index, velocity);
	}	
	
	return index;
}
void FluidSphHfDemo::keyboardCallback(unsigned char key, int x, int y)
{
	switch(key)
	{
		case '/':
			if(m_fluidSph)
			{
				const btScalar SPEED(1.0);
				btVector3 position = getCameraPosition();
				position.setY( position.y() - btScalar(5.0) );
				
				btVector3 direction = (getRayTo(x,y) - position).normalized();
				btVector3 velocity = direction * SPEED;
				
				const btVector3 X_AXIS(1, 0, 0);
				const btVector3 Y_AXIS(0, 1, 0);
				const btVector3 Z_AXIS(0, 0, 1);
				btQuaternion rotation = shortestArcQuat(Z_AXIS, direction);
				btVector3 up = quatRotate(rotation, Y_AXIS);
				btVector3 left = quatRotate(rotation, X_AXIS);
				
				const btScalar SPACING(2.5);
				emitParticle(m_fluidSph, position, velocity);
				emitParticle(m_fluidSph, position + up*SPACING, velocity);
				emitParticle(m_fluidSph, position + -up*SPACING, velocity);
				emitParticle(m_fluidSph, position + left*SPACING, velocity);
				emitParticle(m_fluidSph, position + -left*SPACING, velocity);
			}
			break;
			
		case 'k':
			m_hfFluidShapeDrawer->m_drawHfFluidWithTriangles = !m_hfFluidShapeDrawer->m_drawHfFluidWithTriangles;
			break;
		case 'l':
			m_hfFluidShapeDrawer->m_drawHfGroundWithTriangles = !m_hfFluidShapeDrawer->m_drawHfGroundWithTriangles;
			break;
		case ';':
			m_hfFluidShapeDrawer->m_drawHfFluidAsColumns = !m_hfFluidShapeDrawer->m_drawHfFluidAsColumns;
			break;
		case '\'':
			m_hfFluidShapeDrawer->m_drawHfGroundAsColumns = !m_hfFluidShapeDrawer->m_drawHfGroundAsColumns;
			break;
			
		default:
			PlatformDemoApplication::keyboardCallback(key, x, y);
			break;
	}
}

void FluidSphHfDemo::specialKeyboard(int key, int x, int y)
{
	switch(key)
	{
		case GLUT_KEY_END:
		{
			int numObj = getDynamicsWorld()->getNumCollisionObjects();
			if(numObj)
			{
				btCollisionObject* obj = getDynamicsWorld()->getCollisionObjectArray()[numObj-1];
				if( btFluidSph::upcast(obj) || btFluidHf::upcast(obj) ) return;	//Deleting btFluidSph/btFluidHf will cause crashes in FluidSphHfDemo
				
				getDynamicsWorld()->removeCollisionObject(obj);
				
				btRigidBody* body = btRigidBody::upcast(obj);
				if(body && body->getMotionState()) delete body->getMotionState();
				delete obj;

			}
			return;
		}
	}
	
	PlatformDemoApplication::specialKeyboard(key, x, y);
}

void FluidSphHfDemo::setShootBoxShape()
{
	if (!m_shootBoxShape)
	{
		const btScalar BOX_DIMENSIONS(3.0);
	
		btBoxShape* box = new btBoxShape( btVector3(BOX_DIMENSIONS, BOX_DIMENSIONS, BOX_DIMENSIONS) );
		box->initializePolyhedralFeatures();
		m_shootBoxShape = box;
	}
}

void FluidSphHfDemo::clientResetScene()
{
	exitPhysics();
	initPhysics();
}

