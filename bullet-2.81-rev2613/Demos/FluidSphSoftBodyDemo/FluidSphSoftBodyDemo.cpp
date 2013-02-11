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
#include "FluidSphSoftBodyDemo.h"

#include <stdio.h> 	//printf debugging

#include "GlutStuff.h"
#include "GLDebugFont.h"

#include "btBulletDynamicsCommon.h"
#include "LinearMath/btRandom.h"		//GEN_rand(), GEN_RAND_MAX

#include "BulletSoftBody/btSoftBodyHelpers.h"

#include "BulletFluids/Sph/btFluidSph.h"


FluidSphSoftBodyDemo::FluidSphSoftBodyDemo()
{
	setTexturing(true);
	setShadows(true);
	setCameraDistance(50.0);
	
	initPhysics();
}
FluidSphSoftBodyDemo::~FluidSphSoftBodyDemo() 
{
	exitPhysics(); 
}

void FluidSphSoftBodyDemo::initPhysics()
{
	m_collisionConfiguration = new btFluidSoftRigidCollisionConfiguration();

	m_dispatcher = new btCollisionDispatcher(m_collisionConfiguration);
	m_broadphase = new btDbvtBroadphase();
	m_solver = new btSequentialImpulseConstraintSolver();
	
	m_fluidSphSolver = new btFluidSphSolverDefault();
	
	m_dynamicsWorld = new btFluidSoftRigidDynamicsWorld(m_dispatcher, m_broadphase, m_solver, m_collisionConfiguration, m_fluidSphSolver);
	m_fluidSoftRigidWorld = static_cast<btFluidSoftRigidDynamicsWorld*>(m_dynamicsWorld);
	
	m_fluidSoftRigidWorld->setGravity( btVector3(0.0, -9.8, 0.0) );
	m_fluidSoftRigidWorld->setDebugDrawer(&m_debugDrawer);

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

		m_fluidSoftRigidWorld->addRigidBody(body);
	}
	
	//Create btFluidSph(s), which contain groups of particles
	{
		const int MAX_PARTICLES = 8192;
		btFluidSph* fluidSph = new btFluidSph(m_fluidSoftRigidWorld->getGlobalParameters(), MAX_PARTICLES);
		
		{
			btFluidSphParametersLocal FL = fluidSph->getLocalParameters();
			
			const btScalar AABB_EXTENT(25.0);
			FL.m_aabbBoundaryMin = btVector3(-AABB_EXTENT, -AABB_EXTENT, -AABB_EXTENT);
			FL.m_aabbBoundaryMax = btVector3(AABB_EXTENT, AABB_EXTENT, AABB_EXTENT);
			FL.m_enableAabbBoundary = 1;
			
			fluidSph->setLocalParameters(FL);
		}
		
		const bool ENABLE_CCD = true;
		if(ENABLE_CCD) fluidSph->setCcdMotionThreshold( fluidSph->getLocalParameters().m_particleRadius );
		
		m_fluidSoftRigidWorld->addFluidSph(fluidSph);
		m_fluidSph = fluidSph;
	}
	
	//Create a soft-body cloth patch
	{
		const btScalar POSITION_Y(10.0);
		const btScalar EXTENT(16.0);
		
		const int RESOLUTION = 17;
		
		btSoftBody*	softBody = btSoftBodyHelpers::CreatePatch( m_fluidSoftRigidWorld->getWorldInfo(), btVector3(-EXTENT,POSITION_Y,-EXTENT),
			btVector3(EXTENT,POSITION_Y,-EXTENT),
			btVector3(-EXTENT,POSITION_Y,EXTENT),
			btVector3(EXTENT,POSITION_Y,EXTENT),
			RESOLUTION, RESOLUTION, 1+2+4+8, true);
			
		btSoftBody::Material* material = softBody->appendMaterial();
		material->m_kLST = 0.4;
		material->m_flags -= btSoftBody::fMaterial::DebugDraw;
		
		softBody->m_cfg.kDF	= 1.0;
		softBody->m_cfg.kSRHR_CL = 1.0;
		softBody->m_cfg.kSR_SPLT_CL	= 0.5;
		softBody->m_cfg.collisions = btSoftBody::fCollision::CL_SS + btSoftBody::fCollision::CL_RS;
		
		softBody->setTotalMass(2.0);
		softBody->generateBendingConstraints(2, material);
		softBody->generateClusters(0); 	//Pass zero in generateClusters to create a cluster for each tetrahedron or triangle

		m_fluidSoftRigidWorld->addSoftBody(softBody);
	}
}
void FluidSphSoftBodyDemo::exitPhysics()
{
	//Cleanup in the reverse order of creation/initialization
	for(int i = m_fluidSoftRigidWorld->getNumCollisionObjects() - 1; i >= 0; i--)
	{
		btCollisionObject* obj = m_fluidSoftRigidWorld->getCollisionObjectArray()[i];
		btRigidBody* body = btRigidBody::upcast(obj);
		if( body && body->getMotionState() ) delete body->getMotionState();
		
		m_fluidSoftRigidWorld->removeCollisionObject(obj);
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
	delete m_fluidSoftRigidWorld;
	delete m_solver;
	delete m_broadphase;
	delete m_dispatcher;
	delete m_collisionConfiguration;
	delete m_fluidSphSolver;
}

void FluidSphSoftBodyDemo::clientMoveAndDisplay()
{
	btScalar secondsElapsed = getDeltaTimeMicroseconds() * btScalar(0.000001);
	
	const btFluidSphParametersGlobal& FG = m_fluidSoftRigidWorld->getGlobalParameters();
	if(m_fluidSoftRigidWorld)
	{
		const bool USE_SYNCRONIZED_TIME_STEP = false;	//Default: Rigid bodies == 16ms, Sph particles == 3ms time step
		if(USE_SYNCRONIZED_TIME_STEP)
		{
			btScalar timeStep = m_fluidSoftRigidWorld->getGlobalParameters().m_timeStep;
			m_fluidSoftRigidWorld->stepSimulation(timeStep, 0, timeStep);
		}
		else m_fluidSoftRigidWorld->stepSimulation( btScalar(1.0 / 60.0), 0 );
	}
	
	{
		static int counter = 0;
		if(++counter > 100)
		{
			counter = 0;
			
			for(int i = 0; i < m_fluidSoftRigidWorld->getNumFluidSph(); ++i)
				printf( "m_fluidSoftRigidWorld->getFluidSph(%d)->numParticles(): %d \n", i, m_fluidSoftRigidWorld->getFluidSph(i)->numParticles() );
		}
	}
		
	displayCallback();
}


GLuint generateSphereList(float radius)
{
	//Sphere generation code from FLUIDS v.2
	GLuint glSphereList = glGenLists(1);
	glNewList(glSphereList, GL_COMPILE);
		glBegin(GL_TRIANGLE_STRIP);
			for(float tilt = -90.0f; tilt <= 90.0f; tilt += 10.0f) 
			{
				for(float ang = 0.f; ang <= 360.0f; ang += 30.0f) 
				{
					const float DEGREES_TO_RADIANS = 3.141592f/180.0f;
					
					float ang_radians = ang * DEGREES_TO_RADIANS;
					float tilt_radians = tilt * DEGREES_TO_RADIANS;
					float tilt1_radians = (tilt + 10.0f) * DEGREES_TO_RADIANS;
				
					float x = sin(ang_radians) * cos(tilt_radians);
					float y = cos(ang_radians) * cos(tilt_radians);
					float z = sin(tilt_radians);
					float x1 = sin(ang_radians) * cos(tilt1_radians);
					float y1 = cos(ang_radians) * cos(tilt1_radians);
					float z1 = sin(tilt1_radians);
					
					glNormal3f(x, y, z);	glVertex3f(x*radius, y*radius, z*radius);		
					glNormal3f(x1, y1, z1);	glVertex3f(x1*radius, y1*radius, z1*radius);
				}
			}
		glEnd();
	glEndList();
	
	return glSphereList;
}
inline void drawSphere(GLuint glSphereList, const btVector3& position, float r, float g, float b)
{	
	glPushMatrix();
		glColor3f(r, g, b);
		glTranslatef( position.x(), position.y(), position.z() );
		glCallList(glSphereList);
	glPopMatrix();
}

void FluidSphSoftBodyDemo::displayCallback(void) 
{
	//BT_PROFILE() does not work correctly in this function;
	//timings are captured only when the camera is moving.
	//BT_PROFILE("FluidSphSoftBodyDemo::displayCallback()");

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 
	
	renderme();
	
	static bool areSpheresGenerated = false;
	static GLuint glLargeSphereList;
	static GLuint glSmallSphereList;
	if(!areSpheresGenerated)
	{
		areSpheresGenerated = true;
		
		const float RADIUS(1.0);	//Particle collision radius in btFluidSph::getLocalParameters() should also be changed if this is modified
		glSmallSphereList = generateSphereList(RADIUS*0.333f);
		glLargeSphereList = generateSphereList(RADIUS);
	}
	
	btIDebugDraw* debugDrawer = m_fluidSoftRigidWorld->getDebugDrawer();
	
	GLuint sphereList = glSmallSphereList;
	if(	debugDrawer && !(debugDrawer->getDebugMode() & btIDebugDraw::DBG_DrawWireframe)	) sphereList = glLargeSphereList;
	
	//Draw SPH particles
	{
		for(int i = 0; i < m_fluidSoftRigidWorld->getNumFluidSph(); ++i)
		{
			btFluidSph* fluidSph = m_fluidSoftRigidWorld->getFluidSph(i);
			for(int n = 0; n < fluidSph->numParticles(); ++n)
			{
				const btVector3& position = fluidSph->getPosition(n);
				drawSphere(sphereList, position, 0.6f, 0.9f, 1.0f);
			}
		}
	}
	
	//Draw soft bodies
	if(	debugDrawer && !(debugDrawer->getDebugMode() & btIDebugDraw::DBG_DrawWireframe)	)
	{	
		btAlignedObjectArray<btSoftBody*> softBodies = m_fluidSoftRigidWorld->getSoftBodyArray();
		for(int i = 0; i < softBodies.size(); ++i)
		{
			btSoftBody*	softBody = static_cast<btSoftBody*>(softBodies[i]);
			btSoftBodyHelpers::DrawFrame(softBody, debugDrawer);
			btSoftBodyHelpers::Draw( softBody, debugDrawer, m_fluidSoftRigidWorld->getDrawFlags() );
		}
	}
	
	if(m_fluidSoftRigidWorld) m_fluidSoftRigidWorld->debugDrawWorld();		//Optional but useful: debug drawing to detect problems
	
	if( (getDebugMode() & btIDebugDraw::DBG_NoHelpText) == 0 )
	{
		setOrthographicProjection();
		glDisable(GL_LIGHTING);
		glColor4f(0.0f, 0.0f, 0.0f, 1.0f);
		
		const int LINE_WIDTH = 360;
		int position_x = m_glutScreenWidth - LINE_WIDTH;
		int position_y = 20;
		
		glRasterPos3f(position_x, position_y, 0);
		GLDebugDrawString(position_x, position_y, "Press / to spray SPH particles");
		
		resetPerspectiveProjection();
		glEnable(GL_LIGHTING);
	}
	
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
void FluidSphSoftBodyDemo::keyboardCallback(unsigned char key, int x, int y)
{
	switch(key)
	{
		case ' ':
			clientResetScene();
			break;
	
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
			
		default:
			PlatformDemoApplication::keyboardCallback(key, x, y);
			break;
	}
}

void FluidSphSoftBodyDemo::specialKeyboard(int key, int x, int y)
{
	switch(key)
	{
		case GLUT_KEY_END:
		{
			int numObj = getDynamicsWorld()->getNumCollisionObjects();
			if(numObj)
			{
				btCollisionObject* obj = getDynamicsWorld()->getCollisionObjectArray()[numObj-1];
				if( btFluidSph::upcast(obj)  ) return;	//Deleting btFluidSph will cause crashes in FluidSphSoftBodyDemo
				
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

void FluidSphSoftBodyDemo::setShootBoxShape()
{
	if (!m_shootBoxShape)
	{
		const btScalar BOX_DIMENSIONS(3.0);
	
		btBoxShape* box = new btBoxShape( btVector3(BOX_DIMENSIONS, BOX_DIMENSIONS, BOX_DIMENSIONS) );
		box->initializePolyhedralFeatures();
		m_shootBoxShape = box;
	}
}

void FluidSphSoftBodyDemo::clientResetScene()
{
	exitPhysics();
	initPhysics();
}

