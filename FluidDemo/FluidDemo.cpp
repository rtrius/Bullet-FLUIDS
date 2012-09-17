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
	m_demos.push_back( new Demo_Levee() );
	m_demos.push_back( new Demo_Drain() );
	m_demos.push_back( new Demo_DynamicBox() );
	m_demos.push_back( new Demo_HollowBox() );
	m_demos.push_back( new Demo_FluidToRigidBody() );
	m_demos.push_back( new Demo_MultiFluid() );
	m_demos.push_back( new Demo_ManyParticles() );
	m_demos.push_back( new Demo_Column() );
	
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
	btScalar ms = getDeltaTimeMicroseconds();
	btScalar secondsElapsed = ms * 0.000001f;
	
	const FluidParametersGlobal &FG = m_fluidWorld->getGlobalParameters();
	if(m_dynamicsWorld)
	{
		const bool USE_SYNCRONIZED_TIME_STEP = false;	//Default: btDynamicsWorld == 16ms, FluidWorld == 3ms time step
		if(USE_SYNCRONIZED_TIME_STEP)
		{
			btScalar timeStep = m_fluidWorld->getGlobalParameters().m_timeStep;
			m_dynamicsWorld->stepSimulation(timeStep, 0, timeStep);
			m_fluidWorld->stepSimulation();
		}
		else 
		{
			m_dynamicsWorld->stepSimulation(secondsElapsed);
			m_fluidWorld->stepSimulation();
		}
		
		{
			BT_PROFILE("FluidRigid Interaction");
			m_fluidRigidCollisionDetector.detectCollisions(FG, &m_fluidWorld->internalGetFluids(), m_dynamicsWorld);
			m_fluidRigidConstraintSolver.resolveCollisions( FG, &m_fluidRigidCollisionDetector.internalGetContactGroups() );
		}
		
		if( m_demos.size() ) m_demos[m_currentDemoIndex]->stepSimulation(*m_fluidWorld, &m_fluids);
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
void drawSphere(GLuint glSphereList, const btVector3 &position, float r, float g, float b)
{	
	glPushMatrix();
		//glEnable(GL_BLEND);
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
		//glColor4f(0.9f, 0.9f, 0.9f, 0.6f);
		
		glColor3f(r, g, b);
		glTranslatef( position.x(), position.y(), position.z() );
		glCallList(glSphereList);
	glPopMatrix();
}
void getFluidColors(bool drawFluidsWithMultipleColors, int fluidIndex, FluidSph *fluid, int particleIndex, float *out_r, float *out_g, float *out_b)
{
	const float COLOR_R = 0.3f;
	const float COLOR_G = 0.7f;
	const float COLOR_B = 1.0f;

	if(!drawFluidsWithMultipleColors)
	{
		float brightness = fluid->getVelocity(particleIndex).length() * 2.0f;
		if(brightness < 0.f)brightness = 0.f;
		if(brightness > 1.f)brightness = 1.f;
		
		*out_r = COLOR_R * brightness;
		*out_g = COLOR_G * brightness;
		*out_b = COLOR_B * brightness;
	}
	else
	{
		*out_r = COLOR_R; 
		*out_g = COLOR_G;
		*out_b = COLOR_B;
		
		if(fluidIndex % 2)
		{
			*out_r = 1.0f - COLOR_R;
			*out_g = 1.0f - COLOR_G;
			*out_b = 1.0f - COLOR_B;
		}
	}
}
void FluidDemo::displayCallback(void) 
{
	//BT_PROFILE() does not work correctly in this function;
	//timings are captured only when the camera is moving.
	//BT_PROFILE("FluidDemo::displayCallback()");

	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT); 

	renderme();

	static bool areSpheresGenerated = false;
	static GLuint glSmallSphereList;
	static GLuint glMediumSphereList;
	static GLuint glLargeSphereList;
	if(!areSpheresGenerated)
	{
		const FluidParametersGlobal &FG = m_fluidWorld->getGlobalParameters();
		btScalar particleRadius = FG.m_particleRadius / FG.m_simulationScale;
	
		areSpheresGenerated = true;
		glSmallSphereList = generateSphereList(0.1f);
		glMediumSphereList = generateSphereList(particleRadius * 0.3f);
		glLargeSphereList = generateSphereList(particleRadius);
	}
	bool drawFluidsWithMultipleColors = m_demos[m_currentDemoIndex]->isMultiFluidDemo();
	if(m_fluidRenderMode != FRM_MarchingCubes && m_fluidRenderMode != FRM_ScreenSpace)
	{
		//BT_PROFILE("Draw fluids - spheres");
		
		GLuint glSphereList;
		switch(m_fluidRenderMode)
		{
			case FRM_LargeSpheres:
				glSphereList = glLargeSphereList;
				break;
			case FRM_MediumSpheres:
				glSphereList = glMediumSphereList;
				break;
			case FRM_Points:
			default:
				glSphereList = glSmallSphereList;
				break;
		}
			
		for(int i = 0; i < m_fluids.size(); ++i)
			for(int n = 0; n < m_fluids[i]->numParticles(); ++n)
			{
				float r, g, b;
				getFluidColors(drawFluidsWithMultipleColors, i, m_fluids[i], n, &r, &g, &b);
				
				drawSphere(glSphereList, m_fluids[i]->getPosition(n), r, g, b);
			}
	}
	else if(m_fluidRenderMode == FRM_ScreenSpace)
	{
		if(m_screenSpaceRenderer)
		{
			if(!m_ortho)
			{
				const FluidParametersGlobal &FG = m_fluidWorld->getGlobalParameters();
				btScalar particleRadius = FG.m_particleRadius / FG.m_simulationScale;
				
				for(int i = 0; i < m_fluids.size(); ++i)
				{
					float r = 0.5f;
					float g = 0.8f;
					float b = 1.0f;
					
					//Beer's law constants
					//Controls the darkening of the fluid's color based on its thickness
					//For a constant k, (k > 1) == darkens faster; (k < 1) == darkens slower; (k == 0) == disable
					float absorptionR = 2.5;	
					float absorptionG = 1.0;
					float absorptionB = 0.5;
		
					if(drawFluidsWithMultipleColors)
					{
						r = 0.3f; 
						g = 0.7f;
						b = 1.0f;
						if(i % 2)
						{
							r = 1.0f - r;
							g = 1.0f - g;
							b = 1.0f - b;
						}
						
						absorptionR = 1.0;
						absorptionG = 1.0;
						absorptionB = 1.0;
					}
					
					m_screenSpaceRenderer->render(m_fluids[i]->internalGetFluidParticles().m_pos, particleRadius,
												  r, g, b, absorptionR, absorptionG, absorptionB);
				}
			}
			else printf("Orthogonal rendering not implemented for ScreenSpaceFluidRendererGL.\n");
		}
	}
	else 	//(m_fluidRenderMode == FRM_MarchingCubes)
	{
		//BT_PROFILE("Draw fluids - marching cubes");
			
		const int CELLS_PER_EDGE = 32;
		static MarchingCubes *marchingCubes = 0;
		if(!marchingCubes) 
		{
			marchingCubes = new MarchingCubes;
			marchingCubes->initialize(CELLS_PER_EDGE);
		}
		
		for(int i = 0; i < m_fluids.size(); ++i)
		{
			marchingCubes->generateMesh(*m_fluids[i]);
			
			const btAlignedObjectArray<float> &vertices = marchingCubes->getTriangleVertices();
			if( vertices.size() )
			{
				glEnable(GL_BLEND);
				glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 

				float r = 0.9f;
				float g = 0.9f;
				float b = 0.9f;
				if(drawFluidsWithMultipleColors)
				{
					r = 0.3f; 
					g = 0.7f;
					b = 1.0f;
					if(i % 2)
					{
						r = 1.0f - r;
						g = 1.0f - g;
						b = 1.0f - b;
					}
				}
				
				glEnableClientState(GL_VERTEX_ARRAY);
					glColor4f(r, g, b, 0.6f);
					glVertexPointer( 3, GL_FLOAT, 0, &vertices[0] );
					glDrawArrays( GL_TRIANGLES, 0, vertices.size()/3 );	
				glDisableClientState(GL_VERTEX_ARRAY);
			}
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
		{
			int currentRenderMode = static_cast<int>(m_fluidRenderMode);
			m_fluidRenderMode = static_cast<FluidRenderMode>(currentRenderMode + 1);
		
			if(m_fluidRenderMode > FRM_MarchingCubes)m_fluidRenderMode = FRM_Points;
			return;
		}
		
		case 'e':
			if(m_fluidSolverGPU)
			{
				m_useFluidSolverOpenCL = !m_useFluidSolverOpenCL;
				
				if(m_useFluidSolverOpenCL)m_fluidWorld->setFluidSolver(m_fluidSolverGPU);
				else m_fluidWorld->setFluidSolver(m_fluidSolverCPU );
			}
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
			if(m_maxFluidParticles*2 <= MAX_FLUID_PARTICLES) m_maxFluidParticles *= 2;
			resetCurrentDemo();
			return;
	}
	
	PlatformDemoApplication::keyboardCallback(key, x, y);
}


