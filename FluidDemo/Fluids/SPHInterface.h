/** SPHInterface.h
	Copyright (C) 2012 Jackson Lee

	ZLib license
	This software is provided 'as-is', without any express or implied
	warranty. In no event will the authors be held liable for any damages
	arising from the use of this software.
	
	Permission is granted to anyone to use this software for any purpose,
	including commercial applications, and to alter it and redistribute it
	freely, subject to the following restrictions:
	
	1. The origin of this software must not be misrepresented; you must not
	   claim that you wrote the original software. If you use this software
	   in a product, an acknowledgment in the product documentation would be
	   appreciated but is not required.
	2. Altered source versions must be plainly marked as such, and must not be
	   misrepresented as being the original software.
	3. This notice may not be removed or altered from any source distribution.
*/
#ifndef SPHINTERFACE_H_INCLUDED
#define SPHINTERFACE_H_INCLUDED

#include "fluid_system.h"

#include "btBulletDynamicsCommon.h"
#include "LinearMath/btAlignedObjectArray.h"

struct ParticleResult : public btCollisionWorld::ContactResultCallback
{
	btVector3 m_normal;
	btCollisionObject *m_collisionObject;
	btScalar m_distance;
	
	btCollisionObject *m_particleObject;

	ParticleResult(btCollisionObject *particleObject) : m_particleObject(particleObject), m_collisionObject(0) {}
	
	virtual btScalar addSingleResult( btManifoldPoint &cp, const btCollisionObject *colObj0, int partId0, int index0,
														   const btCollisionObject *colObj1, int partId1, int index1 ) 
	{
		//Value returned from btCollisionWorld::ContactResultCallback::addSingleResult() appears to be unused
		const btScalar UNUSED = 1.0f;
	
		m_distance = cp.getDistance();
	
		if(m_particleObject == colObj0 && colObj1)		//Assume 0 == A
		{
			m_collisionObject = const_cast<btCollisionObject*>(colObj1);
			m_normal = cp.m_normalWorldOnB;
			
			return UNUSED;
		}
		else if(m_particleObject == colObj1 && colObj0)	//Assume 1 == B
		{
			//	this branch is never reached?
		
			m_collisionObject = const_cast<btCollisionObject*>(colObj0);
			m_normal = cp.m_normalWorldOnB;
			m_normal *= -1.0f;		//	verify that m_normal is the normal pushing m_particleObject away from m_collisionObject
			
			return UNUSED;
		}
		
		return UNUSED;
	}

	const bool hasHit() const { return static_cast<bool>(m_collisionObject); }
};

struct ParticleResultMulti : public btCollisionWorld::ContactResultCallback
{
	static const int MAX_COLLISIONS = 4;

	btVector3 m_normals[MAX_COLLISIONS];
	btCollisionObject *m_collisionObjects[MAX_COLLISIONS];
	btScalar m_distances[MAX_COLLISIONS];
	int m_numCollisions;
	
	btCollisionObject *m_particleObject;

	ParticleResultMulti(btCollisionObject *particleObject) : m_particleObject(particleObject), m_numCollisions(0) {}
	
	virtual btScalar addSingleResult( btManifoldPoint &cp, const btCollisionObject *colObj0, int partId0, int index0,
														   const btCollisionObject *colObj1, int partId1, int index1 ) 
	{
		//Value returned from btCollisionWorld::ContactResultCallback::addSingleResult() appears to be unused
		const btScalar UNUSED = 1.0f;
	
		if(m_numCollisions >= MAX_COLLISIONS) return UNUSED;
	
		if(m_particleObject == colObj0 && colObj1)		//Assume 0 == A
		{
			m_collisionObjects[m_numCollisions] = const_cast<btCollisionObject*>(colObj1);
			m_distances[m_numCollisions] = cp.getDistance();
			m_normals[m_numCollisions] = cp.m_normalWorldOnB;
			++m_numCollisions;
			
			return UNUSED;
		}
		else if(m_particleObject == colObj1 && colObj0)	//Assume 1 == B
		{
			//	this branch is never reached?
			btVector3 normal = cp.m_normalWorldOnB;
			normal *= -1.0f;		//	verify that m_normal is the normal pushing m_particleObject away from m_collisionObject
			
			m_collisionObjects[m_numCollisions] = const_cast<btCollisionObject*>(colObj0);
			m_distances[m_numCollisions] = cp.getDistance();
			m_normals[m_numCollisions] = normal;
			++m_numCollisions;
			
			return UNUSED;
		}
		
		return UNUSED;
	}
};

struct AabbTestCallback : public btBroadphaseAabbCallback
{
	btAlignedObjectArray<btCollisionObject*> m_results;

	virtual bool process(const btBroadphaseProxy *proxy) 
	{ 
		m_results.push_back( static_cast<btCollisionObject*>(proxy->m_clientObject) ); 
		
		//	investigate meaning of return value
		return true;
	}
};

class BulletFluidsInterface
{
	static void collideFluidsWithBullet(FluidSystem *fluidSystem, btCollisionWorld *world);
	static void collideFluidsWithBullet2(FluidSystem *fluidSystem, btCollisionWorld *world);
	static void collideFluidsWithBulletCcd(FluidSystem *fluidSystem, btCollisionWorld *world);

	static void resolveCollision(const FluidParameters &FP, Fluid *f, btCollisionObject *object, 
								 const btVector3 &fluidNormal, float distance);
								 
public:
	static void stepSimulation(FluidSystem *fluidSystem, btCollisionWorld *world, float secondsElapsed)
	{
		const bool USE_VARIABLE_DELTA_TIME = false;
		if(USE_VARIABLE_DELTA_TIME)
		{
			// crashes
			FluidParameters FP = fluidSystem->getParameters();
			FP.m_timeStep = secondsElapsed;
			fluidSystem->setParameters(FP);
		}
		
		//collideFluidsWithBullet(fluidSystem, world);
		//collideFluidsWithBullet2(fluidSystem, world);
		collideFluidsWithBulletCcd(fluidSystem, world);
		
		const bool USE_ACCUMULATOR = false;
		if(USE_ACCUMULATOR)
		{
			//Prevent simulation from running too quickly
			static float secondsAccumulated = 0;
			secondsAccumulated += secondsElapsed;
			
			if( secondsAccumulated > fluidSystem->getParameters().m_timeStep * 2.0 )
				secondsAccumulated = fluidSystem->getParameters().m_timeStep;
				
			if( secondsAccumulated >= fluidSystem->getParameters().m_timeStep )
			{
				fluidSystem->stepSimulation();
				secondsAccumulated -= fluidSystem->getParameters().m_timeStep;
			}
		}
		else fluidSystem->stepSimulation();
		
		{
			static int counter = 0;
			if(++counter > 100)
			{
				counter = 0;
				printf( "fluidSystem->numParticles(): %d \n", fluidSystem->numParticles() );
			}
		}
	}
};


struct ParticlesShape
{
	btScalar m_particleRadius;
	btScalar m_particleMass;
	btAlignedObjectArray<btVector3> m_points;
};
class BulletFluidsInterface_P
{
	static void getDynamicRigidBodies(FluidSystem *fluidSystem, btDynamicsWorld *world, btAlignedObjectArray<btRigidBody*> *out_rigidBodies);
	static void convertIntoParticles(btCollisionWorld *world, btCollisionObject *object, btScalar particleRadius, ParticlesShape *out_shape);
	static void runFluidSimulation( FluidSystem *fluidSystem, btDynamicsWorld *world, float secondsElapsed,  
									btAlignedObjectArray<ParticlesShape> *particles,
									btAlignedObjectArray<btRigidBody*> *rigidBodies );
									
public:	
	static void stepSimulation(FluidSystem *fluidSystem, btDynamicsWorld *world, float secondsElapsed)
	{
		//Find all non static rigid bodies
		btAlignedObjectArray<btRigidBody*> dynamicRigidBodies;
		getDynamicRigidBodies(fluidSystem, world, &dynamicRigidBodies);
		
		//Convert rigid bodies into particles
		btAlignedObjectArray<ParticlesShape> rigidBodiesAsParticles;
		{
			const float PARTICLE_RADIUS = 2.0f;
			for(int i = 0; i < dynamicRigidBodies.size(); ++i)
			{
				rigidBodiesAsParticles.push_back( ParticlesShape() );
				int lastIndex = rigidBodiesAsParticles.size() - 1;
				
				convertIntoParticles(world, dynamicRigidBodies[i], PARTICLE_RADIUS, &rigidBodiesAsParticles[lastIndex]);
			}
		}
		
		//Place particles at rigid body transform
		btAlignedObjectArray<ParticlesShape> initialParticles(rigidBodiesAsParticles);
		{
			for(int i = 0; i < initialParticles.size(); ++i)
				for(int n = 0; n < initialParticles[i].m_points.size(); ++n)
				{
					initialParticles[i].m_points[n] = dynamicRigidBodies[i]->getWorldTransform() * initialParticles[i].m_points[n];
				}
		}
	
		//Run fluid simulation
		btAlignedObjectArray<ParticlesShape> collidedParticles(initialParticles);
		runFluidSimulation(fluidSystem, world, secondsElapsed, &collidedParticles, &dynamicRigidBodies);
		
		//Use distance moved for each particle to determine the forces acting on the rigid body
		for(int i = 0; i < collidedParticles.size(); ++i)
		{
			btRigidBody *rigidBody = dynamicRigidBodies[i];
			
			btVector3 motion(0.0, 0.0, 0.0);
			for(int n = 0; n < collidedParticles[i].m_points.size(); ++n)
			{
				const btVector3 &startPosition = initialParticles[i].m_points[n];
				const btVector3 &endPosition = collidedParticles[i].m_points[n];
				motion += endPosition - startPosition;
			}
			
			//	test
			rigidBody->applyCentralForce(motion * 1000);
		}
	}
};

#endif

