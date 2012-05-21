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

#include "btBulletDynamicsCommon.h"
#include "Fluids/fluid_system.h"

#include <ctime>

struct ParticleResult : public btCollisionWorld::ContactResultCallback
{
	btCollisionObject *m_particleObject;
	btCollisionObject *m_collisionObject;
	
	btScalar m_distance;
	btVector3 m_normal;

	ParticleResult(btCollisionObject *particleObject) : m_particleObject(particleObject), m_collisionObject(0) {}
	
	virtual btScalar addSingleResult( btManifoldPoint &cp, const btCollisionObject *colObj0, int partId0, int index0,
														   const btCollisionObject *colObj1, int partId1, int index1 ) 
	{
		const float COLLIDING = 0.5;		//	arbitrary value; verify
		const float NOT_COLLIDING = 1.0;	//	verify 1.0 == invalid
	
		m_distance = cp.getDistance();
	
		if(m_particleObject == colObj0 && colObj1)		//Assume 0 == A
		{
			m_collisionObject = const_cast<btCollisionObject*>(colObj1);
			m_normal = cp.m_normalWorldOnB;
			
			return COLLIDING;
		}
		else if(m_particleObject == colObj1 && colObj0)	//Assume 1 == B
		{
			//	this branch is never reached?
		
			m_collisionObject = const_cast<btCollisionObject*>(colObj0);
			m_normal = cp.m_normalWorldOnB;
			m_normal *= -1.0f;		//	verify that m_normal is the normal pushing m_particleObject away from m_collisionObject
			
			return COLLIDING;
		}
		
		return NOT_COLLIDING;
	}

	const bool hasHit() const { return static_cast<bool>(m_collisionObject); }
	btCollisionObject* getCollidedWith() const { return m_collisionObject; }
	const btVector3& getNormal() const { return m_normal; }
	const btScalar getDistance() const { return m_distance; }
};

struct AabbTestCallback : public btBroadphaseAabbCallback
{
	std::vector<btCollisionObject *const> *m_results;

	AabbTestCallback(std::vector<btCollisionObject* const> *results) : m_results(results) {}
	virtual bool process(const btBroadphaseProxy *proxy) 
	{ 
		m_results->push_back( static_cast<btCollisionObject *const>(proxy->m_clientObject) ); 
		
		//	investigate meaning of return value
		return true;
	}
};

//Demonstrate how to use FluidSystem
struct FluidsWrapper
{
	FluidSystem m_fluidSystem;

	void initialize()
	{
		const int MAX_PARTICLES = 16384;
		const float VOL_BOUND = 30.0f;
		Vector3DF volMin(-VOL_BOUND, -10.0f, -VOL_BOUND);
		Vector3DF volMax(VOL_BOUND, VOL_BOUND*2.0f, VOL_BOUND);
		
		m_fluidSystem.initialize(MAX_PARTICLES, volMin, volMax);
		
		//Spawn particles
		const float INIT_BOUND = 20.0f;
		Vector3DF initMin(-INIT_BOUND, 20.f, -INIT_BOUND);
		//Vector3DF initMax(INIT_BOUND, 30.0f, INIT_BOUND);
		Vector3DF initMax(INIT_BOUND, 55.0f, INIT_BOUND);
		FluidEmitter::addVolume( &m_fluidSystem, initMin, initMax, m_fluidSystem.getEmitterSpacing() * 0.87 );
		
		//Customize parameters
		FluidParameters FP = m_fluidSystem.getParameters();
		FP.m_planeGravity = Vector3DF(0.f, -9.8f, 0.f);
		m_fluidSystem.setParameters(FP);
	}
	
	void toggleOpenCL() { m_fluidSystem.toggleOpenCL(); }
	

	
	void stepSimulation(btCollisionWorld *world, float secondsElapsed)
	{
		processEmittersAndAbsorbers();
		
		const bool USE_VARIABLE_DELTA_TIME = false;
		if(USE_VARIABLE_DELTA_TIME)
		{
			// crashes
			FluidParameters FP = m_fluidSystem.getParameters();
			FP.m_timeStep = secondsElapsed;
			m_fluidSystem.setParameters(FP);
		}
		
		//	implement time accumulator?
		m_fluidSystem.stepSimulation();
		
		collideFluidsWithBullet(world);
		//collideFluidsWithBullet2(world);
		
		{
			static int counter = 0;
			if(++counter > 100)
			{
				counter = 0;
				printf( "m_fluidSystem.numParticles(): %d \n", m_fluidSystem.numParticles() );
			}
		}
	}
	
public:
	void processEmittersAndAbsorbers()
	{
		FluidEmitter emitter;
		emitter.m_pitch = -90.0f;
		emitter.m_velocity = 1.f;
		emitter.m_yawSpread = 1.f;
		emitter.m_pitchSpread = 1.f;
		emitter.m_position.setValue(10.f, 10.f, 10.f);
		
		const int FRAMES_BETWEEN_EMIT = 2;
		static unsigned int frame = 0;

		if( FRAMES_BETWEEN_EMIT > 0 && (++frame) % FRAMES_BETWEEN_EMIT == 0 )
		{
			const int NUM_PARTICLES_EMITTED = 5;
			emitter.emit( &m_fluidSystem, NUM_PARTICLES_EMITTED, m_fluidSystem.getEmitterSpacing() );
		}
		
		//
		const float BOUND = 3.0f;
		FluidAbsorber absorber;
		absorber.m_min.setValue(-BOUND, -BOUND, -BOUND);
		absorber.m_max.setValue(BOUND, BOUND, BOUND);
		//absorber.absorb(&m_fluidSystem);
	}

	void collideFluidsWithBullet(btCollisionWorld *world)
	{
		BT_PROFILE("FluidsWrapper::collideFluidsWithBullet()");
		
		//sph_pradius is at simscale; divide by simscale to transform into world scale
		//btSphereShape particleShape( m_fluidSystem.getParameters().sph_pradius / m_fluidSystem.getParameters().sph_simscale );
		btSphereShape particleShape( m_fluidSystem.getParameters().sph_pradius  );
		
		btCollisionObject particleObject;
		particleObject.setCollisionShape(&particleShape);
			
		btTransform &particleTransform = particleObject.getWorldTransform();
		particleTransform.setRotation( btQuaternion::getIdentity() );
			
		for(int i = 0; i < m_fluidSystem.numParticles(); ++i)
		{
			Fluid *pFluid = m_fluidSystem.getFluid(i);
			btVector3 particlePosition(pFluid->pos.x(), pFluid->pos.y(), pFluid->pos.z());
			
			particleTransform.setOrigin(particlePosition);
				
			ParticleResult result(&particleObject);
			world->contactTest(&particleObject, result);
				
			if( result.hasHit() ) resolveFluidCollision(i, result);
		}
	}
	
	void collideFluidsWithBullet2(btCollisionWorld *world)
	{
		//Sample AABB from all fluid particles, detect intersecting 
		//broadphase proxies, and test each btCollisionObject
		//against all fluid grid cells intersecting with the btCollisionObject's AABB.
		
		BT_PROFILE("FluidsWrapper::collideFluidsWithBullet2()");

		const float LARGE_FLOAT = 1000000.f;
		btVector3 fluidSystemMin(LARGE_FLOAT, LARGE_FLOAT, LARGE_FLOAT);
		btVector3 fluidSystemMax(-LARGE_FLOAT, -LARGE_FLOAT, -LARGE_FLOAT);
		for(int i = 0; i < m_fluidSystem.numParticles(); ++i)
		{
			//	check CLRS for faster method
			const Vector3DF &pos = m_fluidSystem.getFluid(i)->pos;
			if( pos.x() < fluidSystemMin.x() ) fluidSystemMin.m_floats[0] = pos.x();
			if( pos.y() < fluidSystemMin.y() ) fluidSystemMin.m_floats[1] = pos.y();
			if( pos.z() < fluidSystemMin.z() ) fluidSystemMin.m_floats[2] = pos.z();
			
			if( pos.x() > fluidSystemMax.x() ) fluidSystemMax.m_floats[0] = pos.x();
			if( pos.y() > fluidSystemMax.y() ) fluidSystemMax.m_floats[1] = pos.y();
			if( pos.z() > fluidSystemMax.z() ) fluidSystemMax.m_floats[2] = pos.z();		
		}
		
		//
		//float particleRadius = m_fluidSystem.getParameters().sph_pradius / m_fluidSystem.getParameters().sph_simscale;
		float particleRadius = m_fluidSystem.getParameters().sph_pradius;
		btSphereShape particleShape(particleRadius);
		
		btCollisionObject particleObject;
		particleObject.setCollisionShape(&particleShape);
						
		btTransform &particleTransform = particleObject.getWorldTransform();
		particleTransform.setRotation( btQuaternion::getIdentity() );
		
		//
		const Grid &G = m_fluidSystem.getGrid();
		const GridParameters &GP = G.getParameters();
		std::vector<btCollisionObject *const> inFluidsAabb;
		
		AabbTestCallback callback(&inFluidsAabb);
		
		world->getBroadphase()->aabbTest(fluidSystemMin, fluidSystemMax, callback);
		for(int i = 0; i < inFluidsAabb.size(); ++i)
		{
			btCollisionObject *const object = inFluidsAabb[i];
		
			const btVector3 &objectMin = object->getBroadphaseHandle()->m_aabbMin;
			const btVector3 &objectMax = object->getBroadphaseHandle()->m_aabbMax;
			
			int gridMinX, gridMinY, gridMinZ;
			int gridMaxX, gridMaxY, gridMaxZ;
			G.getIndicies( Vector3DF( objectMin.x(), objectMin.y(), objectMin.z() ), &gridMinX, &gridMinY, &gridMinZ );
			G.getIndicies( Vector3DF( objectMax.x(), objectMax.y(), objectMax.z() ), &gridMaxX, &gridMaxY, &gridMaxZ );
			
			for(int z = gridMinZ; z <= gridMaxZ; ++z)
				for(int y = gridMinY; y <= gridMaxY; ++y)
					for(int x = gridMinX; x <= gridMaxX; ++x)
					{
						int currentIndex = G.getLastParticleIndex( (z*GP.m_resolutionY + y)*GP.m_resolutionX + x );
						while(currentIndex != INVALID_PARTICLE_INDEX)
						{
							Fluid *f = m_fluidSystem.getFluid(currentIndex);
							Vector3DF fluidMin( f->pos.x() - particleRadius, f->pos.y() - particleRadius, f->pos.z() - particleRadius );
							Vector3DF fluidMax( f->pos.x() + particleRadius, f->pos.y() + particleRadius, f->pos.z() + particleRadius );
							
							if( fluidMin.x() <= objectMax.x() && objectMin.x() <= fluidMax.x()  
							 && fluidMin.y() <= objectMax.y() && objectMin.y() <= fluidMax.y()
							 && fluidMin.z() <= objectMax.z() && objectMin.z() <= fluidMax.z() )
							{
								btVector3 particlePosition( f->pos.x(), f->pos.y(), f->pos.z() );
								particleTransform.setOrigin(particlePosition);
								
								ParticleResult result(&particleObject);
								world->contactPairTest(&particleObject, const_cast<btCollisionObject*>(object), result);
							
								if( result.hasHit() ) resolveFluidCollision(currentIndex, result);
							}
							
							currentIndex = f->nextFluidIndex;
						}
					}
		}
	}
	
	
	void resolveFluidCollision(int particleIndex, const ParticleResult &collisionResult)
	{
		const float stiff = m_fluidSystem.getParameters().sph_extstiff;
		const float damp = m_fluidSystem.getParameters().sph_extdamp;
		const float radius = m_fluidSystem.getParameters().sph_pradius;
		const float simScale = m_fluidSystem.getParameters().sph_simscale;
		const float COLLISION_EPSILON = 0.00001f;
		
		Fluid *f = m_fluidSystem.getFluid(particleIndex);
		
		float diff = 2.0f * radius - collisionResult.getDistance()*simScale;
		if(diff > COLLISION_EPSILON)
		{
			const btVector3& colNormal = collisionResult.getNormal();
			Vector3DF normal( colNormal.x(), colNormal.y(), colNormal.z() );
			
			//	Alternate method: push the particle out of the object, cancel acceleration
			//	Alternate method: reflect the particle's velocity along the normal (issues with low velocities)
			//Current method: accelerate the particle to avoid collision
			float adj = stiff * diff - damp * normal.dot(f->vel_eval);
			
			//Since the collision is resolved(by pushing the particle out of btCollisionObject) in 1 frame,
			//particles 'jump' when they collide. On average 2-3 times the needed acceleration is applied,
			//so we scale the acceleration to make the collision appear more smooth.
			//const float SCALE = 0.3f;
			const float SCALE = 0.5f;
			Vector3DF acceleration( adj * normal.x() * SCALE,
									adj * normal.y() * SCALE,
									adj * normal.z() * SCALE );
				
			//if externalAcceleration is very high, the fluid simulation will explode
			f->externalAcceleration = acceleration;
			
			/*
			//btRigidBody-Fluid interaction
			btRigidBody *pRigidBody = btRigidBody::upcast( collisionResult.getCollidedWith() );
			if(pRigidBody)
			{
				const float fluidMass = m_fluidSystem.getParameters().sph_pmass;
				btVector3 fluidPosition(f->pos.x(), f->pos.y(), f->pos.z());
				btVector3 fluidForce(f->sph_force.x(), f->sph_force.y(), f->sph_force.z());
			
				const float invertedMass = pRigidBody->getInvMass();
			
				//Rigid body to particle
				//const btVector3& linearVelocity = pRigidBody->getLinearVelocity();
				//acceleration.x() += ;
				//acceleration.y() += ;
				//acceleration.z() += ;
				
				//Particle to rigid body
				//	probably incorrect
				//pRigidBody->applyForce( fluidForce, fluidPosition - pRigidBody->getWorldTransform().getOrigin() );
			}
			*/
		}
	}

};


#endif

