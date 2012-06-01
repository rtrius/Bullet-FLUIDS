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

struct BulletFluidsInterface
{
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
		
		//	implement time accumulator?
		fluidSystem->stepSimulation();
		
		collideFluidsWithBullet(fluidSystem, world);
		//collideFluidsWithBullet2(fluidSystem, world);
		//collideFluidsWithBulletCcd(fluidSystem, world);
		
		{
			static int counter = 0;
			if(++counter > 100)
			{
				counter = 0;
				printf( "fluidSystem->numParticles(): %d \n", fluidSystem->numParticles() );
			}
		}
	}

	static void collideFluidsWithBulletCcd(FluidSystem *fluidSystem, btCollisionWorld *world)
	{
		BT_PROFILE("BulletFluidsInterface::collideFluidsWithBulletCcd()");
		
		//from collideFluidsWithBullet()
		btSphereShape particleShape( fluidSystem->getParameters().sph_pradius );
		
		btCollisionObject particleObject;
		particleObject.setCollisionShape(&particleShape);
			
		btTransform &particleTransform = particleObject.getWorldTransform();
		particleTransform.setRotation( btQuaternion::getIdentity() );
		//from collideFluidsWithBullet()
		
		//int numCollided = 0;
		for(int i = 0; i < fluidSystem->numParticles(); ++i)
		{
			Fluid *f = fluidSystem->getFluid(i);
			
			//	investigate:
			//	Calling btCollisionWorld::rayTest() with the ray's origin inside a btCollisionObject
			//	reports no collisions; is this expected behavior when using btCollisionWorld::ClosestRayResultCallback?
			//btCollisionWorld::ClosestRayResultCallback result( btVector3(0.f, -1.0f, 0.f), btVector3(0.f, -2.0f, 0.f) );
			btCollisionWorld::ClosestRayResultCallback result(f->prev_pos, f->pos);
			result.m_collisionFilterGroup = btBroadphaseProxy::DefaultFilter;
			result.m_collisionFilterMask = btBroadphaseProxy::AllFilter;
			
			world->rayTest(result.m_rayFromWorld, result.m_rayToWorld, result);
			
			//printf( "result.m_closestHitFraction, result.hasHit(): %f, %d \n", result.m_closestHitFraction, result.hasHit() );
			
			if( result.hasHit() )
			{
				//++numCollided;
				resolveFluidCollisionCcd( fluidSystem->getParameters(), result, fluidSystem->getFluid(i) );
			}
			else
			{
				//Since btCollisionWorld::rayTest() does not report a collision if the ray's origin
				//is inside a btCollsionObject, use btCollisionWorld::contactTest() in case
				//the particle is inside.
				
				//from collideFluidsWithBullet()
				particleTransform.setOrigin(f->pos);
					
				ParticleResult result(&particleObject);
				world->contactTest(&particleObject, result);
					
				if( result.hasHit() ) resolveFluidCollision( fluidSystem->getParameters(), result, fluidSystem->getFluid(i) );
				//from collideFluidsWithBullet()
			}
		}
		
		//printf("numCollided: %d \n", numCollided);
	}	
	static void resolveFluidCollisionCcd(const FluidParameters &FP, 
										 const btCollisionWorld::ClosestRayResultCallback &result, Fluid *f)
	{
		const float stiff = FP.sph_extstiff;
		const float damp = FP.sph_extdamp;
		const float radius = FP.sph_pradius;
		const float simScale = FP.sph_simscale;
		const float COLLISION_EPSILON = 0.00001f;
		
		float distanceMoved = result.m_rayFromWorld.distance(result.m_rayToWorld);
		float distanceCollided = result.m_rayFromWorld.distance(result.m_hitPointWorld);
		
		float depthOfPenetration = distanceMoved - distanceCollided;
		
		float diff = 2.0f * radius - depthOfPenetration*simScale;
		if(diff > COLLISION_EPSILON)
		{
			const btVector3& normal = result.m_hitNormalWorld;
			
			//	Alternate method: push the particle out of the object, cancel acceleration
			//	Alternate method: reflect the particle's velocity along the normal (issues with low velocities)
			//Current method: accelerate the particle to avoid collision
			float adj = stiff * diff - damp * normal.dot(f->vel_eval);
			
			//Since the collision is resolved(by pushing the particle out of btCollisionObject) in 1 frame,
			//particles 'jump' when they collide. On average 2-3 times the needed acceleration is applied,
			//so we scale the acceleration to make the collision appear more smooth.
			//const float SCALE = 0.3f;
			const float SCALE = 0.5f;
			btVector3 acceleration( adj * normal.x() * SCALE,
									adj * normal.y() * SCALE,
									adj * normal.z() * SCALE );
				
			//if externalAcceleration is very high, the fluid simulation will explode
			f->externalAcceleration = acceleration;
			
			f->pos = result.m_hitPointWorld;
		}
	}
	
	
	static void collideFluidsWithBullet(FluidSystem *fluidSystem, btCollisionWorld *world)
	{
		BT_PROFILE("BulletFluidsInterface::collideFluidsWithBullet()");
		
		//sph_pradius is at simscale; divide by simscale to transform into world scale
		//btSphereShape particleShape( fluidSystem->getParameters().sph_pradius / fluidSystem->getParameters().sph_simscale );
		btSphereShape particleShape( fluidSystem->getParameters().sph_pradius  );
		
		btCollisionObject particleObject;
		particleObject.setCollisionShape(&particleShape);
			
		btTransform &particleTransform = particleObject.getWorldTransform();
		particleTransform.setRotation( btQuaternion::getIdentity() );
			
		for(int i = 0; i < fluidSystem->numParticles(); ++i)
		{
			Fluid *pFluid = fluidSystem->getFluid(i);
			
			particleTransform.setOrigin(pFluid->pos);
				
			ParticleResult result(&particleObject);
			world->contactTest(&particleObject, result);
				
			if( result.hasHit() ) resolveFluidCollision( fluidSystem->getParameters(), result, fluidSystem->getFluid(i) );
		}
	}
	
	static void collideFluidsWithBullet2(FluidSystem *fluidSystem, btCollisionWorld *world)
	{
		//Sample AABB from all fluid particles, detect intersecting 
		//broadphase proxies, and test each btCollisionObject
		//against all fluid grid cells intersecting with the btCollisionObject's AABB.
		
		BT_PROFILE("BulletFluidsInterface::collideFluidsWithBullet2()");

		const float LARGE_FLOAT = 1000000.f;
		btVector3 fluidSystemMin(LARGE_FLOAT, LARGE_FLOAT, LARGE_FLOAT);
		btVector3 fluidSystemMax(-LARGE_FLOAT, -LARGE_FLOAT, -LARGE_FLOAT);
		for(int i = 0; i < fluidSystem->numParticles(); ++i)
		{
			//	check CLRS for faster method
			const btVector3 &pos = fluidSystem->getFluid(i)->pos;
			if( pos.x() < fluidSystemMin.x() ) fluidSystemMin.m_floats[0] = pos.x();
			if( pos.y() < fluidSystemMin.y() ) fluidSystemMin.m_floats[1] = pos.y();
			if( pos.z() < fluidSystemMin.z() ) fluidSystemMin.m_floats[2] = pos.z();
			
			if( pos.x() > fluidSystemMax.x() ) fluidSystemMax.m_floats[0] = pos.x();
			if( pos.y() > fluidSystemMax.y() ) fluidSystemMax.m_floats[1] = pos.y();
			if( pos.z() > fluidSystemMax.z() ) fluidSystemMax.m_floats[2] = pos.z();		
		}
		
		//
		//float particleRadius = fluidSystem->getParameters().sph_pradius / fluidSystem->getParameters().sph_simscale;
		float particleRadius = fluidSystem->getParameters().sph_pradius;
		btSphereShape particleShape(particleRadius);
		
		btCollisionObject particleObject;
		particleObject.setCollisionShape(&particleShape);
						
		btTransform &particleTransform = particleObject.getWorldTransform();
		particleTransform.setRotation( btQuaternion::getIdentity() );
		
		//
		const Grid &G = fluidSystem->getGrid();
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
			G.getIndicies(objectMin, &gridMinX, &gridMinY, &gridMinZ);
			G.getIndicies(objectMax, &gridMaxX, &gridMaxY, &gridMaxZ);
			
			for(int z = gridMinZ; z <= gridMaxZ; ++z)
				for(int y = gridMinY; y <= gridMaxY; ++y)
					for(int x = gridMinX; x <= gridMaxX; ++x)
					{
						int currentIndex = G.getLastParticleIndex( (z*GP.m_resolutionY + y)*GP.m_resolutionX + x );
						while(currentIndex != INVALID_PARTICLE_INDEX)
						{
							Fluid *f = fluidSystem->getFluid(currentIndex);
							btVector3 fluidMin( f->pos.x() - particleRadius, f->pos.y() - particleRadius, f->pos.z() - particleRadius );
							btVector3 fluidMax( f->pos.x() + particleRadius, f->pos.y() + particleRadius, f->pos.z() + particleRadius );
							
							if( fluidMin.x() <= objectMax.x() && objectMin.x() <= fluidMax.x()  
							 && fluidMin.y() <= objectMax.y() && objectMin.y() <= fluidMax.y()
							 && fluidMin.z() <= objectMax.z() && objectMin.z() <= fluidMax.z() )
							{
								particleTransform.setOrigin(f->pos);
								
								ParticleResult result(&particleObject);
								world->contactPairTest(&particleObject, const_cast<btCollisionObject*>(object), result);
							
								if( result.hasHit() ) resolveFluidCollision( fluidSystem->getParameters(), result, 
																			fluidSystem->getFluid(currentIndex) );
							}
							
							currentIndex = f->nextFluidIndex;
						}
					}
		}
	}
	
	
	static void resolveFluidCollision(const FluidParameters &FP, const ParticleResult &collisionResult, Fluid *f)
	{
		const float stiff = FP.sph_extstiff;
		const float damp = FP.sph_extdamp;
		const float radius = FP.sph_pradius;
		const float simScale = FP.sph_simscale;
		const float COLLISION_EPSILON = 0.00001f;
		
		float diff = 2.0f * radius - collisionResult.getDistance()*simScale;
		if(diff > COLLISION_EPSILON)
		{
			const btVector3& normal = collisionResult.getNormal();
			
			//	Alternate method: push the particle out of the object, cancel acceleration
			//	Alternate method: reflect the particle's velocity along the normal (issues with low velocities)
			//Current method: accelerate the particle to avoid collision
			float adj = stiff * diff - damp * normal.dot(f->vel_eval);
			
			//Since the collision is resolved(by pushing the particle out of btCollisionObject) in 1 frame,
			//particles 'jump' when they collide. On average 2-3 times the needed acceleration is applied,
			//so we scale the acceleration to make the collision appear more smooth.
			//const float SCALE = 0.3f;
			const float SCALE = 0.5f;
			btVector3 acceleration( adj * normal.x() * SCALE,
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
			
				const float invertedMass = pRigidBody->getInvMass();
			
				//Rigid body to particle
				//const btVector3& linearVelocity = pRigidBody->getLinearVelocity();
				//acceleration.x() += ;
				//acceleration.y() += ;
				//acceleration.z() += ;
				
				//Particle to rigid body
				//	probably incorrect
				//pRigidBody->applyForce( f->sph_force, f->pos - pRigidBody->getWorldTransform().getOrigin() );
			}
			*/
		}
	}

};


#endif

