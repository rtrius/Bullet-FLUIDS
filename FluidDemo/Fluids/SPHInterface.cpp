/** SPHInterface.cpp
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

#include "SPHInterface.h"

void getAabb(const FluidParametersGlobal &FG, const FluidSph &F, btVector3 *out_min, btVector3 *out_max)
{
	btVector3 min(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);
	btVector3 max(-BT_LARGE_FLOAT, -BT_LARGE_FLOAT, -BT_LARGE_FLOAT);
	
	for(int i = 0; i < F.numParticles(); ++i)
	{
		const btVector3 &pos = F.getPosition(i);
		if( pos.x() < min.x() ) min.setX( pos.x() );
		if( pos.y() < min.y() ) min.setY( pos.y() );
		if( pos.z() < min.z() ) min.setZ( pos.z() );
		
		if( pos.x() > max.x() ) max.setX( pos.x() );
		if( pos.y() > max.y() ) max.setY( pos.y() );
		if( pos.z() > max.z() ) max.setZ( pos.z() );
	}
	
	//
	const btScalar particleRadius = FG.sph_pradius / FG.sph_simscale;
	for(int i = 0; i < 3; ++i)min.m_floats[i] -= particleRadius;
	for(int i = 0; i < 3; ++i)max.m_floats[i] += particleRadius;
	
	//
	*out_min = min;
	*out_max = max;
}
void BulletFluidsInterface::collideFluidsWithBullet(const FluidParametersGlobal &FG, FluidSph *fluid, btCollisionWorld *world)
{
	BT_PROFILE("BulletFluidsInterface::collideFluidsWithBullet()");
	
	//sph_pradius is at simscale; divide by simscale to transform into world scale
	btSphereShape particleShape( FG.sph_pradius / FG.sph_simscale );
	
	btCollisionObject particleObject;
	particleObject.setCollisionShape(&particleShape);
		
	btTransform &particleTransform = particleObject.getWorldTransform();
	particleTransform.setRotation( btQuaternion::getIdentity() );
		
	for(int i = 0; i < fluid->numParticles(); ++i)
	{
		particleTransform.setOrigin( fluid->getPosition(i) );
			
		ParticleResultMulti result(&particleObject);
		world->contactTest(&particleObject, result);
			
		for(int n = 0; n < result.m_numCollisions; ++n)
		{
			resolveCollision(FG, fluid, i, result.m_collidedWith[n], 
							 result.m_normals[n], result.m_collidedWithHitPointsWorld[n], result.m_distances[n]);
		}
	}
}

void BulletFluidsInterface::collideFluidsWithBullet2(const FluidParametersGlobal &FG, FluidSph *fluid, btCollisionWorld *world)
{
	//Sample AABB from all fluid particles, detect intersecting 
	//broadphase proxies, and test each btCollisionObject
	//against all fluid grid cells intersecting with the btCollisionObject's AABB.
	
	BT_PROFILE("BulletFluidsInterface::collideFluidsWithBullet2()");

	btVector3 fluidSystemMin, fluidSystemMax;
	getAabb(FG, *fluid, &fluidSystemMin, &fluidSystemMax);
	
	//sph_pradius is at simscale; divide by simscale to transform into world scale
	btScalar particleRadius = FG.sph_pradius / FG.sph_simscale;
	btSphereShape particleShape(particleRadius);
	
	btCollisionObject particleObject;
	particleObject.setCollisionShape(&particleShape);
					
	btTransform &particleTransform = particleObject.getWorldTransform();
	particleTransform.setRotation( btQuaternion::getIdentity() );
	
	//
	const FluidGrid *grid = fluid->getGrid();
	const bool isLinkedList = (grid->getGridType() == FT_LinkedList);
	
	AabbTestCallback aabbTest;
	world->getBroadphase()->aabbTest(fluidSystemMin, fluidSystemMax, aabbTest);
	
	for(int i = 0; i < aabbTest.m_results.size(); ++i)
	{
		btCollisionObject *const object = aabbTest.m_results[i];
	
		const btVector3 &objectMin = object->getBroadphaseHandle()->m_aabbMin;
		const btVector3 &objectMax = object->getBroadphaseHandle()->m_aabbMax;
		
		btAlignedObjectArray<int> gridCellIndicies;
		grid->getGridCellIndiciesInAabb(objectMin, objectMax, &gridCellIndicies);
		
		for(int j = 0; j < gridCellIndicies.size(); ++j)
		{
			FluidGridIterator FI = grid->getGridCell( gridCellIndicies[j] );
			
			for( int n = FI.m_firstIndex; FluidGridIterator::isIndexValid(n, FI.m_lastIndex);
					 n = FluidGridIterator::getNextIndex(n, isLinkedList, fluid->getNextFluidIndicies()) )
			{
				const btVector3 &fluidPos = fluid->getPosition(n);
				
				btVector3 fluidMin( fluidPos.x() - particleRadius, fluidPos.y() - particleRadius, fluidPos.z() - particleRadius );
				btVector3 fluidMax( fluidPos.x() + particleRadius, fluidPos.y() + particleRadius, fluidPos.z() + particleRadius );
						
				if( fluidMin.x() <= objectMax.x() && objectMin.x() <= fluidMax.x()  
				 && fluidMin.y() <= objectMax.y() && objectMin.y() <= fluidMax.y()
				 && fluidMin.z() <= objectMax.z() && objectMin.z() <= fluidMax.z() )
				{
					particleTransform.setOrigin(fluidPos);
					
					ParticleResult result(&particleObject);
					world->contactPairTest(&particleObject, const_cast<btCollisionObject*>(object), result);
					
					if( result.hasHit() ) 
					{
						resolveCollision(FG, fluid, n, result.m_collidedWith, 
										 result.m_normal, result.m_collidedWithHitPointWorld, result.m_distance);
					}
				}
			}
		}
	}
}

void BulletFluidsInterface::collideFluidsWithBulletCcd(const FluidParametersGlobal &FG, FluidSph *fluid, btCollisionWorld *world)
{
	BT_PROFILE("BulletFluidsInterface::collideFluidsWithBulletCcd()");
	
	//from collideFluidsWithBullet()
	//sph_pradius is at simscale; divide by simscale to transform into world scale
	btSphereShape particleShape( FG.sph_pradius / FG.sph_simscale );
	
	btCollisionObject particleObject;
	particleObject.setCollisionShape(&particleShape);
		
	btTransform &particleTransform = particleObject.getWorldTransform();
	particleTransform.setRotation( btQuaternion::getIdentity() );
	//from collideFluidsWithBullet()
	
	//int numCollided = 0;
	for(int i = 0; i < fluid->numParticles(); ++i)
	{
		const btVector3& pos = fluid->getPosition(i);
		const btVector3& prev_pos = fluid->getPrevPosition(i);
		
		//	investigate:
		//	Calling btCollisionWorld::rayTest() with the ray's origin inside a btCollisionObject
		//	reports no collisions; is this expected behavior when using btCollisionWorld::ClosestRayResultCallback?
		//btCollisionWorld::ClosestRayResultCallback result( btVector3(0.f, -1.0f, 0.f), btVector3(0.f, -2.0f, 0.f) );
		btCollisionWorld::ClosestRayResultCallback result(prev_pos, pos);
		result.m_collisionFilterGroup = btBroadphaseProxy::DefaultFilter;
		result.m_collisionFilterMask = btBroadphaseProxy::AllFilter;
		
		world->rayTest(result.m_rayFromWorld, result.m_rayToWorld, result);
		
		//printf( "result.m_closestHitFraction, result.hasHit(): %f, %d \n", result.m_closestHitFraction, result.hasHit() );
		
		if( result.hasHit() )
		{
			btScalar distanceMoved = result.m_rayFromWorld.distance(result.m_rayToWorld);
			btScalar distanceCollided = result.m_rayFromWorld.distance(result.m_hitPointWorld);
	
			btScalar distance = distanceMoved - distanceCollided;
		
			//++numCollided;
			resolveCollision(FG, fluid, i, result.m_collisionObject, result.m_hitNormalWorld, result.m_hitPointWorld, distance);
		}
		else
		{
			//Since btCollisionWorld::rayTest() does not report a collision if the ray's origin
			//is inside a btCollsionObject, use btCollisionWorld::contactTest() in case
			//the particle is inside.
			
			//from collideFluidsWithBullet()
			particleTransform.setOrigin(pos);
				
			ParticleResultMulti result(&particleObject);
			world->contactTest(&particleObject, result);
				
			for(int n = 0; n < result.m_numCollisions; ++n)
				resolveCollision(FG, fluid, i, result.m_collidedWith[n], 
								 result.m_normals[n], result.m_collidedWithHitPointsWorld[n], result.m_distances[n]);
			//from collideFluidsWithBullet()
		}
	}
	
	//printf("numCollided: %d \n", numCollided);
}	


void BulletFluidsInterface::resolveCollision(const FluidParametersGlobal &FG, FluidSph *fluid, int fluidIndex, btCollisionObject *object, 
											 const btVector3 &fluidNormal, const btVector3 &hitPointWorld, btScalar distance)
{
	const btScalar COLLISION_EPSILON = 0.00001f;

	const FluidParametersLocal &FL = fluid->getLocalParameters();
	
	btScalar depthOfPenetration = btFabs(distance)*FG.sph_simscale;
	if(depthOfPenetration > COLLISION_EPSILON)
	{
		btScalar adj = FL.m_extstiff * depthOfPenetration - FL.m_extdamp * fluidNormal.dot( fluid->getEvalVelocity(fluidIndex) );
		
		btVector3 acceleration = fluidNormal;
		acceleration *= adj;
			
		//if acceleration is very high, the fluid simulation will explode
		fluid->applyAcceleration(fluidIndex, acceleration);
		
		//btRigidBody-Fluid interaction
		//	test
		btRigidBody *rigidBody = btRigidBody::upcast(object);
		if( rigidBody && rigidBody->getInvMass() != 0.0f )
		{
			const btScalar fluidMass = FL.m_particleMass;
			const btScalar invertedMass = rigidBody->getInvMass();
			
			//Rigid body to particle
			//const btVector3& linearVelocity = rigidBody->getLinearVelocity();
			
			//Particle to rigid body
			btVector3 force = -acceleration * fluidMass;
			rigidBody->applyForce( force, hitPointWorld - rigidBody->getWorldTransform().getOrigin() );
			rigidBody->activate(true);
		}
	}
}

