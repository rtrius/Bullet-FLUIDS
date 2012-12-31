/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#include "btFluidSphRigidConstraintSolver.h"

#include "LinearMath/btVector3.h"

#include "BulletDynamics/Dynamics/btRigidBody.h"
#include "btFluidSph.h"

void btFluidSphRigidConstraintSolver::resolveCollisionsForce(const btFluidSphParametersGlobal& FG, btFluidSph *fluid)
{
	BT_PROFILE("resolveCollisionsForce()");
	
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	const btAlignedObjectArray<btFluidSphRigidContactGroup>& contactGroups = fluid->internalGetRigidContacts();
	
	m_accumulatedRigidForces.resize( contactGroups.size() );
	m_accumulatedRigidTorques.resize( contactGroups.size() );
	
	for(int i = 0; i < m_accumulatedRigidForces.size(); ++i) m_accumulatedRigidForces[i].setValue(0,0,0);
	for(int i = 0; i < m_accumulatedRigidTorques.size(); ++i) m_accumulatedRigidTorques[i].setValue(0,0,0);
	
	if(FL.m_enableAabbBoundary)
		btFluidSphRigidConstraintSolver::applyBoundaryForcesSingleFluid(FG, fluid);
	
	//Accumulate forces on rigid bodies, apply forces to fluid particles
	for(int i = 0; i < contactGroups.size(); ++i)
	{
		const btFluidSphRigidContactGroup& current = contactGroups[i];
	
		for(int n = 0; n < current.numContacts(); ++n)
		{
			resolveContactPenaltyForce(FG, fluid, const_cast<btCollisionObject*>(current.m_object),
										 current.m_contacts[n], m_accumulatedRigidForces[i], m_accumulatedRigidTorques[i]);
		}
	}
	
	//Apply forces to rigid bodies
	const btScalar timeStep = FG.m_timeStep;
	for(int i = 0; i < contactGroups.size(); ++i)
	{
		btRigidBody* rigidBody = btRigidBody::upcast( const_cast<btCollisionObject*>(contactGroups[i].m_object) );
		if( rigidBody && rigidBody->getInvMass() != btScalar(0.0) )
		{
			rigidBody->activate(true);
			
			btVector3 linearVelocity = rigidBody->getLinearVelocity();
			btVector3 angularVelocity = rigidBody->getAngularVelocity();
			
			linearVelocity += m_accumulatedRigidForces[i] * (rigidBody->getInvMass() * timeStep);
			angularVelocity += rigidBody->getInvInertiaTensorWorld() * m_accumulatedRigidTorques[i] * timeStep;
			
			const btScalar MAX_ANGVEL = SIMD_HALF_PI;
			btScalar angVel = angularVelocity.length();
			if(angVel*timeStep > MAX_ANGVEL) angularVelocity *= (MAX_ANGVEL/timeStep) / angVel;
			
			rigidBody->setLinearVelocity(linearVelocity);
			rigidBody->setAngularVelocity(angularVelocity);
		}
	}
}
void btFluidSphRigidConstraintSolver::resolveCollisionsImpulse(const btFluidSphParametersGlobal& FG, btFluidSph *fluid)
{
	BT_PROFILE("resolveCollisionsImpulse()");
	
	//Apply fluid particle-static impulses last to prevent penetration into static geometry
	const bool SEPARATE_STATIC_AND_DYNAMIC_RESPONSE = true;
	
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	const btAlignedObjectArray<btFluidSphRigidContactGroup>& contactGroups = fluid->internalGetRigidContacts();
	
	m_accumulatedFluidImpulses.resize( fluid->numParticles() );
	m_accumulatedRigidForces.resize( contactGroups.size() );
	m_accumulatedRigidTorques.resize( contactGroups.size() );
	
	for(int i = 0; i < m_accumulatedRigidForces.size(); ++i) m_accumulatedRigidForces[i].setValue(0,0,0);
	for(int i = 0; i < m_accumulatedRigidTorques.size(); ++i) m_accumulatedRigidTorques[i].setValue(0,0,0);
	for(int i = 0; i < fluid->numParticles(); ++i) m_accumulatedFluidImpulses[i].setValue(0,0,0);
	
	if(!SEPARATE_STATIC_AND_DYNAMIC_RESPONSE && FL.m_enableAabbBoundary) 
		btFluidSphRigidConstraintSolver::accumulateBoundaryImpulsesSingleFluid(FG, fluid, m_accumulatedFluidImpulses);
	
	//Accumulate forces on rigid bodies, impulses on fluids
	for(int i = 0; i < contactGroups.size(); ++i)
	{
		const btFluidSphRigidContactGroup& current = contactGroups[i];
		
		if(SEPARATE_STATIC_AND_DYNAMIC_RESPONSE)
		{
			const btRigidBody* rigidBody = btRigidBody::upcast(current.m_object);
			bool isDynamicRigidBody = ( rigidBody && rigidBody->getInvMass() != btScalar(0.0) );
			if(!isDynamicRigidBody) continue;
		}
		
		for(int n = 0; n < current.numContacts(); ++n)
		{
			resolveContactImpulse(FG, fluid, const_cast<btCollisionObject*>(current.m_object), current.m_contacts[n], 
									m_accumulatedRigidForces[i], m_accumulatedRigidTorques[i], m_accumulatedFluidImpulses);
		}
	}
	
	//Apply forces to rigid bodies
	const btScalar timeStep = FG.m_timeStep;
	for(int i = 0; i < contactGroups.size(); ++i)
	{
		btRigidBody* rigidBody = btRigidBody::upcast( const_cast<btCollisionObject*>(contactGroups[i].m_object) );
		if( rigidBody && rigidBody->getInvMass() != btScalar(0.0) )
		{
			rigidBody->activate(false);
			
			btVector3 linearVelocity = rigidBody->getLinearVelocity();
			btVector3 angularVelocity = rigidBody->getAngularVelocity();
			
			linearVelocity += m_accumulatedRigidForces[i] * (rigidBody->getInvMass() * timeStep);
			angularVelocity += rigidBody->getInvInertiaTensorWorld() * m_accumulatedRigidTorques[i] * timeStep;
			
			const btScalar MAX_ANGVEL = SIMD_HALF_PI;
			btScalar angVel = angularVelocity.length();
			if(angVel*timeStep > MAX_ANGVEL) angularVelocity *= (MAX_ANGVEL/timeStep) / angVel;
			
			rigidBody->setLinearVelocity(linearVelocity);
			rigidBody->setAngularVelocity(angularVelocity);
		}
	}
	
	//Apply impulses to fluid particles
	btFluidParticles& particles = fluid->internalGetParticles();
	
	for(int i = 0; i < fluid->numParticles(); ++i)
	{
		btVector3& vel = particles.m_vel[i];
		btVector3& vel_eval = particles.m_vel_eval[i];
		
		//Leapfrog integration
		btVector3 vnext = vel + m_accumulatedFluidImpulses[i];
		vel_eval = (vel + vnext) * btScalar(0.5);
		vel = vnext;
	}
	
	if(SEPARATE_STATIC_AND_DYNAMIC_RESPONSE)
	{
		for(int i = 0; i < fluid->numParticles(); ++i) m_accumulatedFluidImpulses[i].setValue(0,0,0);
		
		//Accumulate impulses on fluid particles
		for(int i = 0; i < contactGroups.size(); ++i)
		{
			const btFluidSphRigidContactGroup& current = contactGroups[i];

			const btRigidBody* rigidBody = btRigidBody::upcast(current.m_object);
			bool isDynamicRigidBody = ( rigidBody && rigidBody->getInvMass() != btScalar(0.0) );
			if(isDynamicRigidBody) continue;
			
			for(int n = 0; n < current.numContacts(); ++n)
			{
				resolveContactImpulse(FG, fluid, const_cast<btCollisionObject*>(current.m_object), current.m_contacts[n], 
										m_accumulatedRigidForces[i], m_accumulatedRigidTorques[i], m_accumulatedFluidImpulses);
			}
		}
		
		if(FL.m_enableAabbBoundary) 
			btFluidSphRigidConstraintSolver::accumulateBoundaryImpulsesSingleFluid(FG, fluid, m_accumulatedFluidImpulses);
		
		//Apply impulses to fluid particles
		for(int i = 0; i < fluid->numParticles(); ++i)
		{
			btVector3& vel = particles.m_vel[i];
			btVector3& vel_eval = particles.m_vel_eval[i];
			
			//Leapfrog integration
			btVector3 vnext = vel + m_accumulatedFluidImpulses[i];
			vel_eval = (vel + vnext) * btScalar(0.5);
			vel = vnext;
		}
	}
}
void btFluidSphRigidConstraintSolver::resolveContactPenaltyForce(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, 
																btCollisionObject *object, const btFluidSphRigidContact& contact,
																btVector3 &accumulatedRigidForce, btVector3 &accumulatedRigidTorque)
{
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	
	if( contact.m_distance < btScalar(0.0) )
	{
		btRigidBody* rigidBody = btRigidBody::upcast(object);
		bool isDynamicRigidBody = ( rigidBody && rigidBody->getInvMass() != btScalar(0.0) );
		
		btVector3 rigidLocalHitPoint = contact.m_hitPointWorldOnObject - object->getWorldTransform().getOrigin();
		
		btVector3 fluidVelocity = fluid->getEvalVelocity(contact.m_fluidParticleIndex);
		btVector3 rigidVelocity = (isDynamicRigidBody) ? rigidBody->getVelocityInLocalPoint(rigidLocalHitPoint) : btVector3(0,0,0); 
		rigidVelocity *= FG.m_simulationScale;
	
		btVector3 relativeVelocity = fluidVelocity - rigidVelocity;
		btScalar relativeNormalMagnitude = contact.m_normalOnObject.dot(relativeVelocity);
		
		btScalar penetrationDepth = -contact.m_distance * FG.m_simulationScale;
		btScalar forceMagnitude = FL.m_boundaryStiff * penetrationDepth - FL.m_boundaryDamp * relativeNormalMagnitude;
		
		btVector3 acceleration = contact.m_normalOnObject * forceMagnitude;
		
		if( FL.m_boundaryFriction != btScalar(0.0) )
		{
			//Since an acceleration is used to resolve the collision,
			//the change in velocity per frame is (acceleration * FG.timeStep).
			const btScalar tangentRemovedPerFrame = FL.m_boundaryFriction / FG.m_timeStep; 
			
			btVector3 relativeNormalVelocity = contact.m_normalOnObject * relativeNormalMagnitude;
			btVector3 relativeTangentialVelocity = relativeVelocity - relativeNormalVelocity;
			acceleration -= relativeTangentialVelocity * tangentRemovedPerFrame;
		}
		
		btVector3 force = acceleration * FL.m_particleMass;
		
		if(isDynamicRigidBody)
		{
			btVector3 worldScaleForce = -force / FG.m_simulationScale;
		
			const btVector3& linearFactor = rigidBody->getLinearFactor();
			accumulatedRigidForce += worldScaleForce * linearFactor;
			accumulatedRigidTorque += rigidLocalHitPoint.cross(worldScaleForce * linearFactor) * rigidBody->getAngularFactor();
		}
		
		//if acceleration is very high, the fluid simulation will explode
		fluid->applyForce(contact.m_fluidParticleIndex, force);
	}
}

void btFluidSphRigidConstraintSolver::resolveContactImpulse(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, 
																btCollisionObject *object, const btFluidSphRigidContact& contact,
																btVector3 &accumulatedRigidForce, btVector3 &accumulatedRigidTorque,
																btAlignedObjectArray<btVector3>& accumulatedImpulses)
{
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	
	if( contact.m_distance < btScalar(0.0) )
	{
		int i = contact.m_fluidParticleIndex;
		btFluidParticles& particles = fluid->internalGetParticles();
		
		btRigidBody* rigidBody = btRigidBody::upcast(object);
		bool isDynamicRigidBody = ( rigidBody && rigidBody->getInvMass() != btScalar(0.0) );
		
		btVector3 rigidLocalHitPoint = contact.m_hitPointWorldOnObject - object->getWorldTransform().getOrigin();
		
		btVector3 fluidVelocity = fluid->getVelocity(i);
		btVector3 rigidVelocity = (isDynamicRigidBody) ? rigidBody->getVelocityInLocalPoint(rigidLocalHitPoint) : btVector3(0,0,0); 
		rigidVelocity *= FG.m_simulationScale;
	
		btVector3 relativeVelocity = fluidVelocity - rigidVelocity;
		btScalar penetratingMagnitude = relativeVelocity.dot(-contact.m_normalOnObject);
		if( penetratingMagnitude < btScalar(0.0) ) penetratingMagnitude = btScalar(0.0);
		
		btVector3 penetratingVelocity = -contact.m_normalOnObject * penetratingMagnitude;
		btVector3 tangentialVelocity = relativeVelocity - penetratingVelocity;
		
		penetratingVelocity *= btScalar(1.0) + FL.m_boundaryRestitution;
		
		btScalar positionError = (-contact.m_distance) * (FG.m_simulationScale/FG.m_timeStep) * FL.m_boundaryErp;
		btVector3 particleImpulse = -(penetratingVelocity + (-contact.m_normalOnObject*positionError) + tangentialVelocity*FL.m_boundaryFriction);
		
		if(isDynamicRigidBody)
		{
			btScalar inertiaParticle = btScalar(1.0) / FL.m_particleMass;
			
			btVector3 relPosCrossNormal = rigidLocalHitPoint.cross(contact.m_normalOnObject);
			btScalar inertiaRigid = rigidBody->getInvMass() + ( relPosCrossNormal * rigidBody->getInvInertiaTensorWorld() ).dot(relPosCrossNormal);
			
			particleImpulse *= btScalar(1.0) / (inertiaParticle + inertiaRigid);
			
			btVector3 worldScaleImpulse = -particleImpulse / FG.m_simulationScale;
			worldScaleImpulse /= FG.m_timeStep;		//Impulse is accumulated as force
			
			const btVector3& linearFactor = rigidBody->getLinearFactor();
			accumulatedRigidForce += worldScaleImpulse * linearFactor;
			accumulatedRigidTorque += rigidLocalHitPoint.cross(worldScaleImpulse * linearFactor) * rigidBody->getAngularFactor();
			
			particleImpulse *= inertiaParticle;
		}
		
		accumulatedImpulses[i] += particleImpulse;
	}
}


inline void resolveAabbCollision(const btFluidSphParametersLocal& FL, const btVector3& vel_eval,
								 btVector3* acceleration, const btVector3& normal, btScalar distance)
{
	if( distance < btScalar(0.0) )	//Negative distance indicates penetration
	{
		btScalar penetrationDepth = -distance;
	
		btScalar accelerationMagnitude = FL.m_boundaryStiff * penetrationDepth - FL.m_boundaryDamp * normal.dot(vel_eval);
		
		*acceleration += normal * accelerationMagnitude;
	}
}
void applyBoundaryForceToParticle(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, btScalar simScaleParticleRadius,
									btFluidParticles& particles, int particleIndex)
{
	int i = particleIndex;
	
	const btScalar radius = simScaleParticleRadius;
	const btScalar simScale = FG.m_simulationScale;
	
	const btVector3& min = FL.m_aabbBoundaryMin;
	const btVector3& max = FL.m_aabbBoundaryMax;
	
	const btVector3& pos = particles.m_pos[i];
	const btVector3& vel_eval = particles.m_vel_eval[i];
	
	btVector3 acceleration(0,0,0);
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3( 1.0, 0.0, 0.0), ( pos.x() - min.x() )*simScale - radius );
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3(-1.0, 0.0, 0.0), ( max.x() - pos.x() )*simScale - radius );
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3(0.0,  1.0, 0.0), ( pos.y() - min.y() )*simScale - radius );
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3(0.0, -1.0, 0.0), ( max.y() - pos.y() )*simScale - radius );
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3(0.0, 0.0,  1.0), ( pos.z() - min.z() )*simScale - radius );
	resolveAabbCollision( FL, vel_eval, &acceleration, btVector3(0.0, 0.0, -1.0), ( max.z() - pos.z() )*simScale - radius );
	
	particles.m_accumulatedForce[i] += acceleration * FL.m_particleMass;
}
void btFluidSphRigidConstraintSolver::applyBoundaryForcesSingleFluid(const btFluidSphParametersGlobal& FG, btFluidSph* fluid)
{
	BT_PROFILE("applyBoundaryForcesSingleFluid()");
	
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	btFluidParticles& particles = fluid->internalGetParticles();
	
	const btScalar simScaleParticleRadius = FL.m_particleRadius * FG.m_simulationScale;

	for(int i = 0; i < particles.size(); ++i) applyBoundaryForceToParticle(FG, FL, simScaleParticleRadius, particles, i);
}


inline void resolveAabbCollision_impulse(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, const btVector3& velocity, 
										const btVector3& normal, btScalar distance, btVector3& out_impulse)
{
	if( distance < btScalar(0.0) )	//Negative distance indicates penetration
	{
		btScalar penetratingMagnitude = velocity.dot(-normal);
		if( penetratingMagnitude < btScalar(0.0) ) penetratingMagnitude = btScalar(0.0);
		
		btVector3 penetratingVelocity = -normal * penetratingMagnitude;
		btVector3 tangentialVelocity = velocity - penetratingVelocity;
		
		penetratingVelocity *= btScalar(1.0) + FL.m_boundaryRestitution;
		
		btScalar positionError = (-distance) * (FG.m_simulationScale/FG.m_timeStep) * FL.m_boundaryErp;
		
		out_impulse += -( penetratingVelocity + (-normal*positionError) + tangentialVelocity * FL.m_boundaryFriction );
	}
}
void applyBoundaryImpulseToParticle(const btFluidSphParametersGlobal& FG, btScalar simScaleParticleRadius,
									const btFluidSphParametersLocal& FL, btFluidParticles& particles, int particleIndex,
									btVector3& out_accumulatedImpulse)
{
	int i = particleIndex;
	
	const btScalar radius = simScaleParticleRadius;
	const btScalar simScale = FG.m_simulationScale;
	
	const btVector3& boundaryMin = FL.m_aabbBoundaryMin;
	const btVector3& boundaryMax = FL.m_aabbBoundaryMax;
	
	const btVector3& pos = particles.m_pos[i];
	const btVector3& vel = particles.m_vel[i];
	
	btVector3& impulse = out_accumulatedImpulse;
	resolveAabbCollision_impulse( FG, FL, vel, btVector3( 1.0, 0.0, 0.0), ( pos.x() - boundaryMin.x() )*simScale - radius, impulse );
	resolveAabbCollision_impulse( FG, FL, vel, btVector3(-1.0, 0.0, 0.0), ( boundaryMax.x() - pos.x() )*simScale - radius, impulse );
	resolveAabbCollision_impulse( FG, FL, vel, btVector3(0.0,  1.0, 0.0), ( pos.y() - boundaryMin.y() )*simScale - radius, impulse );
	resolveAabbCollision_impulse( FG, FL, vel, btVector3(0.0, -1.0, 0.0), ( boundaryMax.y() - pos.y() )*simScale - radius, impulse );
	resolveAabbCollision_impulse( FG, FL, vel, btVector3(0.0, 0.0,  1.0), ( pos.z() - boundaryMin.z() )*simScale - radius, impulse );
	resolveAabbCollision_impulse( FG, FL, vel, btVector3(0.0, 0.0, -1.0), ( boundaryMax.z() - pos.z() )*simScale - radius, impulse );
}
void btFluidSphRigidConstraintSolver::accumulateBoundaryImpulsesSingleFluid(const btFluidSphParametersGlobal& FG, btFluidSph* fluid,
																		btAlignedObjectArray<btVector3>& accumulatedFluidImpulses)
{
	BT_PROFILE("accumulateBoundaryImpulsesSingleFluid()");
	
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	btFluidParticles& particles = fluid->internalGetParticles();
	
	const btScalar simScaleParticleRadius = FL.m_particleRadius * FG.m_simulationScale;
	
	for(int i = 0; i < particles.size(); ++i)
		applyBoundaryImpulseToParticle(FG, simScaleParticleRadius, FL, particles, i, accumulatedFluidImpulses[i]);
}



