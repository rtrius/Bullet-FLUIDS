/* 	FluidRigidConstraintSolver.cpp
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
#include "FluidRigidConstraintSolver.h"

void FluidRigidConstraintSolver::resolveCollisionPenaltyForce(const FluidParametersGlobal &FG, FluidSph *fluid, FluidRigidContact *contact)
{
	const FluidParametersLocal &FL = fluid->getLocalParameters();
	
	if( contact->m_distance < btScalar(0.0) )
	{
		btRigidBody *rigidBody = btRigidBody::upcast(contact->m_object);
		bool isDynamicRigidBody = ( rigidBody && rigidBody->getInvMass() != btScalar(0.0) );
		
		btVector3 rigidLocalHitPoint = contact->m_hitPointWorldOnObject - contact->m_object->getWorldTransform().getOrigin();
		
		btVector3 fluidVelocity = fluid->getEvalVelocity(contact->m_fluidParticleIndex);
		btVector3 rigidVelocity = (isDynamicRigidBody) ? rigidBody->getVelocityInLocalPoint(rigidLocalHitPoint) : btVector3(0,0,0); 
		rigidVelocity *= FG.m_simulationScale;
	
		btVector3 relativeVelocity = fluidVelocity - rigidVelocity;
		
		btScalar depthOfPenetration = -contact->m_distance * FG.m_simulationScale;
		btScalar forceMagnitude = FL.m_extstiff * depthOfPenetration - FL.m_extdamp * contact->m_normalOnObject.dot(relativeVelocity);
		
		btVector3 acceleration = contact->m_normalOnObject * forceMagnitude;
		
		if(isDynamicRigidBody)
		{
			btVector3 force = -acceleration * (FL.m_particleMass / FG.m_simulationScale);
			rigidBody->applyForce(force, rigidLocalHitPoint);
			rigidBody->activate(true);
		}
		
		//if acceleration is very high, the fluid simulation will explode
		fluid->applyAcceleration(contact->m_fluidParticleIndex, acceleration);
	}
}
