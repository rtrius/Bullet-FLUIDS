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
#ifndef BT_FLUID_SPH_RIGID_CONSTRAINT_SOLVER_H
#define BT_FLUID_SPH_RIGID_CONSTRAINT_SOLVER_H

#include "LinearMath/btAlignedObjectArray.h"

class btVector3;
class btCollisionObject;
struct btFluidSphParametersGlobal;
struct btFluidSphRigidContact;
class btFluidSph;

///Resolves collisions between btFluidSph and btCollisionObject / btRigidBody.
class btFluidSphRigidConstraintSolver
{
	btAlignedObjectArray<btVector3> m_accumulatedFluidImpulses;

	btAlignedObjectArray<btVector3> m_accumulatedRigidForces;	//Each element corresponds to a btCollisionObject / btRigidBody
	btAlignedObjectArray<btVector3> m_accumulatedRigidTorques;

public:
	void resolveParticleCollisions(const btFluidSphParametersGlobal& FG, btFluidSph *fluid, bool useImpulses);
	
	static void applyBoundaryForcesSingleFluid(const btFluidSphParametersGlobal& FG, btFluidSph* fluid);
	static void accumulateBoundaryImpulsesSingleFluid(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, 
												btAlignedObjectArray<btVector3>& accumulatedFluidImpulses);
	
private:
	void resolveContactPenaltyForce(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, 
									btCollisionObject *object, const btFluidSphRigidContact& contact,
									btVector3& accumulatedRigidForce, btVector3& accumulatedRigidTorque);
									
	void resolveContactImpulse(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, 
									btCollisionObject *object, const btFluidSphRigidContact& contact,
									btVector3& accumulatedRigidForce, btVector3& accumulatedRigidTorque,
									btAlignedObjectArray<btVector3>& accumulatedFluidImpulses);
};

#endif