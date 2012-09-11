/* 	FluidRigidConstraintSolver.h
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
#ifndef FLUID_RIGID_CONSTRAINT_SOLVER_H
#define FLUID_RIGID_CONSTRAINT_SOLVER_H

#include "BulletDynamics/Dynamics/btRigidBody.h"

#include "FluidRigidCollisionDetector.h"

///Resolves collisions between FluidSph and btCollisionObject/btRigidBody.
class FluidRigidConstraintSolver
{
public:
	void resolveCollisions(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidRigidContactGroup> *contactGroups)
	{
		BT_PROFILE("FluidRigidConstraintSolver::resolveCollisions()");
	
		for(int i = 0; i < contactGroups->size(); ++i)
			resolveCollisionsSingleFluid(FG, &(*contactGroups)[i]);
	}
	
private:
	void resolveCollisionsSingleFluid(const FluidParametersGlobal &FG, FluidRigidContactGroup *contactGroup)
	{
		FluidSph *fluid = contactGroup->m_fluid;
		for(int i = 0; i < contactGroup->m_contacts.size(); ++i)
			resolveCollisionPenaltyForce(FG, fluid, &contactGroup->m_contacts[i]);
	}
	
	void resolveCollisionPenaltyForce(const FluidParametersGlobal &FG, FluidSph *fluid, FluidRigidContact *contact);
};

#endif