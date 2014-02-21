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
#ifndef BT_FLUID_RIGID_DYNAMICS_WORLD_H
#define BT_FLUID_RIGID_DYNAMICS_WORLD_H

#include "BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h"

#include "Sph/btFluidSphParameters.h"
#include "Sph/btFluidSphRigidCollisionDetector.h"
#include "Sph/btFluidSphRigidConstraintSolver.h"

class btFluidSph;
class btFluidSphSolver;
class btFluidRigidDynamicsWorld;

class btFluidEmitter;
class btFluidAbsorber;

typedef void (*btInternalFluidTickCallback)(btFluidRigidDynamicsWorld* world, btScalar timeStep);

///Adds particle fluid simulation support on top of btDiscreteDynamicsWorld.
class btFluidRigidDynamicsWorld : public btDiscreteDynamicsWorld
{
protected:
	btFluidSphParametersGlobal m_globalParameters;

	btAlignedObjectArray<btFluidSph*> m_fluids;
	btAlignedObjectArray<btFluidSph*> m_tempOverrideFluids;	//Contains the subset of m_fluids with (getOverrideSolver/Parameters() != 0)
	btAlignedObjectArray<btFluidSph*> m_tempDefaultFluids;	//Contains the subset of m_fluids with no override set
	
	btFluidSphSolver* m_fluidSolver;
	
	btFluidSphRigidCollisionDetector m_fluidRigidCollisionDetector;
	btFluidSphRigidConstraintSolver m_fluidRigidConstraintSolver;
	
	btInternalFluidTickCallback m_internalFluidPreTickCallback;
	btInternalFluidTickCallback m_internalFluidPostTickCallback;
	btInternalFluidTickCallback m_internalFluidMidTickCallback;
	
	btAlignedObjectArray<btFluidEmitter*> m_emitters;
	//btAlignedObjectArray<btFluidAbsorber*> m_absorbers;
	
	
public:
	btFluidRigidDynamicsWorld(btDispatcher* dispatcher, btBroadphaseInterface* pairCache, btConstraintSolver* constraintSolver, 
							btCollisionConfiguration* collisionConfiguration, btFluidSphSolver* fluidSolver);
	virtual ~btFluidRigidDynamicsWorld() {}
	
	virtual int stepSimulation( btScalar timeStep, int maxSubSteps = 1, btScalar fixedTimeStep = btScalar(1.0/60.0) );
	
	//
	void addFluidSph(btFluidSph* fluid, short int collisionFilterGroup = btBroadphaseProxy::DefaultFilter,
									short int collisionFilterMask = btBroadphaseProxy::AllFilter);
	void removeFluidSph(btFluidSph* fluid);
	
	virtual void addCollisionObject(btCollisionObject* collisionObject, short int collisionFilterGroup = btBroadphaseProxy::DefaultFilter,
																		short int collisionFilterMask = btBroadphaseProxy::AllFilter);
	virtual void removeCollisionObject(btCollisionObject* collisionObject);
	
	void addSphEmitter(btFluidEmitter* emitter);
	void removeSphEmitter(btFluidEmitter* emitter);
	
	int getNumSphEmitters() const { return m_emitters.size(); }
	btFluidEmitter* getSphEmitter(int index) { return m_emitters[index]; }
	
	//
	int getNumFluidSph() const { return m_fluids.size(); }
	btFluidSph* getFluidSph(int index) { return m_fluids[index]; }
	
	const btFluidSphParametersGlobal& getGlobalParameters() const { return m_globalParameters; }
	void setGlobalParameters(const btFluidSphParametersGlobal& FG) { m_globalParameters = FG; }
	
	btFluidSphSolver* getFluidSolver() const { return m_fluidSolver; }
	void setFluidSolver(btFluidSphSolver* solver) { m_fluidSolver = solver; }
	
	btAlignedObjectArray<btFluidSph*>& internalGetFluids() { return m_fluids; }
	
	//virtual btDynamicsWorldType getWorldType() const { return BT_FLUID_RIGID_DYNAMICS_WORLD; }
	
	virtual void debugDrawWorld();
	
	///@param isPreTick Sets pre-tick callback if true, sets post-tick callback otherwise 
	void setInternalFluidTickCallback(btInternalFluidTickCallback cb, bool isPreTick = false) 
	{
		if(isPreTick) m_internalFluidPreTickCallback = cb;
		else m_internalFluidPostTickCallback = cb;
	}
	
	void setInternalFluidMidTickcallback(btInternalFluidTickCallback cb) { m_internalFluidMidTickCallback = cb; }
	
protected:
	virtual void internalSingleStepSimulation(btScalar timeStep) ;
};


#endif
