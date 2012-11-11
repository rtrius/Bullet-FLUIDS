/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. 
   If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_FLUID_RIGID_DYNAMICS_WORLD_H
#define BT_FLUID_RIGID_DYNAMICS_WORLD_H

#include "BulletDynamics/Dynamics/btDiscreteDynamicsWorld.h"

#include "btFluidSph.h"
#include "btFluidSolver.h"
#include "btFluidRigidCollisionDetector.h"
#include "btFluidRigidConstraintSolver.h"

///@brief Coordinates several btFluidSph and global fluid properities.
///@remarks
///Terminology:
/// - World - a set of fluids.
/// - Fluid - a set of particles.
/// - Particle - a point with position and velocity that may be influenced by other particles using SPH.
class btFluidRigidDynamicsWorld : public btDiscreteDynamicsWorld
{
	btFluidParametersGlobal m_globalParameters;

	btAlignedObjectArray<btFluidSph*> m_fluids;

	btFluidSolver* m_fluidSolver;
	
	btFluidRigidCollisionDetector m_fluidRigidCollisionDetector;
	btFluidRigidConstraintSolver m_fluidRigidConstraintSolver;
	
public:
	btFluidRigidDynamicsWorld(btDispatcher* dispatcher, btBroadphaseInterface* pairCache, btConstraintSolver* constraintSolver, 
							btCollisionConfiguration* collisionConfiguration, btFluidSolver* fluidSolver) 
								: btDiscreteDynamicsWorld(dispatcher, pairCache, constraintSolver, collisionConfiguration),
								m_fluidSolver(fluidSolver) {}
	virtual ~btFluidRigidDynamicsWorld() {}
	
	const btFluidParametersGlobal& getGlobalParameters() const { return m_globalParameters; }
	void setGlobalParameters(const btFluidParametersGlobal& FG) { m_globalParameters = FG; }
	
	void addFluid(btFluidSph* fluid, short int collisionFilterGroup = btBroadphaseProxy::DefaultFilter,
									short int collisionFilterMask = btBroadphaseProxy::AllFilter)
	{
		btAssert(fluid);
		btAssert( m_fluids.findLinearSearch(fluid) == m_fluids.size() );
		
		m_fluids.push_back(fluid);
		btCollisionWorld::addCollisionObject(fluid, collisionFilterGroup, collisionFilterMask);
	}
	void removeFluid(btFluidSph* fluid)
	{
		m_fluids.remove(fluid);  //Swaps elements if fluid is not the last
		btCollisionWorld::removeCollisionObject(fluid);
	}
	
	virtual void addCollisionObject(btCollisionObject* collisionObject, short int collisionFilterGroup = btBroadphaseProxy::DefaultFilter,
																		short int collisionFilterMask = btBroadphaseProxy::AllFilter)
	{
		btFluidSph* fluid = btFluidSph::upcast(collisionObject);
		if(fluid) addFluid(fluid, collisionFilterGroup, collisionFilterMask);
		else btDiscreteDynamicsWorld::addCollisionObject(collisionObject, collisionFilterGroup, collisionFilterMask);
	}

	virtual void removeCollisionObject(btCollisionObject* collisionObject)
	{
		btFluidSph* fluid = btFluidSph::upcast(collisionObject);
		if(fluid) removeFluid(fluid);
		else btDiscreteDynamicsWorld::removeCollisionObject(collisionObject);
	}
	

	int getNumFluids() const { return m_fluids.size(); }
	btFluidSph* getFluid(int index) { return m_fluids[index]; }
	
	
	void setFluidSolver(btFluidSolver* solver) { m_fluidSolver = solver; }
	btFluidSolver* getFluidSolver() { return m_fluidSolver; }
	
	btAlignedObjectArray<btFluidSph*>& internalGetFluids() { return m_fluids; }
	
	//virtual btDynamicsWorldType getWorldType() const { return BT_FLUID_RIGID_DYNAMICS_WORLD; }
protected:
	virtual void internalSingleStepSimulation(btScalar timeStep) 
	{
		for(int i = 0; i < m_fluids.size(); ++i) m_fluids[i]->internalClearRigidContacts();
		
		btDiscreteDynamicsWorld::internalSingleStepSimulation(timeStep);
		
		for(int i = 0; i < m_fluids.size(); ++i) m_fluids[i]->removeMarkedParticles();
		
		if( m_fluids.size() )m_fluidSolver->stepSimulation( m_globalParameters, &m_fluids[0], m_fluids.size() ); 
		
		//	fix: AABB/broadphase may not be syncronized when determining fluid-rigid interaction
		{
			BT_PROFILE("FluidRigid Interaction");
			
			for(int i = 0; i < m_fluids.size(); ++i)
				m_fluidRigidCollisionDetector.detectCollisionsSingleFluid(m_dispatcher1, m_dispatchInfo, m_fluids[i]);
				
			for(int i = 0; i < m_fluids.size(); ++i)
				m_fluidRigidConstraintSolver.resolveCollisionsSingleFluid(m_globalParameters, m_fluids[i]);
		}
	}
};


#endif
