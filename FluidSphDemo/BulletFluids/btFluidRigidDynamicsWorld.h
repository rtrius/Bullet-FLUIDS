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

#include "Sph/btFluidSph.h"
#include "Sph/btFluidSphSolver.h"
#include "Sph/btFluidSphRigidCollisionDetector.h"
#include "Sph/btFluidSphRigidConstraintSolver.h"

///@brief Coordinates several btFluidSph and global fluid properities.
///@remarks
///Terminology:
/// - World - a set of fluids.
/// - Fluid - a set of particles.
/// - Particle - a point with position and velocity that may be influenced by other particles using SPH.
class btFluidRigidDynamicsWorld : public btDiscreteDynamicsWorld
{
	btFluidSphParametersGlobal m_globalParameters;

	btAlignedObjectArray<btFluidSph*> m_fluids;
	btAlignedObjectArray<btFluidSph*> m_tempOverrideSolverFluids;	//Contains the subset of m_fluids with (getOverrideSolver() != 0)
	btAlignedObjectArray<btFluidSph*> m_tempDefaultSolverFluids;	//Contains the subset of m_fluids with (getOverrideSolver() == 0)
	
	btFluidSphSolver* m_fluidSolver;
	
	btFluidSphRigidCollisionDetector m_fluidRigidCollisionDetector;
	btFluidSphRigidConstraintSolver m_fluidRigidConstraintSolver;
	
	
public:
	btFluidRigidDynamicsWorld(btDispatcher* dispatcher, btBroadphaseInterface* pairCache, btConstraintSolver* constraintSolver, 
							btCollisionConfiguration* collisionConfiguration, btFluidSphSolver* fluidSolver) 
								: btDiscreteDynamicsWorld(dispatcher, pairCache, constraintSolver, collisionConfiguration),
								m_fluidSolver(fluidSolver) {}
	virtual ~btFluidRigidDynamicsWorld() {}
	
	virtual int stepSimulation( btScalar timeStep, int maxSubSteps = 1, btScalar fixedTimeStep = btScalar(1.0)/btScalar(60.0) )
	{
		m_tempOverrideSolverFluids.resize(0);
		m_tempDefaultSolverFluids.resize(0);
		for(int i = 0; i < m_fluids.size(); ++i)
		{
			if( m_fluids[i]->getOverrideSolver() ) m_tempOverrideSolverFluids.push_back(m_fluids[i]);
			else m_tempDefaultSolverFluids.push_back(m_fluids[i]);
		}
	
		//
		return btDiscreteDynamicsWorld::stepSimulation(timeStep, maxSubSteps, fixedTimeStep);
	}
	
	const btFluidSphParametersGlobal& getGlobalParameters() const { return m_globalParameters; }
	void setGlobalParameters(const btFluidSphParametersGlobal& FG) { m_globalParameters = FG; }
	
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
	
	
	void setFluidSolver(btFluidSphSolver* solver) { m_fluidSolver = solver; }
	btFluidSphSolver* getFluidSolver() const { return m_fluidSolver; }
	
	btAlignedObjectArray<btFluidSph*>& internalGetFluids() { return m_fluids; }
	
	//virtual btDynamicsWorldType getWorldType() const { return BT_FLUID_RIGID_DYNAMICS_WORLD; }
	
protected:
	virtual void internalSingleStepSimulation(btScalar timeStep) 
	{
		BT_PROFILE("FluidRigidWorld - singleStepSimulation()");
	
		for(int i = 0; i < m_fluids.size(); ++i) m_fluids[i]->internalClearRigidContacts();
		
		//btFluidSph-btRigidBody/btCollisionObject AABB intersections are 
		//detected here(not midphase/narrowphase), so calling removeMarkedParticles() 
		//below does not invalidate the collisions.
		btDiscreteDynamicsWorld::internalSingleStepSimulation(timeStep);
		
		for(int i = 0; i < m_fluids.size(); ++i) m_fluids[i]->removeMarkedParticles();
		
		//SPH forces
		{
			int numDefaultSolverFluids = m_tempDefaultSolverFluids.size();
			if(numDefaultSolverFluids)
			{
				m_fluidSolver->updateGridAndCalculateSphForces(m_globalParameters, &m_tempDefaultSolverFluids[0], numDefaultSolverFluids);
			}
			
			for(int i = 0; i < m_tempOverrideSolverFluids.size(); ++i)
			{
				btFluidSphSolver* overrideSolver = m_tempOverrideSolverFluids[i]->getOverrideSolver();
				
				overrideSolver->updateGridAndCalculateSphForces( m_globalParameters, &m_tempOverrideSolverFluids[i], 1 );
			}
		}
		
		for(int i = 0; i < m_fluids.size(); ++i)
			m_fluidRigidCollisionDetector.performNarrowphase(m_dispatcher1, m_dispatchInfo, m_fluids[i]);
		
		//Resolve collisions, integrate
		{
			const bool USE_IMPULSE_BOUNDARY = 1; 	//Penalty force otherwise
		
			for(int i = 0; i < m_fluids.size(); ++i) 
			{
				btFluidSph* fluid = m_fluids[i];
				const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
				
				if(!USE_IMPULSE_BOUNDARY)
				{
					m_fluidRigidConstraintSolver.resolveParticleCollisions(m_globalParameters, m_fluids[i], USE_IMPULSE_BOUNDARY);
				}
				btFluidSphSolver::applyForcesSingleFluid(m_globalParameters, fluid);
				
				if(USE_IMPULSE_BOUNDARY)
				{
					m_fluidRigidConstraintSolver.resolveParticleCollisions(m_globalParameters, m_fluids[i], USE_IMPULSE_BOUNDARY);
				}
				btFluidSphSolver::integratePositionsSingleFluid( m_globalParameters, fluid->internalGetParticles() );
			}
		}
	}
};


#endif
