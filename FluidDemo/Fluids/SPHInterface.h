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

#include "FluidSph.h"
#include "FluidWorld.h"

#include "btBulletDynamicsCommon.h"
#include "LinearMath/btAlignedObjectArray.h"

struct ParticleResult : public btCollisionWorld::ContactResultCallback
{
	btVector3 m_normal;
	btVector3 m_collidedWithHitPointWorld;
	btCollisionObject *m_collidedWith;
	btScalar m_distance;
	
	btCollisionObject *m_particleObject;

	ParticleResult(btCollisionObject *particleObject) : m_particleObject(particleObject), m_collidedWith(0) {}
	
	virtual btScalar addSingleResult( btManifoldPoint &cp, const btCollisionObject *colObj0, int partId0, int index0,
														   const btCollisionObject *colObj1, int partId1, int index1 ) 
	{
		m_distance = cp.getDistance();
	
		//Assume 0 == A, 1 == B
		if(m_particleObject == colObj0 && colObj1)		
		{
			m_collidedWith = const_cast<btCollisionObject*>(colObj1);
			m_normal = cp.m_normalWorldOnB;
			m_collidedWithHitPointWorld = cp.getPositionWorldOnB();
		}
		else if(m_particleObject == colObj1 && colObj0)
		{
			//	this branch is never reached?
		
			m_collidedWith = const_cast<btCollisionObject*>(colObj0);
			m_normal = -cp.m_normalWorldOnB;
			m_collidedWithHitPointWorld = cp.getPositionWorldOnA();
		}
		
		//Value returned from btCollisionWorld::ContactResultCallback::addSingleResult() appears to be unused
		const btScalar UNUSED = btScalar(1.0);
		return UNUSED;
	}

	const bool hasHit() const { return static_cast<bool>(m_collidedWith); }
};

struct ParticleResultMulti : public btCollisionWorld::ContactResultCallback
{
	static const int MAX_COLLISIONS = 4;

	btVector3 m_normals[MAX_COLLISIONS];
	btVector3 m_collidedWithHitPointsWorld[MAX_COLLISIONS];
	btCollisionObject *m_collidedWith[MAX_COLLISIONS];
	btScalar m_distances[MAX_COLLISIONS];
	int m_numCollisions;
	
	btCollisionObject *m_particleObject;

	ParticleResultMulti(btCollisionObject *particleObject) : m_particleObject(particleObject), m_numCollisions(0) {}
	
	virtual btScalar addSingleResult( btManifoldPoint &cp, const btCollisionObject *colObj0, int partId0, int index0,
														   const btCollisionObject *colObj1, int partId1, int index1 ) 
	{
		//Value returned from btCollisionWorld::ContactResultCallback::addSingleResult() appears to be unused
		const btScalar UNUSED = 1.0f;
	
		if(m_numCollisions >= MAX_COLLISIONS) return UNUSED;
	
		//Assume 0 == A, 1 == B
		if(m_particleObject == colObj0 && colObj1)
		{
			m_collidedWith[m_numCollisions] = const_cast<btCollisionObject*>(colObj1);
			m_normals[m_numCollisions] = cp.m_normalWorldOnB;
			m_collidedWithHitPointsWorld[m_numCollisions] = cp.getPositionWorldOnB();
			m_distances[m_numCollisions] = cp.getDistance();
			++m_numCollisions;
		}
		else if(m_particleObject == colObj1 && colObj0)
		{
			//	this branch is never reached?
			
			m_collidedWith[m_numCollisions] = const_cast<btCollisionObject*>(colObj0);
			m_normals[m_numCollisions] = -cp.m_normalWorldOnB;
			m_collidedWithHitPointsWorld[m_numCollisions] = cp.getPositionWorldOnA();
			m_distances[m_numCollisions] = cp.getDistance();
			++m_numCollisions;
		}
		
		return UNUSED;
	}
};

struct AabbTestCallback : public btBroadphaseAabbCallback
{
	btAlignedObjectArray<btCollisionObject*> m_results;

	virtual bool process(const btBroadphaseProxy *proxy) 
	{ 
		m_results.push_back( static_cast<btCollisionObject*>(proxy->m_clientObject) ); 
		
		//	investigate meaning of return value
		return true;
	}
};
class BulletFluidsInterface
{
	static void collideFluidsWithBullet(const FluidParametersGlobal &FG, FluidSph *fluid, btCollisionWorld *world);
	static void collideFluidsWithBullet2(const FluidParametersGlobal &FG, FluidSph *fluid, btCollisionWorld *world);
	static void collideFluidsWithBulletCcd(const FluidParametersGlobal &FG, FluidSph *fluid, btCollisionWorld *world);

	static void resolveCollision(const FluidParametersGlobal &FG, FluidSph *fluid, int fluidIndex, btCollisionObject *object, 
								 const btVector3 &fluidNormal, const btVector3 &hitPointWorld, btScalar distance);
								 
public:
	static void stepSimulation(FluidWorld *fluidWorld, btCollisionWorld *world, btScalar secondsElapsed)
	{
		const bool USE_VARIABLE_DELTA_TIME = false;
		if(USE_VARIABLE_DELTA_TIME)
		{
			// crashes
			FluidParametersGlobal FP = fluidWorld->getGlobalParameters();
			FP.m_timeStep = secondsElapsed;
			fluidWorld->setGlobalParameters(FP);
		}
		
		const FluidParametersGlobal &FG = fluidWorld->getGlobalParameters();
		
		for(int i = 0; i < fluidWorld->getNumFluids(); ++i)
		{
			FluidSph *fluid = fluidWorld->getFluid(i);
		
			//collideFluidsWithBullet(FG, fluid, world);
			//collideFluidsWithBullet2(FG, fluid, world);
			collideFluidsWithBulletCcd(FG, fluid, world);
		}
		
		const bool USE_ACCUMULATOR = false;
		if(USE_ACCUMULATOR)
		{
			//Prevent simulation from running too quickly
			static btScalar secondsAccumulated = 0;
			secondsAccumulated += secondsElapsed;
			
			if(secondsAccumulated > FG.m_timeStep * 2.0) secondsAccumulated = FG.m_timeStep;
				
			if(secondsAccumulated >= FG.m_timeStep)
			{
				fluidWorld->stepSimulation();
				secondsAccumulated -= FG.m_timeStep;
			}
		}
		else fluidWorld->stepSimulation();
		
		{
			static int counter = 0;
			if(++counter > 100)
			{
				counter = 0;
				
				for(int i = 0; i < fluidWorld->getNumFluids(); ++i)
					printf( "fluidWorld->getFluid(%d)->numParticles(): %d \n", i, fluidWorld->getFluid(i)->numParticles() );
			}
		}
	}
};

#endif

