/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2009 Erwin Coumans  http://bulletphysics.com

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.

Experimental Buoyancy fluid demo written by John McCutchan
*/
//This is an altered source version based on the HeightFieldFluidDemo included with Bullet Physics 2.80(bullet-2.80-rev2531).

#include "btFluidHfRigidDynamicsWorld.h"

#include <stdio.h>
#include "LinearMath/btQuickprof.h"
#include "LinearMath/btIDebugDraw.h"
#include "BulletCollision/CollisionShapes/btCollisionShape.h"

#include "Hf/btFluidHf.h"
#include "Hf/btFluidHfBuoyantConvexShape.h"


btFluidHfRigidDynamicsWorld::btFluidHfRigidDynamicsWorld(btDispatcher* dispatcher, btBroadphaseInterface* pairCache,
														btConstraintSolver* constraintSolver, btCollisionConfiguration* collisionConfiguration,
														btFluidSphSolver* sphSolver)
: btFluidRigidDynamicsWorld(dispatcher, pairCache, constraintSolver, collisionConfiguration, sphSolver)
{
	m_drawMode = DRAWMODE_NORMAL;
	m_bodyDrawMode = BODY_DRAWMODE_NORMAL;
}
		
void btFluidHfRigidDynamicsWorld::internalSingleStepSimulation(btScalar timeStep)
{
	btFluidRigidDynamicsWorld::internalSingleStepSimulation(timeStep);

	{
		BT_PROFILE("updateHfFluids");
	
		for(int i = 0; i < m_hfFluids.size(); i++)
		{
			btFluidHf* hfFluid = m_hfFluids[i];
			hfFluid->stepSimulation(timeStep);
		}
	}
}


void btFluidHfRigidDynamicsWorld::addFluidHf(btFluidHf* body)
{
	m_hfFluids.push_back(body);

	btCollisionWorld::addCollisionObject(body, btBroadphaseProxy::DefaultFilter, btBroadphaseProxy::AllFilter);
}

void btFluidHfRigidDynamicsWorld::removeFluidHf(btFluidHf* body)
{
	m_hfFluids.remove(body);

	btCollisionWorld::removeCollisionObject(body);
}

void btFluidHfRigidDynamicsWorld::drawFluidHfGround (btIDebugDraw* debugDraw, btFluidHf* fluid)
{
	btVector3 com = fluid->getWorldTransform().getOrigin();
	btVector3 color = btVector3(btScalar(0.13f), btScalar(0.13f), btScalar(0.0));
	for (int i = 1; i < fluid->getNumNodesX()-1; i++)
	{
		for (int j = 1; j < fluid->getNumNodesZ()-1; j++)
		{
			int sw = fluid->arrayIndex (i, j);
			int se = fluid->arrayIndex (i+1, j);
			int nw = fluid->arrayIndex (i, j+1);
			int ne = fluid->arrayIndex (i+1, j+1);
			btVector3 swV = btVector3 (fluid->getCellPosX (i), fluid->getGroundHeight(sw), fluid->getCellPosZ (j));
			btVector3 seV = btVector3 (fluid->getCellPosX (i+1), fluid->getGroundHeight(se), fluid->getCellPosZ (j));
			btVector3 nwV = btVector3 (fluid->getCellPosX (i), fluid->getGroundHeight(nw), fluid->getCellPosZ (j+1));
			btVector3 neV = btVector3 (fluid->getCellPosX (i+1), fluid->getGroundHeight(ne), fluid->getCellPosZ (j+1));
			debugDraw->drawTriangle (swV+com, seV+com, nwV+com, color, btScalar(1.0f));
			debugDraw->drawTriangle (seV+com, neV+com, nwV+com, color, btScalar(1.0f));
		}
	}
}

void btFluidHfRigidDynamicsWorld::drawFluidHfVelocity (btIDebugDraw* debugDraw, btFluidHf* fluid)
{
	const btVector3& origin = fluid->getWorldTransform().getOrigin();
	const btVector3 red = btVector3(btScalar(1.0f), btScalar(0.0f), btScalar(0.0));
	const btVector3 green = btVector3(btScalar(0.0f), btScalar(1.0f), btScalar(0.0));
	for (int i = 1; i < fluid->getNumNodesX()-1; i++)
	{
		for (int j = 1; j < fluid->getNumNodesZ()-1; j++)
		{
			int index = fluid->arrayIndex (i, j);
			if( !fluid->isActive(index) ) continue;
			btVector3 from = btVector3 ( fluid->getCellPosX(i), fluid->getCombinedHeight(index) + btScalar(0.1f), fluid->getCellPosZ(j) );

			btVector3 velocity( fluid->getVelocityX(index), btScalar(0.0), fluid->getVelocityZ(index) );
			velocity.normalize();
			btVector3 to = from + velocity;
			
			debugDraw->drawLine (from+origin, to+origin, red, green);
		}
	}
}

void btFluidHfRigidDynamicsWorld::drawFluidHfBuoyantConvexShape (btIDebugDraw* debugDrawer, btCollisionObject* object, btFluidHfBuoyantConvexShape* buoyantShape, int voxelDraw)
{
	if (voxelDraw)
	{
		const btTransform& xform = object->getWorldTransform();
		for (int i = 0; i < buoyantShape->getNumVoxels(); i++)
		{
			btVector3 p = xform * buoyantShape->getVoxelPosition(i);
			debugDrawer->drawSphere( p, buoyantShape->getVoxelRadius(), btVector3(1.0, 0.0, 0.0) );
		}
	} 
	else 
	{
		btVector3 color(btScalar(255.),btScalar(255.),btScalar(255.));
		switch(object->getActivationState())
		{
			case  ACTIVE_TAG:
				color = btVector3(btScalar(255.),btScalar(255.),btScalar(255.)); break;
			case ISLAND_SLEEPING:
				color =  btVector3(btScalar(0.),btScalar(255.),btScalar(0.));break;
			case WANTS_DEACTIVATION:
				color = btVector3(btScalar(0.),btScalar(255.),btScalar(255.));break;
			case DISABLE_DEACTIVATION:
				color = btVector3(btScalar(255.),btScalar(0.),btScalar(0.));break;
			case DISABLE_SIMULATION:
				color = btVector3(btScalar(255.),btScalar(255.),btScalar(0.));break;
			default:
			{
				color = btVector3(btScalar(255.),btScalar(0.),btScalar(0.));
			}
		};

		btConvexShape* convexShape = ((btFluidHfBuoyantConvexShape*)object->getCollisionShape())->getConvexShape();
		debugDrawObject(object->getWorldTransform(),(btCollisionShape*)convexShape,color);
	}
}

void btFluidHfRigidDynamicsWorld::drawFluidHfNormal (btIDebugDraw* debugDraw, btFluidHf* fluid)
{
	const btVector3& com = fluid->getWorldTransform().getOrigin();
	for (int i = 0; i < fluid->getNumNodesX()-1; i++)
	{
		for (int j = 0; j < fluid->getNumNodesZ()-1; j++)
		{
			int sw = fluid->arrayIndex (i, j);

			btScalar h = fluid->getFluidHeight(sw);
			btScalar g = fluid->getGroundHeight(sw);

			if( h < btScalar(0.03f) ) continue;

			btVector3 boxMin = btVector3(fluid->getCellPosX (i), g, fluid->getCellPosZ(j));
			btVector3 boxMax = btVector3(fluid->getCellPosX(i+1), g+h, fluid->getCellPosZ(j+1));
			boxMin += com;
			boxMax += com;
			
			debugDraw->drawBox (boxMin, boxMax, btVector3(btScalar(0.0f), btScalar(0.0f), btScalar(1.0f)));
		}
	}
}

void btFluidHfRigidDynamicsWorld::debugDrawWorld()
{
	if (getDebugDrawer())
	{
		int i;
		for (  i=0;i<this->m_hfFluids.size();i++)
		{
			btFluidHf*	phh=(btFluidHf*)this->m_hfFluids[i];
			switch (m_drawMode)
			{
			case DRAWMODE_NORMAL:
				//drawFluidHfGround (m_debugDrawer, phh);
				//drawFluidHfNormal (m_debugDrawer, phh);
			break;
			case DRAWMODE_VELOCITY:
				//drawFluidHfGround (m_debugDrawer, phh);
				//drawFluidHfNormal (m_debugDrawer, phh);
				drawFluidHfVelocity (m_debugDrawer, phh);
			break;
			default:
				btAssert (0);
			break;
			}
		}	
		for (i = 0; i < this->m_collisionObjects.size(); i++)
		{
			btCollisionShape* shape = m_collisionObjects[i]->getCollisionShape();
			if (shape->getShapeType() == HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE)
			{
				btFluidHfBuoyantConvexShape* buoyantShape = (btFluidHfBuoyantConvexShape*)shape;
				drawFluidHfBuoyantConvexShape (m_debugDrawer, m_collisionObjects[i], buoyantShape, m_bodyDrawMode);
			}
		}
	}
	
	btFluidRigidDynamicsWorld::debugDrawWorld();
}
