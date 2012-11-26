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
#include "hfDemos.h"

#include "LinearMath/btRandom.h"
#include "LinearMath/btDefaultMotionState.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btBoxShape.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"

#include "BulletFluids/btFluidHfRigidDynamicsWorld.h"
#include "BulletFluids/btFluidHfRigidCollisionConfiguration.h"
#include "BulletFluids/Hf/btFluidHf.h"
#include "BulletFluids/Hf/btFluidHfBuoyantConvexShape.h"

#define ARRAY_SIZE_X 1
#define ARRAY_SIZE_Y 1
#define ARRAY_SIZE_Z 1

#define START_POS_X 5
#define START_POS_Y -5
#define START_POS_Z 3

void Init_Floatyness (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(0.25), 100, 100 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-10.0), btScalar(-5.0), btScalar(-10.0)) );
		
		fluid->setWorldTransform(transform);
		fluid->setHorizontalVelocityScale( btScalar(0.0f) );
		fluid->setVolumeDisplacementScale( btScalar(0.0f) );

		for(int i = 0; i < fluid->getNumNodesZ()*fluid->getNumNodesX(); i++) fluid->setFluidHeight( i, btScalar(5.0f) );
		fluid->prep();
	}
	world->addFluidHf(fluid);
	
	btConvexShape* sphereShape = new btSphereShape( btScalar(1.0) );
	const btScalar MASS = btScalar(1.0);
	
	btVector3 localInertia(0,0,0);
	sphereShape->calculateLocalInertia(MASS, localInertia);
		
	const int numObjects = 5;
	btScalar floatyness = btScalar(1.0f);
	btScalar dfloatyness = btScalar(0.25f);
	btScalar start_x = btScalar(-5.0f);
	btScalar step_x = btScalar(3.0f);
	btScalar start_z = btScalar(-5.0f);
	for (int i = 0; i < numObjects; i++)
	{
		btFluidHfBuoyantConvexShape* buoyantShape = new btFluidHfBuoyantConvexShape(sphereShape);
		buoyantShape->generateShape (btScalar(0.25f), btScalar(0.05f));
		buoyantShape->setFloatyness (floatyness + dfloatyness * i);
		collisionShapes.push_back (buoyantShape);


		btTransform transform( btQuaternion::getIdentity(), btVector3(step_x * i + start_x, 7.5f, start_z) );

		//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
		btDefaultMotionState* myMotionState = new btDefaultMotionState(transform);
		btRigidBody::btRigidBodyConstructionInfo rbInfo(MASS,myMotionState,buoyantShape,localInertia);
		btRigidBody* body = new btRigidBody(rbInfo);
		world->addRigidBody(body);
	}
	
	floatyness = btScalar(2.0f);
	start_z = btScalar(5.0f);
	for (int i = 0; i < numObjects; i++)
	{
		btFluidHfBuoyantConvexShape* buoyantShape = new btFluidHfBuoyantConvexShape(sphereShape);
		buoyantShape->generateShape (btScalar(0.25f), btScalar(0.05f));
		buoyantShape->setFloatyness (floatyness + dfloatyness * i);
		collisionShapes.push_back (buoyantShape);
		
		btTransform transform( btQuaternion::getIdentity(), btVector3(step_x * i + start_x, -4.0f, start_z) );
		
		//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
		btDefaultMotionState* myMotionState = new btDefaultMotionState(transform);
		btRigidBody::btRigidBodyConstructionInfo rbInfo(MASS,myMotionState,buoyantShape,localInertia);
		btRigidBody* body = new btRigidBody(rbInfo);
		world->addRigidBody(body);
	}
}

void Init_Bowl (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf (btScalar(1.0), 50, 50);
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-10.0), btScalar(-5.0), btScalar(-10.0)) );
		fluid->setWorldTransform(transform);
		
		btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray();
		btAlignedObjectArray<btScalar>& eta = fluid->getFluidArray();
		btScalar amplitude = btScalar(200.0);
		for (int i = 0; i < fluid->getNumNodesX(); i++)
		{
			btScalar x = btScalar(i - fluid->getNumNodesX()/2)/btScalar(fluid->getNumNodesX()*2);
			btScalar xh = amplitude * (x * x) + btScalar(5.0);
			for (int j = 0; j < fluid->getNumNodesZ(); j++)
			{
				btScalar y = btScalar(j - fluid->getNumNodesZ()/2)/btScalar(fluid->getNumNodesZ()*2);
				btScalar yh = amplitude * (y * y) + btScalar(5.0);
				btScalar groundHeight = btMax(xh,yh);
				
				int index = fluid->arrayIndex(i, j);
				ground[index] = groundHeight;
				
				btScalar waterHeight = btScalar(0.0f);
				if(groundHeight > 14.0) waterHeight = btScalar(0.0f); 
				else waterHeight = btScalar(14.0f) - groundHeight;
				
				eta[index] = waterHeight;		
			}
		}

		fluid->prep ();
	}
	world->addFluidHf (fluid);
	
	{
		//create a few dynamic rigidbodies
		// Re-using the same collision is better for memory usage and performance

		btConvexShape* colShape = new btBoxShape(btVector3(1,1,1));
		btFluidHfBuoyantConvexShape* buoyantShape = new btFluidHfBuoyantConvexShape(colShape);
		buoyantShape->generateShape (btScalar(0.25f), btScalar(0.05f));
		//btCollisionShape* colShape = new btSphereShape(btScalar(1.));
		collisionShapes.push_back(colShape);
		collisionShapes.push_back (buoyantShape);

		/// Create Dynamic Objects
		btTransform startTransform;
		startTransform.setIdentity();

		btScalar	mass(1.f);

		//rigidbody is dynamic if and only if mass is non zero, otherwise static
		bool isDynamic = (mass != 0.f);

		btVector3 localInertia(0,0,0);
		if (isDynamic)
			colShape->calculateLocalInertia(mass,localInertia);

		float start_x = START_POS_X - ARRAY_SIZE_X/2;
		float start_y = START_POS_Y;
		float start_z = START_POS_Z - ARRAY_SIZE_Z/2;

		for (int k=0;k<ARRAY_SIZE_Y;k++)
		{
			for (int i=0;i<ARRAY_SIZE_X;i++)
			{
				for(int j = 0;j<ARRAY_SIZE_Z;j++)
				{
					btVector3 origin( btScalar(2.0)*i + start_x, 
										btScalar(10.0)+btScalar(2.0)*k + start_y, 
										btScalar(2.0)*j + start_z );

					startTransform.setOrigin(origin);
					
					//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
					btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
					btRigidBody::btRigidBodyConstructionInfo rbInfo(mass,myMotionState,buoyantShape,localInertia);
					btRigidBody* body = new btRigidBody(rbInfo);
					world->addRigidBody(body);
				}
			}
		}
	}
}

void Init_Drops (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(0.5), 50, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-10.0), btScalar(-5.0), btScalar(-10.0)) );
		fluid->setWorldTransform(transform);

		for (int i = 0; i < fluid->getNumNodesZ()*fluid->getNumNodesX(); i++) fluid->setFluidHeight(i, btScalar(5.0f));
		fluid->prep();
	}
	world->addFluidHf (fluid);
	
	{
		//create a few dynamic rigidbodies
		// Re-using the same collision is better for memory usage and performance

		btConvexShape* colShape = new btBoxShape( btVector3(5, 0.5, 5) );
		btFluidHfBuoyantConvexShape* buoyantShape = new btFluidHfBuoyantConvexShape(colShape);
		buoyantShape->generateShape( btScalar(0.25f), btScalar(0.05f) );
		
		collisionShapes.push_back(colShape);
		collisionShapes.push_back(buoyantShape);

		/// Create Dynamic Objects
		btTransform startTransform;
		startTransform.setIdentity();

		const btScalar MASS(1.f);

		//rigidbody is dynamic if and only if mass is non zero, otherwise static
		bool isDynamic = (MASS != 0.f);

		btVector3 localInertia(0,0,0);
		if (isDynamic) colShape->calculateLocalInertia(MASS,localInertia);

		float start_x = START_POS_X - ARRAY_SIZE_X/2;
		float start_y = START_POS_Y;
		float start_z = START_POS_Z - ARRAY_SIZE_Z/2;

		for (int k=0;k<ARRAY_SIZE_Y;k++)
		{
			for (int i=0;i<ARRAY_SIZE_X;i++)
			{
				for(int j = 0;j<ARRAY_SIZE_Z;j++)
				{
					btVector3 origin( btScalar(2.0)*i + start_x, 
										btScalar(10.0)+btScalar(2.0)*k + start_y, 
										btScalar(2.0)*j + start_z );

					startTransform.setOrigin(origin);

			
					//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
					btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
					btRigidBody::btRigidBodyConstructionInfo rbInfo(MASS,myMotionState,buoyantShape,localInertia);
					btRigidBody* body = new btRigidBody(rbInfo);
					world->addRigidBody(body);
				}
			}
		}
	}
}

void Init_Wave (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0f), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);

		btAlignedObjectArray<btScalar>& eta = fluid->getFluidArray();
		for (int i = 0; i < fluid->getNumNodesZ()*fluid->getNumNodesX(); i++) eta[i] = btScalar(10.0f);

		for (int i = 1; i < fluid->getNumNodesX()-1; i++)
		{
			eta[ fluid->arrayIndex(i, fluid->getNumNodesZ()/2-1) ] = btScalar(2.0);
			eta[ fluid->arrayIndex(i, fluid->getNumNodesZ()/2) ] = btScalar(2.0);
			eta[ fluid->arrayIndex(i, fluid->getNumNodesZ()/2+1) ] = btScalar(2.0);
		}

		fluid->prep();
	}
	world->addFluidHf(fluid);
}

void Init_RandomDrops (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		for (int i = 0; i < fluid->getNumNodesZ()*fluid->getNumNodesX(); i++) fluid->getFluidArray()[i] = btScalar(0.0f);
		
		fluid->prep();
	}
	world->addFluidHf (fluid);
}
void Run_RandomDrops (btFluidHfRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0f);
	btScalar dt = btScalar(1.0/60.);	
	
	if (dtSinceLastDrop > btScalar(0.5f))
	{
		btAlignedObjectArray<btFluidHf*>& fluids = world->getFluidHfArray();
		btFluidHf* fluid = fluids[0];
		
		dtSinceLastDrop = btScalar(0.0f);
		int randomXNode = GEN_rand () % (fluid->getNumNodesX()-2);
		int randomZNode = GEN_rand () % (fluid->getNumNodesZ()-2);
		if (randomXNode <= 1)
			randomXNode = 2;
		if (randomZNode <= 1)
			randomZNode = 2;

		btAlignedObjectArray<btScalar>& eta = fluid->getFluidArray ();
		btAlignedObjectArray<btScalar>& height = fluid->getCombinedHeightArray ();
		const btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray ();
		btAlignedObjectArray<bool>& flags = fluid->getFlagsArray();
		int index = fluid->arrayIndex (randomXNode, randomZNode);
		eta[index] += btScalar(4.5f);
		eta[index-1] += btScalar(2.25f);
		eta[index+1] += btScalar(2.25f);
		eta[index+fluid->getNumNodesX()] += btScalar(2.25f);
		eta[index-fluid->getNumNodesX()] += btScalar(2.25f);
		height[index] = eta[index] + ground[index];
		height[index-1] = eta[index-1] + ground[index-1];
		height[index+1] = eta[index+1] + ground[index+1];
		height[index+fluid->getNumNodesX()] = eta[index+fluid->getNumNodesX()] + ground[index+fluid->getNumNodesX()];
		height[index-fluid->getNumNodesX()] = eta[index-fluid->getNumNodesX()] + ground[index-fluid->getNumNodesX()];
		flags[index] = true;
		flags[index-1] = true;
		flags[index+1] = true;
		flags[index+fluid->getNumNodesX()] = true;
		flags[index-fluid->getNumNodesX()] = true;
	} 
	else 
	{
		dtSinceLastDrop += dt;
	}
}

void Init_FillPool (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	const int gridLength = 50;
	const int gridWidth = 50;
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), gridLength, gridWidth );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-20.0), btScalar(-5.0), btScalar(-20.0)) );
		fluid->setWorldTransform(transform);
		
		btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray();

		const btScalar poolEdgeHeight = btScalar(10.0f);
		const btScalar poolBottomHeight = btScalar(1.0f);
		const btScalar poolPourerHeight = btScalar(6.0f);
		for (int j = 1; j < fluid->getNumNodesZ()-1; j++)
		{
			for (int i = 1; i < fluid->getNumNodesX()-1; i++)
			{
				int index = fluid->arrayIndex (i, j);
				
				// pool edge
				//if (j == 1 || i == 1 || j == fluid->getNumNodesZ()-2 || i == fluid->getNumNodesX()-2)
				if (j <= 1+3 || i <= 1+8 || j >= fluid->getNumNodesZ()-2 || i >= fluid->getNumNodesX()-2-3)
				{
					ground[index] = poolEdgeHeight;
					continue;
				}
				if (j > 35)
				{
					if (i <= 25 || i >= 30) ground[index] = poolEdgeHeight;
					else ground[index] = poolPourerHeight;
					
					continue;
				}
				ground[index] = poolBottomHeight;
			}
		}
		
		fluid->prep();
	}
	world->addFluidHf(fluid);
	
	const int SPAWN_BOXES = 0;
	if(SPAWN_BOXES)
	{
		btConvexShape* colShape = new btBoxShape(btVector3(btScalar(1.), btScalar(1.), btScalar(1.)));
		btFluidHfBuoyantConvexShape* buoyantShape = new btFluidHfBuoyantConvexShape(colShape);
		buoyantShape->generateShape (btScalar(0.25f), btScalar(0.05f));
		collisionShapes.push_back(colShape);
		collisionShapes.push_back(buoyantShape);

		/// Create Dynamic Objects
		btTransform startTransform;
		startTransform.setIdentity();

		btScalar	mass(1.f);

		//rigidbody is dynamic if and only if mass is non zero, otherwise static
		bool isDynamic = (mass != 0.f);

		btVector3 localInertia(0,0,0);
		if (isDynamic)
			colShape->calculateLocalInertia(mass,localInertia);

		int gridSize = 2;
		btScalar startPosX = btScalar(-10.0f);
		btScalar startPosY = btScalar(2.0f);
		btScalar startPosZ = btScalar(-10.f);
		float start_x = startPosX - gridSize/2;
		float start_y = startPosY;
		float start_z = startPosZ - gridSize/2;
		
		for (int k=0;k<gridSize;k++)
		{
			for (int i=0;i<gridSize;i++)
			{
				for(int j = 0;j<gridSize;j++)
				{
					btVector3 origin( btScalar(2.0)*i + start_x, 
										btScalar(10.0)+btScalar(2.0)*k + start_y, 
										btScalar(2.0)*j + start_z );

					startTransform.setOrigin(origin);

			
					//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
					btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
					btRigidBody::btRigidBodyConstructionInfo rbInfo(mass,myMotionState,buoyantShape,localInertia);
					btRigidBody* body = new btRigidBody(rbInfo);
					world->addRigidBody(body);
				}
			}
		}
	}
}
void Run_FillPool (btFluidHfRigidDynamicsWorld* world)
{
	btAlignedObjectArray<btFluidHf*>& fluids = world->getFluidHfArray();
	btFluidHf* fluid = fluids[0];

	for (int i = 26; i < 30; i++) fluid->setFluidHeight( fluid->arrayIndex(i, fluid->getNumNodesZ()-3), btScalar(3.0f) );
	
	const int DEBUG = 0;
	if(DEBUG)
	{
		for (int j = 1; j < fluid->getNumNodesZ()-1; j++)
			for (int i = 1; i < fluid->getNumNodesX()-1; i++)
			{
				if(j <= 1+15) 
				{
					int index = fluid->arrayIndex(i, j);
					fluid->getVelocityXArray()[index] = btScalar(0.0);
					fluid->getVelocityZArray()[index] = btScalar(0.0);
				}
			}
	}
}


void Init_Fill (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0f), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for (int i = 0; i < fluid->getNumNodesZ()*fluid->getNumNodesX(); i++) fluid->getFluidArray()[i] = btScalar(0.0f);
		
		fluid->prep ();
	}
	world->addFluidHf (fluid);
}

void Run_Fill (btFluidHfRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0f);
	btScalar dt = btScalar(1.0/60.);

	btAlignedObjectArray<btFluidHf*>& fluids = world->getFluidHfArray ();
	btFluidHf* fluid = fluids[0];

	if (dtSinceLastDrop > btScalar(0.25f))
	{
		dtSinceLastDrop = btScalar(0.0f);

		btAlignedObjectArray<btScalar>& eta = fluid->getFluidArray ();
		btAlignedObjectArray<btScalar>& velocityU = fluid->getVelocityXArray ();
		btAlignedObjectArray<btScalar>& velocityV = fluid->getVelocityZArray ();
		btAlignedObjectArray<btScalar>& height = fluid->getCombinedHeightArray ();
		const btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray ();
		btAlignedObjectArray<bool>& flags = fluid->getFlagsArray();
		int index = fluid->arrayIndex (fluid->getNumNodesX()/2, fluid->getNumNodesZ()/2);
		eta[index] += btScalar(4.5f);
		eta[index-1] += btScalar(2.25f);
		eta[index+1] += btScalar(2.25f);
		eta[index+fluid->getNumNodesX()] += btScalar(2.25f);
		eta[index-fluid->getNumNodesX()] += btScalar(2.25f);

		velocityU[index] = btScalar(0.0f);
		velocityU[index-1] = btScalar(-10.0f);
		velocityU[index+1] = btScalar(10.0f);
		velocityU[index+fluid->getNumNodesX()] = btScalar(0.0f);
		velocityU[index-fluid->getNumNodesX()] = btScalar(0.0f);

		velocityV[index] = btScalar(0.0f);
		velocityV[index-1] = btScalar(0.0f);
		velocityV[index+1] = btScalar(0.0f);
		velocityV[index+fluid->getNumNodesX()] = btScalar(10.0f);
		velocityV[index-fluid->getNumNodesX()] = btScalar(-10.0f);

		height[index] = eta[index] + ground[index];
		height[index-1] = eta[index-1] + ground[index-1];
		height[index+1] = eta[index+1] + ground[index+1];
		height[index+fluid->getNumNodesX()] = eta[index+fluid->getNumNodesX()] + ground[index+fluid->getNumNodesX()];
		height[index-fluid->getNumNodesX()] = eta[index-fluid->getNumNodesX()] + ground[index-fluid->getNumNodesX()];
		flags[index] = true;
		flags[index-1] = true;
		flags[index+1] = true;
		flags[index+fluid->getNumNodesX()] = true;
		flags[index-fluid->getNumNodesX()] = true;
	} 
	else 
	{
		dtSinceLastDrop += dt;
	}
	
}

void Init_BlockWave (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf (btScalar(1.0), 75, 50);
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);

		btAlignedObjectArray<btScalar>& eta = fluid->getFluidArray();

		for (int i = 0; i < fluid->getNumNodesZ() * fluid->getNumNodesX(); i++) eta[i] = btScalar(12.0f);

		for (int i = fluid->getNumNodesX()/8; i < fluid->getNumNodesX()/4; i++)
		{
			for (int j = fluid->getNumNodesZ()/8; j < fluid->getNumNodesZ()/4; j++)
			{
				int index = fluid->arrayIndex(i, j);
				eta[index] = btScalar(4.0f);
			}
		}
		fluid->prep ();
	}
	world->addFluidHf (fluid);
	
	{
		btConvexShape* colShape = new btSphereShape(btScalar(1.));
		btFluidHfBuoyantConvexShape* buoyantShape = new btFluidHfBuoyantConvexShape(colShape);
		buoyantShape->generateShape (btScalar(0.25f), btScalar(0.05f));
		collisionShapes.push_back(buoyantShape);
		collisionShapes.push_back(colShape);

		/// Create Dynamic Objects
		btTransform startTransform;
		startTransform.setIdentity();

		btScalar	mass(1.f);

		//rigidbody is dynamic if and only if mass is non zero, otherwise static
		bool isDynamic = (mass != 0.f);

		btVector3 localInertia(0,0,0);
		if (isDynamic)
			colShape->calculateLocalInertia(mass,localInertia);

		int gridSize = 2;
		btScalar startPosX = btScalar(-10.0f);
		btScalar startPosY = btScalar(2.0f);
		btScalar startPosZ = btScalar(-10.f);
		float start_x = startPosX - gridSize/2;
		float start_y = startPosY;
		float start_z = startPosZ - gridSize/2;
		
		for (int k=0;k<gridSize;k++)
		{
			for (int i=0;i<gridSize;i++)
			{
				for(int j = 0;j<gridSize;j++)
				{
					btVector3 origin( btScalar(2.0)*i + start_x, 
										btScalar(10.0)+btScalar(2.0)*k + start_y, 
										btScalar(2.0)*j + start_z );

					startTransform.setOrigin(origin);

			
					//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
					btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
					btRigidBody::btRigidBodyConstructionInfo rbInfo(mass,myMotionState,buoyantShape,localInertia);
					btRigidBody* body = new btRigidBody(rbInfo);
					world->addRigidBody(body);
				}
			}
		}
	}
}

void Init_Ground (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0f), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);

		btAlignedObjectArray<btScalar>& eta = fluid->getFluidArray();

		for (int i = 0; i < fluid->getNumNodesZ() * fluid->getNumNodesX(); i++) eta[i] = btScalar(4.0f);

		btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray();
		for (int i = 0; i < fluid->getNumNodesX(); i++)
		{
			for (int j = 0; j < fluid->getNumNodesZ(); j++)
			{
				int index = fluid->arrayIndex(i, j);

				if (j <= fluid->getNumNodesZ()/2) ground[index] = btScalar(5.0f);
				else if (j > (fluid->getNumNodesZ()/8*6)) ground[index] = btScalar(0.0f);
				else ground[index] = btScalar(6.5f);
				
				if (j <= fluid->getNumNodesZ()/4 && j > fluid->getNumNodesZ()/8) eta[index] = btScalar(8.0f);
				else if (j <= fluid->getNumNodesZ()/8) eta[index] = btScalar(20.0f);
				else eta[index] = btScalar(0.0f);
			}
		}
		fluid->prep();
	}
	
	world->addFluidHf(fluid);
}

void Init_Ground2 (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf (btScalar(1.0f), 75, 50);
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(5.0), btScalar(-50.0)) );
		fluid->setWorldTransform (transform);
		
		btAlignedObjectArray<btScalar>& eta = fluid->getFluidArray();
		btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray();
		
		for (int i = 0; i < fluid->getNumNodesX(); i++)
		{
			for (int j = 0; j < fluid->getNumNodesZ(); j++)
			{
				int index = fluid->arrayIndex(i, j);

				ground[index] = (btScalar(j)/fluid->getNumNodesZ()-1)*btScalar(8.0f);
			}
		}

		for (int i = 0; i < fluid->getNumNodesZ() * fluid->getNumNodesX(); i++) eta[i] = btScalar(2.0f);

		for (int i = fluid->getNumNodesX()/8; i < fluid->getNumNodesX()/4; i++)
		{
			for (int j = fluid->getNumNodesZ()/8; j < fluid->getNumNodesZ()/4; j++)
			{
				int index = fluid->arrayIndex(i, j);
				eta[index] = btScalar(8.0f);
			}
		}
		fluid->prep();
	}
	world->addFluidHf(fluid);
}

void Init_Fill2 (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), 100, 100 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for (int i = 0; i < fluid->getNumNodesZ()*fluid->getNumNodesX(); i++) fluid->getFluidArray()[i] = btScalar(0.0f);
		
		fluid->prep();
	}
	world->addFluidHf (fluid);
}
void Run_Fill2 (btFluidHfRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0f);
	btScalar dt = btScalar(1.0/60.);

	btAlignedObjectArray<btFluidHf*>& fluids = world->getFluidHfArray ();
	btFluidHf* fluid = fluids[0];

	if (dtSinceLastDrop > btScalar(0.25f))
	{
		dtSinceLastDrop = btScalar(0.0f);

		btAlignedObjectArray<btScalar>& eta = fluid->getFluidArray ();
		btAlignedObjectArray<btScalar>& velocityU = fluid->getVelocityXArray ();
		btAlignedObjectArray<btScalar>& height = fluid->getCombinedHeightArray ();
		const btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray ();
		btAlignedObjectArray<bool>& flags = fluid->getFlagsArray();

		for (int i = 1; i < fluid->getNumNodesX()-1; i++)
		{
			int index = fluid->arrayIndex (i, 1);
			eta[index] += btScalar(3.0f);
			velocityU[index] = btScalar(4.0f);
			height[index] = ground[index] + eta[index];
			flags[index] = true;
		}
	} 
	else 
	{
		dtSinceLastDrop += dt;
	}
	
}

void Init_MovingPour (btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for (int i = 0; i < fluid->getNumNodesZ()*fluid->getNumNodesX(); i++) fluid->getFluidArray()[i] = btScalar(5.0f);
		
		fluid->prep();
	}
	world->addFluidHf(fluid);
	
	{
		//create a few dynamic rigidbodies
		// Re-using the same collision is better for memory usage and performance

		btCollisionShape* colShape = new btBoxShape(btVector3(1,1,1));
		//btCollisionShape* colShape = new btSphereShape(btScalar(1.));
		collisionShapes.push_back(colShape);

		/// Create Dynamic Objects
		btTransform startTransform;
		startTransform.setIdentity();

		btScalar	mass(1.f);

		//rigidbody is dynamic if and only if mass is non zero, otherwise static
		bool isDynamic = (mass != 0.f);

		btVector3 localInertia(0,0,0);
		if (isDynamic)
			colShape->calculateLocalInertia(mass,localInertia);

		float start_x = START_POS_X - ARRAY_SIZE_X/2;
		float start_y = START_POS_Y;
		float start_z = START_POS_Z - ARRAY_SIZE_Z/2;

		for (int k=0;k<ARRAY_SIZE_Y;k++)
		{
			for (int i=0;i<ARRAY_SIZE_X;i++)
			{
				for(int j = 0;j<ARRAY_SIZE_Z;j++)
				{
					btVector3 origin( btScalar(2.0)*i + start_x, 
										btScalar(10.0)+btScalar(2.0)*k + start_y, 
										btScalar(2.0)*j + start_z );

					startTransform.setOrigin(origin);

			
					//using motionstate is recommended, it provides interpolation capabilities, and only synchronizes 'active' objects
					btDefaultMotionState* myMotionState = new btDefaultMotionState(startTransform);
					btRigidBody::btRigidBodyConstructionInfo rbInfo(mass,myMotionState,colShape,localInertia);
					btRigidBody* body = new btRigidBody(rbInfo);
					world->addRigidBody(body);
				}
			}
		}
	}
}

void Run_MovingPour(btFluidHfRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0f);
	static btScalar x = 4;
	static btScalar z = 4;
	static btScalar dx = btScalar(20.0f);
	static btScalar dz = btScalar(30.0f);
	btScalar dt = btScalar(1.0/60.);

	btAlignedObjectArray<btFluidHf*>& fluids = world->getFluidHfArray ();
	btFluidHf* fluid = fluids[0];

	int minX = 2;
	int minZ = 2;
	int maxX = fluid->getNumNodesX() - 2;
	int maxZ = fluid->getNumNodesZ() - 2;

	x += dx * dt;
	
	if (x <= minX)
	{
		dx *= btScalar(-1.0f);
		x = static_cast<btScalar>(minX);
	} 
	else if (x >= maxX) 
	{
		dx *= btScalar(-1.0f);
		x = static_cast<btScalar>(maxX);
	}
	z += dz * dt;
	
	if (z <= minZ)
	{
		dz *= btScalar(-1.0f);
		z =static_cast<btScalar>(minZ);
	} 
	else if (z >= maxZ) 
	{
		dz *= btScalar(-1.0f);
		z = static_cast<btScalar>(maxZ);
	}

	const btScalar dropHeight = btScalar(3.0f);
	
	{
		int iX = (int)x;
		int iZ = (int)z;
		fluid->addFluidHeight (iX,iZ, dropHeight);
		//fluid->addFluidHeight (x, z+1, dropHeight);
		//fluid->addFluidHeight (x+1, z, dropHeight);
		//fluid->addFluidHeight (x+1, z+1, dropHeight);
	}
}