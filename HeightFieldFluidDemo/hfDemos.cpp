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

#include "BulletHfFluid/btHfFluidRigidDynamicsWorld.h"
#include "BulletHfFluid/btHfFluid.h"
#include "BulletHfFluid/btHfFluidRigidCollisionConfiguration.h"
#include "BulletHfFluid/btHfFluidBuoyantConvexShape.h"

#define ARRAY_SIZE_X 1
#define ARRAY_SIZE_Y 1
#define ARRAY_SIZE_Z 1

#define START_POS_X 5
#define START_POS_Y -5
#define START_POS_Z 3

void Init_Floatyness (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btHfFluid* fluid = new btHfFluid( btScalar(0.25), 100, 100 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-10.0), btScalar(-5.0), btScalar(-10.0)) );
		
		fluid->setWorldTransform(transform);
		fluid->setHorizontalVelocityScale( btScalar(0.0f) );
		fluid->setVolumeDisplacementScale( btScalar(0.0f) );

		for(int i = 0; i < fluid->getNumNodesLength()*fluid->getNumNodesWidth(); i++) fluid->setFluidHeight( i, btScalar(5.0f) );
		fluid->prep();
	}
	world->addHfFluid(fluid);
	
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
		btHfFluidBuoyantConvexShape* buoyantShape = new btHfFluidBuoyantConvexShape(sphereShape);
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
		btHfFluidBuoyantConvexShape* buoyantShape = new btHfFluidBuoyantConvexShape(sphereShape);
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

void Init_Bowl (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btHfFluid* fluid = new btHfFluid (btScalar(1.0), 50, 50);
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-10.0), btScalar(-5.0), btScalar(-10.0)) );
		fluid->setWorldTransform(transform);
		
		btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray();
		btAlignedObjectArray<btScalar>& eta = fluid->getEtaArray();
		btScalar amplitude = btScalar(200.0);
		for (int i = 0; i < fluid->getNumNodesWidth(); i++)
		{
			btScalar x = btScalar(i - fluid->getNumNodesWidth()/2)/btScalar(fluid->getNumNodesWidth()*2);
			btScalar xh = amplitude * (x * x) + btScalar(5.0);
			for (int j = 0; j < fluid->getNumNodesLength(); j++)
			{
				btScalar y = btScalar(j - fluid->getNumNodesLength()/2)/btScalar(fluid->getNumNodesLength()*2);
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
	world->addHfFluid (fluid);
	
	{
		//create a few dynamic rigidbodies
		// Re-using the same collision is better for memory usage and performance

		btConvexShape* colShape = new btBoxShape(btVector3(1,1,1));
		btHfFluidBuoyantConvexShape* buoyantShape = new btHfFluidBuoyantConvexShape(colShape);
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
					startTransform.setOrigin( btVector3(2.0*i + start_x, 10+2.0*k + start_y, 2.0*j + start_z) );
					
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

void Init_Drops (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btHfFluid* fluid = new btHfFluid( btScalar(0.5), 50, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-10.0), btScalar(-5.0), btScalar(-10.0)) );
		fluid->setWorldTransform(transform);

		for (int i = 0; i < fluid->getNumNodesLength()*fluid->getNumNodesWidth(); i++) fluid->setFluidHeight(i, btScalar(5.0f));
		fluid->prep();
	}
	world->addHfFluid (fluid);
	
	{
		//create a few dynamic rigidbodies
		// Re-using the same collision is better for memory usage and performance

		btConvexShape* colShape = new btBoxShape( btVector3(5, 0.5, 5) );
		btHfFluidBuoyantConvexShape* buoyantShape = new btHfFluidBuoyantConvexShape(colShape);
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
					startTransform.setOrigin( btVector3(2.0*i + start_x, 10+2.0*k + start_y, 2.0*j + start_z) );

			
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

void Init_Wave (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btHfFluid* fluid = new btHfFluid( btScalar(1.0f), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);

		btAlignedObjectArray<btScalar>& eta = fluid->getEtaArray();
		for (int i = 0; i < fluid->getNumNodesLength()*fluid->getNumNodesWidth(); i++) eta[i] = btScalar(10.0f);

		for (int i = 1; i < fluid->getNumNodesWidth()-1; i++)
		{
			eta[ fluid->arrayIndex(i, fluid->getNumNodesLength()/2-1) ] = btScalar(2.0);
			eta[ fluid->arrayIndex(i, fluid->getNumNodesLength()/2) ] = btScalar(2.0);
			eta[ fluid->arrayIndex(i, fluid->getNumNodesLength()/2+1) ] = btScalar(2.0);
		}

		fluid->prep();
	}
	world->addHfFluid(fluid);
}

void Init_RandomDrops (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btHfFluid* fluid = new btHfFluid( btScalar(1.0), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		for (int i = 0; i < fluid->getNumNodesLength()*fluid->getNumNodesWidth(); i++) fluid->getEtaArray()[i] = btScalar(0.0f);
		
		fluid->prep();
	}
	world->addHfFluid (fluid);
}
void Run_RandomDrops (btHfFluidRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0f);
	btScalar dt = btScalar(1.0/60.);	
	
	if (dtSinceLastDrop > btScalar(0.5f))
	{
		btAlignedObjectArray<btHfFluid*>& fluids = world->getHfFluidArray();
		btHfFluid* fluid = fluids[0];
		
		dtSinceLastDrop = btScalar(0.0f);
		int randomXNode = GEN_rand () % (fluid->getNumNodesWidth()-2);
		int randomZNode = GEN_rand () % (fluid->getNumNodesLength()-2);
		if (randomXNode <= 1)
			randomXNode = 2;
		if (randomZNode <= 1)
			randomZNode = 2;

		btAlignedObjectArray<btScalar>& eta = fluid->getEtaArray ();
		btAlignedObjectArray<btScalar>& height = fluid->getHeightArray ();
		const btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray ();
		btAlignedObjectArray<bool>& flags = fluid->getFlagsArray();
		int index = fluid->arrayIndex (randomXNode, randomZNode);
		eta[index] += btScalar(4.5f);
		eta[index-1] += btScalar(2.25f);
		eta[index+1] += btScalar(2.25f);
		eta[index+fluid->getNumNodesWidth()] += btScalar(2.25f);
		eta[index-fluid->getNumNodesWidth()] += btScalar(2.25f);
		height[index] = eta[index] + ground[index];
		height[index-1] = eta[index-1] + ground[index-1];
		height[index+1] = eta[index+1] + ground[index+1];
		height[index+fluid->getNumNodesWidth()] = eta[index+fluid->getNumNodesWidth()] + ground[index+fluid->getNumNodesWidth()];
		height[index-fluid->getNumNodesWidth()] = eta[index-fluid->getNumNodesWidth()] + ground[index-fluid->getNumNodesWidth()];
		flags[index] = true;
		flags[index-1] = true;
		flags[index+1] = true;
		flags[index+fluid->getNumNodesWidth()] = true;
		flags[index-fluid->getNumNodesWidth()] = true;
	} 
	else 
	{
		dtSinceLastDrop += dt;
	}
}

void Init_FillPool (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	const int gridLength = 50;
	const int gridWidth = 50;
	btHfFluid* fluid = new btHfFluid( btScalar(1.0), gridLength, gridWidth );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-20.0), btScalar(-5.0), btScalar(-20.0)) );
		fluid->setWorldTransform(transform);
		
		btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray();

		const btScalar poolEdgeHeight = btScalar(10.0f);
		const btScalar poolBottomHeight = btScalar(1.0f);
		const btScalar poolPourerHeight = btScalar(6.0f);
		for (int j = 1; j < fluid->getNumNodesLength()-1; j++)
		{
			for (int i = 1; i < fluid->getNumNodesWidth()-1; i++)
			{
				int index = fluid->arrayIndex (i, j);
				
				// pool edge
				if (j == 1 || i == 1 || j == fluid->getNumNodesLength()-2 || i == fluid->getNumNodesWidth()-2)
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
	world->addHfFluid(fluid);
	
	{
		btConvexShape* colShape = new btBoxShape(btVector3(btScalar(1.), btScalar(1.), btScalar(1.)));
		btHfFluidBuoyantConvexShape* buoyantShape = new btHfFluidBuoyantConvexShape(colShape);
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
					startTransform.setOrigin(btVector3(
										2.0*i + start_x,
										10+2.0*k + start_y,
										2.0*j + start_z));

			
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
void Run_FillPool (btHfFluidRigidDynamicsWorld* world)
{
	btAlignedObjectArray<btHfFluid*>& fluids = world->getHfFluidArray();
	btHfFluid* fluid = fluids[0];

	for (int i = 26; i < 30; i++) fluid->setFluidHeight( fluid->arrayIndex(i, fluid->getNumNodesLength()-3), btScalar(3.0f) );
}


void Init_Fill (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btHfFluid* fluid = new btHfFluid( btScalar(1.0f), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for (int i = 0; i < fluid->getNumNodesLength()*fluid->getNumNodesWidth(); i++) fluid->getEtaArray()[i] = btScalar(0.0f);
		
		fluid->prep ();
	}
	world->addHfFluid (fluid);
}

void Run_Fill (btHfFluidRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0f);
	btScalar dt = btScalar(1.0/60.);

	btAlignedObjectArray<btHfFluid*>& fluids = world->getHfFluidArray ();
	btHfFluid* fluid = fluids[0];

	if (dtSinceLastDrop > btScalar(0.25f))
	{
		dtSinceLastDrop = btScalar(0.0f);

		btAlignedObjectArray<btScalar>& eta = fluid->getEtaArray ();
		btAlignedObjectArray<btScalar>& velocityU = fluid->getVelocityUArray ();
		btAlignedObjectArray<btScalar>& velocityV = fluid->getVelocityVArray ();
		btAlignedObjectArray<btScalar>& height = fluid->getHeightArray ();
		const btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray ();
		btAlignedObjectArray<bool>& flags = fluid->getFlagsArray();
		int index = fluid->arrayIndex (fluid->getNumNodesWidth()/2, fluid->getNumNodesLength()/2);
		eta[index] += btScalar(4.5f);
		eta[index-1] += btScalar(2.25f);
		eta[index+1] += btScalar(2.25f);
		eta[index+fluid->getNumNodesWidth()] += btScalar(2.25f);
		eta[index-fluid->getNumNodesWidth()] += btScalar(2.25f);

		velocityU[index] = btScalar(0.0f);
		velocityU[index-1] = btScalar(-10.0f);
		velocityU[index+1] = btScalar(10.0f);
		velocityU[index+fluid->getNumNodesWidth()] = btScalar(0.0f);
		velocityU[index-fluid->getNumNodesWidth()] = btScalar(0.0f);

		velocityV[index] = btScalar(0.0f);
		velocityV[index-1] = btScalar(0.0f);
		velocityV[index+1] = btScalar(0.0f);
		velocityV[index+fluid->getNumNodesWidth()] = btScalar(10.0f);
		velocityV[index-fluid->getNumNodesWidth()] = btScalar(-10.0f);

		height[index] = eta[index] + ground[index];
		height[index-1] = eta[index-1] + ground[index-1];
		height[index+1] = eta[index+1] + ground[index+1];
		height[index+fluid->getNumNodesWidth()] = eta[index+fluid->getNumNodesWidth()] + ground[index+fluid->getNumNodesWidth()];
		height[index-fluid->getNumNodesWidth()] = eta[index-fluid->getNumNodesWidth()] + ground[index-fluid->getNumNodesWidth()];
		flags[index] = true;
		flags[index-1] = true;
		flags[index+1] = true;
		flags[index+fluid->getNumNodesWidth()] = true;
		flags[index-fluid->getNumNodesWidth()] = true;
	} 
	else 
	{
		dtSinceLastDrop += dt;
	}
	
}

void Init_BlockWave (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btHfFluid* fluid = new btHfFluid (btScalar(1.0), 75, 50);
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);

		btAlignedObjectArray<btScalar>& eta = fluid->getEtaArray();

		for (int i = 0; i < fluid->getNumNodesLength() * fluid->getNumNodesWidth(); i++) eta[i] = btScalar(12.0f);

		for (int i = fluid->getNumNodesWidth()/8; i < fluid->getNumNodesWidth()/4; i++)
		{
			for (int j = fluid->getNumNodesLength()/8; j < fluid->getNumNodesLength()/4; j++)
			{
				int index = fluid->arrayIndex(i, j);
				eta[index] = btScalar(4.0f);
			}
		}
		fluid->prep ();
	}
	world->addHfFluid (fluid);
	
	{
		btConvexShape* colShape = new btSphereShape(btScalar(1.));
		btHfFluidBuoyantConvexShape* buoyantShape = new btHfFluidBuoyantConvexShape(colShape);
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
					startTransform.setOrigin(btVector3(
										2.0*i + start_x,
										10+2.0*k + start_y,
										2.0*j + start_z));

			
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

void Init_Ground (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btHfFluid* fluid = new btHfFluid( btScalar(1.0f), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);

		btAlignedObjectArray<btScalar>& eta = fluid->getEtaArray();

		for (int i = 0; i < fluid->getNumNodesLength() * fluid->getNumNodesWidth(); i++) eta[i] = btScalar(4.0f);

		btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray();
		for (int i = 0; i < fluid->getNumNodesWidth(); i++)
		{
			for (int j = 0; j < fluid->getNumNodesLength(); j++)
			{
				int index = fluid->arrayIndex(i, j);

				if (j <= fluid->getNumNodesLength()/2) ground[index] = btScalar(5.0f);
				else if (j > (fluid->getNumNodesLength()/8*6)) ground[index] = btScalar(0.0f);
				else ground[index] = btScalar(6.5f);
				
				if (j <= fluid->getNumNodesLength()/4 && j > fluid->getNumNodesLength()/8) eta[index] = btScalar(8.0f);
				else if (j <= fluid->getNumNodesLength()/8) eta[index] = btScalar(20.0f);
				else eta[index] = btScalar(0.0f);
			}
		}
		fluid->prep();
	}
	
	world->addHfFluid(fluid);
}

void Init_Ground2 (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btHfFluid* fluid = new btHfFluid (btScalar(1.0f), 75, 50);
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(5.0), btScalar(-50.0)) );
		fluid->setWorldTransform (transform);
		
		btAlignedObjectArray<btScalar>& eta = fluid->getEtaArray();
		btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray();
		
		for (int i = 0; i < fluid->getNumNodesWidth(); i++)
		{
			for (int j = 0; j < fluid->getNumNodesLength(); j++)
			{
				int index = fluid->arrayIndex(i, j);

				ground[index] = (btScalar(j)/fluid->getNumNodesLength()-1)*btScalar(8.0f);
			}
		}

		for (int i = 0; i < fluid->getNumNodesLength() * fluid->getNumNodesWidth(); i++) eta[i] = btScalar(2.0f);

		for (int i = fluid->getNumNodesWidth()/8; i < fluid->getNumNodesWidth()/4; i++)
		{
			for (int j = fluid->getNumNodesLength()/8; j < fluid->getNumNodesLength()/4; j++)
			{
				int index = fluid->arrayIndex(i, j);
				eta[index] = btScalar(8.0f);
			}
		}
		fluid->prep();
	}
	world->addHfFluid(fluid);
}

void Init_Fill2 (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btHfFluid* fluid = new btHfFluid( btScalar(1.0), 100, 100 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for (int i = 0; i < fluid->getNumNodesLength()*fluid->getNumNodesWidth(); i++) fluid->getEtaArray()[i] = btScalar(0.0f);
		
		fluid->prep();
	}
	world->addHfFluid (fluid);
}
void Run_Fill2 (btHfFluidRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0f);
	btScalar dt = btScalar(1.0/60.);

	btAlignedObjectArray<btHfFluid*>& fluids = world->getHfFluidArray ();
	btHfFluid* fluid = fluids[0];

	if (dtSinceLastDrop > btScalar(0.25f))
	{
		dtSinceLastDrop = btScalar(0.0f);

		btAlignedObjectArray<btScalar>& eta = fluid->getEtaArray ();
		btAlignedObjectArray<btScalar>& velocityU = fluid->getVelocityUArray ();
		btAlignedObjectArray<btScalar>& height = fluid->getHeightArray ();
		const btAlignedObjectArray<btScalar>& ground = fluid->getGroundArray ();
		btAlignedObjectArray<bool>& flags = fluid->getFlagsArray();

		for (int i = 1; i < fluid->getNumNodesWidth()-1; i++)
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

void Init_MovingPour (btHfFluidRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btHfFluid* fluid = new btHfFluid( btScalar(1.0), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for (int i = 0; i < fluid->getNumNodesLength()*fluid->getNumNodesWidth(); i++) fluid->getEtaArray()[i] = btScalar(5.0f);
		
		fluid->prep();
	}
	world->addHfFluid(fluid);
	
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
					startTransform.setOrigin(btVector3(
										2.0*i + start_x,
										10+2.0*k + start_y,
										2.0*j + start_z));

			
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

void Run_MovingPour(btHfFluidRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0f);
	static btScalar x = 4;
	static btScalar z = 4;
	static btScalar dx = btScalar(20.0f);
	static btScalar dz = btScalar(30.0f);
	btScalar dt = btScalar(1.0/60.);

	btAlignedObjectArray<btHfFluid*>& fluids = world->getHfFluidArray ();
	btHfFluid* fluid = fluids[0];

	int minX = 2;
	int minZ = 2;
	int maxX = fluid->getNumNodesWidth() - 2;
	int maxZ = fluid->getNumNodesLength() - 2;

	x += dx * dt;
	
	if (x <= minX)
	{
		dx *= btScalar(-1.0f);
		x = minX;
	} 
	else if (x >= maxX) 
	{
		dx *= btScalar(-1.0f);
		x = maxX;
	}
	z += dz * dt;
	
	if (z <= minZ)
	{
		dz *= btScalar(-1.0f);
		z = minZ;
	} 
	else if (z >= maxZ) 
	{
		dz *= btScalar(-1.0f);
		z = maxZ;
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