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

inline btRigidBody* createRigidBody(const btTransform& transform, btScalar mass, btCollisionShape* shape)
{
	//Rigid bodies are dynamic if and only if mass is non zero, otherwise static
	btVector3 localInertia(0,0,0);
	if(mass != 0.f) shape->calculateLocalInertia(mass, localInertia);
	
	btDefaultMotionState* motionState = new btDefaultMotionState(transform);
	btRigidBody::btRigidBodyConstructionInfo rbInfo(mass, motionState, shape, localInertia);
	btRigidBody* rigid = new btRigidBody(rbInfo);
		
	return rigid;
}

void Init_Floatyness(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(0.25), 100, 100 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-10.0), btScalar(-5.0), btScalar(-10.0)) );
		
		fluid->setWorldTransform(transform);
		
		btFluidHfParameters HP = fluid->getParameters();
		HP.m_volumeDisplacementScale = btScalar(0.0);
		HP.m_horizontalVelocityScale = btScalar(0.0);
		fluid->setParameters(HP);

		for(int i = 0; i < fluid->getNumNodes(); i++) fluid->setFluidHeight( i, btScalar(5.0) );
		fluid->prep();
	}
	world->addFluidHf(fluid);
	
	btConvexShape* sphereShape = new btSphereShape( btScalar(1.0) );
	
	
	const int NUM_OBJECTS = 5;
	btScalar floatyness = btScalar(1.0);
	btScalar floatynessDelta = btScalar(0.25);
	btScalar start_x = btScalar(-5.0);
	btScalar start_y = btScalar(7.5);
	btScalar start_z = btScalar(-5.0);
	btScalar step_x = btScalar(3.0);
	
	for(int j = 0; j < 2; ++j)
	{
		if(j >= 1)
		{
			floatyness = btScalar(2.0);
			start_y = btScalar(-4.0);
			start_z = btScalar(5.0);
		}
		
		for(int i = 0; i < NUM_OBJECTS; i++)
		{
			btFluidHfBuoyantConvexShape* buoyantShape = new btFluidHfBuoyantConvexShape(sphereShape);
			buoyantShape->generateShape( btScalar(0.25), btScalar(0.05) );
			buoyantShape->setFloatyness(floatyness + floatynessDelta * i);
			collisionShapes.push_back(buoyantShape);
			
			btTransform transform( btQuaternion::getIdentity(), btVector3(step_x * i + start_x, start_y, start_z) );
			
			const btScalar MASS = btScalar(1.0);
			btRigidBody* body = createRigidBody(transform, MASS, buoyantShape);
			world->addRigidBody(body);
		}
	}
}

void Init_Bowl(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), 50, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-10.0), btScalar(-5.0), btScalar(-10.0)) );
		fluid->setWorldTransform(transform);
		
		btScalar amplitude = btScalar(200.0);
		for(int i = 0; i < fluid->getNumNodesX(); i++)
		{
			btScalar x = btScalar(i - fluid->getNumNodesX()/2) / btScalar(fluid->getNumNodesX()*2);
			btScalar xh = amplitude * (x * x) + btScalar(5.0);
			for(int j = 0; j < fluid->getNumNodesZ(); j++)
			{
				btScalar y = btScalar(j - fluid->getNumNodesZ()/2) / btScalar(fluid->getNumNodesZ()*2);
				btScalar yh = amplitude * (y * y) + btScalar(5.0);
				
				btScalar groundHeight = btMax(xh,yh);
				
				btScalar waterHeight = btScalar(0.0);
				if(groundHeight > 14.0) waterHeight = btScalar(0.0); 
				else waterHeight = btScalar(14.0) - groundHeight;
				
				int index = fluid->arrayIndex(i, j);
				fluid->setGroundHeight(index, groundHeight);
				fluid->setFluidHeight(index, waterHeight);
			}
		}

		fluid->prep();
	}
	world->addFluidHf(fluid);
	
	//Create a few dynamic rigidbodies
	{
		btConvexShape* colShape = new btBoxShape( btVector3(1,1,1) );
		btFluidHfBuoyantConvexShape* buoyantShape = new btFluidHfBuoyantConvexShape(colShape);
		buoyantShape->generateShape( btScalar(0.25), btScalar(0.05) );
		collisionShapes.push_back(colShape);
		collisionShapes.push_back(buoyantShape);
		
		for(int k=0;k<ARRAY_SIZE_Y;k++)
		{
			for(int i=0;i<ARRAY_SIZE_X;i++)
			{
				for(int j = 0;j<ARRAY_SIZE_Z;j++)
				{
					const btVector3 START( btScalar(5.0 - ARRAY_SIZE_X/2), btScalar(5.0), btScalar(3.0 - ARRAY_SIZE_Z/2) );
					btVector3 position = START + btVector3( btScalar(i), btScalar(k), btScalar(j) ) * btScalar(2.0);
					
					btTransform transform( btQuaternion::getIdentity(), position );
					
					const btScalar MASS(1.0);
					btRigidBody* body = createRigidBody(transform, MASS, buoyantShape);
					world->addRigidBody(body);
				}
			}
		}
	}
}

void Init_Drops(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(0.5), 50, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-10.0), btScalar(-5.0), btScalar(-10.0)) );
		fluid->setWorldTransform(transform);

		for(int i = 0; i < fluid->getNumNodes(); i++) fluid->setFluidHeight(i, btScalar(5.0));
		fluid->prep();
	}
	world->addFluidHf(fluid);
	
	//Create a few dynamic rigidbodies
	{
		btConvexShape* colShape = new btBoxShape( btVector3(5, 0.5, 5) );
		btFluidHfBuoyantConvexShape* buoyantShape = new btFluidHfBuoyantConvexShape(colShape);
		buoyantShape->generateShape( btScalar(0.25), btScalar(0.05) );
		
		collisionShapes.push_back(colShape);
		collisionShapes.push_back(buoyantShape);
		
		for(int k=0;k<ARRAY_SIZE_Y;k++)
		{
			for(int i=0;i<ARRAY_SIZE_X;i++)
			{
				for(int j = 0;j<ARRAY_SIZE_Z;j++)
				{
					const btVector3 START( btScalar(5.0 - ARRAY_SIZE_X/2), btScalar(5.0), btScalar(3.0 - ARRAY_SIZE_Z/2) );
					btVector3 position = START + btVector3( btScalar(i), btScalar(k), btScalar(j) ) * btScalar(2.0);
					
					btTransform transform( btQuaternion::getIdentity(), position );
					
					const btScalar MASS(1.0);
					btRigidBody* body = createRigidBody(transform, MASS, buoyantShape);
					world->addRigidBody(body);
				}
			}
		}
	}
}

void Init_Wave(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for(int i = 0; i < fluid->getNumNodes(); i++) fluid->setFluidHeight( i, btScalar(10.0) );

		int midZ = fluid->getNumNodesZ() / 2;
		for(int i = 1; i < fluid->getNumNodesX()-1; i++)
		{
			fluid->setFluidHeight( fluid->arrayIndex(i, midZ - 1), btScalar(2.0) );
			fluid->setFluidHeight( fluid->arrayIndex(i, midZ), btScalar(2.0) );
			fluid->setFluidHeight( fluid->arrayIndex(i, midZ + 1), btScalar(2.0) );
		}

		fluid->prep();
	}
	world->addFluidHf(fluid);
}

void Init_RandomDrops(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		for(int i = 0; i < fluid->getNumNodes(); i++) fluid->setFluidHeight( i, btScalar(0.0) );
		
		fluid->prep();
	}
	world->addFluidHf(fluid);
}
void Run_RandomDrops(btFluidHfRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0);
	btScalar dt = btScalar(1.0/60.0);	
	
	if( dtSinceLastDrop > btScalar(0.5) )
	{
		dtSinceLastDrop = btScalar(0.0);
		
		btAlignedObjectArray<btFluidHf*>& fluids = world->getFluidHfArray();
		btFluidHf* fluid = fluids[0];
		
		int randomXNode = GEN_rand() % ( fluid->getNumNodesX()-2 );
		int randomZNode = GEN_rand() % ( fluid->getNumNodesZ()-2 );
		if(randomXNode <= 1) randomXNode = 2;
		if(randomZNode <= 1) randomZNode = 2;

		int index = fluid->arrayIndex(randomXNode, randomZNode);
		
		fluid->addFluidHeight( index, btScalar(4.5) );
		fluid->addFluidHeight( index-1, btScalar(2.25) );
		fluid->addFluidHeight( index+1, btScalar(2.25) );
		fluid->addFluidHeight( index-fluid->getNumNodesX(), btScalar(2.25) );
		fluid->addFluidHeight( index+fluid->getNumNodesX(), btScalar(2.25) );
	} 
	else dtSinceLastDrop += dt;
}

void Init_FillPool(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	const int gridLength = 50;
	const int gridWidth = 50;
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), gridLength, gridWidth );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-20.0), btScalar(-5.0), btScalar(-20.0)) );
		fluid->setWorldTransform(transform);
		
		const btScalar poolEdgeHeight = btScalar(10.0);
		const btScalar poolBottomHeight = btScalar(1.0);
		const btScalar poolPourerHeight = btScalar(6.0);
		for(int j = 1; j < fluid->getNumNodesZ()-1; j++)
		{
			for(int i = 1; i < fluid->getNumNodesX()-1; i++)
			{
				int index = fluid->arrayIndex (i, j);
				
				// pool edge
				//if (j == 1 || i == 1 || j == fluid->getNumNodesZ()-2 || i == fluid->getNumNodesX()-2)
				if (j <= 1+3 || i <= 1+8 || j >= fluid->getNumNodesZ()-2 || i >= fluid->getNumNodesX()-2-3)
				{
					fluid->setGroundHeight(index, poolEdgeHeight);
					continue;
				}
				if (j > 35)
				{
					if (i <= 25 || i >= 30) fluid->setGroundHeight(index, poolEdgeHeight);
					else fluid->setGroundHeight(index, poolPourerHeight);
					
					continue;
				}
				fluid->setGroundHeight(index, poolBottomHeight);
			}
		}
		
		fluid->prep();
	}
	world->addFluidHf(fluid);
	
	const int SPAWN_BOXES = 0;
	if(SPAWN_BOXES)
	{
		btConvexShape* colShape = new btBoxShape( btVector3(1, 1, 1) );
		btFluidHfBuoyantConvexShape* buoyantShape = new btFluidHfBuoyantConvexShape(colShape);
		buoyantShape->generateShape( btScalar(0.25), btScalar(0.05) );
		collisionShapes.push_back(colShape);
		collisionShapes.push_back(buoyantShape);

		/// Create Dynamic Objects
		const int GRID_SIZE = 2;
		for(int k=0;k<GRID_SIZE;k++)
		{
			for(int i=0;i<GRID_SIZE;i++)
			{
				for(int j = 0;j<GRID_SIZE;j++)
				{
					const btVector3 START( btScalar(-10.0 - GRID_SIZE/2), btScalar(2.0), btScalar(-10.0 - GRID_SIZE/2) );
					btVector3 position = START + btVector3( btScalar(i), btScalar(k), btScalar(j) ) * btScalar(2.0);

					btTransform transform( btQuaternion::getIdentity(), position );
					
					const btScalar MASS(1.0);
					btRigidBody* body = createRigidBody(transform, MASS, buoyantShape);
					world->addRigidBody(body);
				}
			}
		}
	}
}
void Run_FillPool(btFluidHfRigidDynamicsWorld* world)
{
	btAlignedObjectArray<btFluidHf*>& fluids = world->getFluidHfArray();
	btFluidHf* fluid = fluids[0];

	for(int i = 26; i < 30; i++) fluid->setFluidHeight( fluid->arrayIndex(i, fluid->getNumNodesZ()-3), btScalar(3.0) );
	
	const int DEBUG = 0;
	if(DEBUG)
	{
		for(int j = 1; j < fluid->getNumNodesZ()-1; j++)
			for(int i = 1; i < fluid->getNumNodesX()-1; i++)
			{
				if(j <= 1+15) 
				{
					int index = fluid->arrayIndex(i, j);
					fluid->setVelocity( index, btScalar(0.0), btScalar(0.0) );
				}
			}
	}
}


void Init_Fill(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for(int i = 0; i < fluid->getNumNodes(); i++) fluid->setFluidHeight( i, btScalar(0.0) );
		
		fluid->prep ();
	}
	world->addFluidHf(fluid);
}

void Run_Fill(btFluidHfRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0);
	btScalar dt = btScalar(1.0/60.0);

	btAlignedObjectArray<btFluidHf*>& fluids = world->getFluidHfArray();
	btFluidHf* fluid = fluids[0];

	if( dtSinceLastDrop > btScalar(0.25) )
	{
		dtSinceLastDrop = btScalar(0.0);

		int index = fluid->arrayIndex( fluid->getNumNodesX()/2, fluid->getNumNodesZ()/2 );
		
		fluid->addFluidHeight( index, btScalar(4.5) );
		fluid->addFluidHeight( index-1, btScalar(2.25) );
		fluid->addFluidHeight( index+1, btScalar(2.25) );
		fluid->addFluidHeight( index-fluid->getNumNodesX(), btScalar(2.25) );
		fluid->addFluidHeight( index+fluid->getNumNodesX(), btScalar(2.25) );
		
		fluid->setVelocity( index, btScalar(0.0), btScalar(0.0) );
		fluid->setVelocity( index-1, btScalar(-10.0), btScalar(0.0) );
		fluid->setVelocity( index+1, btScalar(10.0), btScalar(0.0) );
		fluid->setVelocity( index-fluid->getNumNodesX(), btScalar(0.0), btScalar(-10.0) );
		fluid->setVelocity( index+fluid->getNumNodesX(), btScalar(0.0), btScalar(10.0) );
	} 
	else dtSinceLastDrop += dt;
	
}

void Init_BlockWave(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf(btScalar(1.0), 75, 50);
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for(int i = 0; i < fluid->getNumNodes(); i++) fluid->setFluidHeight( i, btScalar(12.0) );

		for(int i = fluid->getNumNodesX()/8; i < fluid->getNumNodesX()/4; i++)
		{
			for(int j = fluid->getNumNodesZ()/8; j < fluid->getNumNodesZ()/4; j++)
			{
				int index = fluid->arrayIndex(i, j);
				fluid->setFluidHeight( index, btScalar(4.0) );
			}
		}
		fluid->prep ();
	}
	world->addFluidHf(fluid);
	
	{
		btConvexShape* colShape = new btSphereShape( btScalar(1.0) );
		btFluidHfBuoyantConvexShape* buoyantShape = new btFluidHfBuoyantConvexShape(colShape);
		buoyantShape->generateShape( btScalar(0.25), btScalar(0.05) );
		collisionShapes.push_back(buoyantShape);
		collisionShapes.push_back(colShape);
		
		const int GRID_SIZE = 2;
		for(int k=0;k<GRID_SIZE;k++)
		{
			for(int i=0;i<GRID_SIZE;i++)
			{
				for(int j = 0;j<GRID_SIZE;j++)
				{
					const btVector3 START( btScalar(-10.0 - GRID_SIZE/2), btScalar(2.0), btScalar(-10.0 - GRID_SIZE/2) );
					btVector3 position = START + btVector3( btScalar(i), btScalar(k), btScalar(j) ) * btScalar(2.0);
					
					btTransform transform( btQuaternion::getIdentity(), position );
					
					const btScalar MASS(1.0);
					btRigidBody* body = createRigidBody(transform, MASS, buoyantShape);
					world->addRigidBody(body);
				}
			}
		}
	}
}

void Init_Ground(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);

		for(int i = 0; i < fluid->getNumNodes(); i++) fluid->setFluidHeight( i, btScalar(4.0) );

		for(int i = 0; i < fluid->getNumNodesX(); i++)
		{
			for(int j = 0; j < fluid->getNumNodesZ(); j++)
			{
				int index = fluid->arrayIndex(i, j);

				btScalar groundHeight;
				if (j <= fluid->getNumNodesZ()/2) groundHeight = btScalar(5.0);
				else if (j > (fluid->getNumNodesZ()/8*6)) groundHeight = btScalar(0.0);
				else groundHeight = btScalar(6.5);
				
				btScalar fluidHeight;
				if (j <= fluid->getNumNodesZ()/4 && j > fluid->getNumNodesZ()/8) fluidHeight = btScalar(8.0);
				else if (j <= fluid->getNumNodesZ()/8) fluidHeight = btScalar(20.0);
				else fluidHeight = btScalar(0.0);
				
				fluid->setGroundHeight(index, groundHeight);
				fluid->setFluidHeight(index, fluidHeight);
			}
		}
		fluid->prep();
	}
	
	world->addFluidHf(fluid);
}

void Init_Ground2(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf(btScalar(1.0), 75, 50);
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for(int i = 0; i < fluid->getNumNodesX(); i++)
		{
			for(int j = 0; j < fluid->getNumNodesZ(); j++)
			{
				int index = fluid->arrayIndex(i, j);

				fluid->setGroundHeight( index, (btScalar(j)/fluid->getNumNodesZ()-1)*btScalar(8.0) );
			}
		}

		for(int i = 0; i < fluid->getNumNodes(); i++) fluid->setFluidHeight( i, btScalar(2.0) );

		for(int i = fluid->getNumNodesX()/8; i < fluid->getNumNodesX()/4; i++)
		{
			for(int j = fluid->getNumNodesZ()/8; j < fluid->getNumNodesZ()/4; j++)
			{
				int index = fluid->arrayIndex(i, j);
				fluid->setFluidHeight( index, btScalar(8.0) );
			}
		}
		fluid->prep();
	}
	world->addFluidHf(fluid);
}

void Init_Fill2(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), 100, 100 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for(int i = 0; i < fluid->getNumNodes(); i++) fluid->setFluidHeight( i, btScalar(0.0) );
		
		fluid->prep();
	}
	world->addFluidHf(fluid);
}
void Run_Fill2(btFluidHfRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0);
	btScalar dt = btScalar(1.0/60.0);

	btAlignedObjectArray<btFluidHf*>& fluids = world->getFluidHfArray();
	btFluidHf* fluid = fluids[0];

	if( dtSinceLastDrop > btScalar(0.25) )
	{
		dtSinceLastDrop = btScalar(0.0);

		for(int i = 1; i < fluid->getNumNodesX()-1; i++)
		{
			int index = fluid->arrayIndex(i, 1);
			fluid->setVelocityX( index, btScalar(4.0) );
			fluid->addFluidHeight( index, btScalar(3.0) );
		}
	} 
	else dtSinceLastDrop += dt;
}

void Init_MovingPour(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes)
{
	btFluidHf* fluid = new btFluidHf( btScalar(1.0), 75, 50 );
	{
		btTransform transform( btQuaternion::getIdentity(), btVector3(btScalar(-50.0), btScalar(-5.0), btScalar(-50.0)) );
		fluid->setWorldTransform(transform);
		
		for(int i = 0; i < fluid->getNumNodes(); i++) fluid->setFluidHeight( i, btScalar(5.0) );
		
		fluid->prep();
	}
	world->addFluidHf(fluid);
	
	//Create a few dynamic rigidbodies
	{
		btCollisionShape* boxShape = new btBoxShape( btVector3(1,1,1) );
		collisionShapes.push_back(boxShape);
		
		for(int k=0;k<ARRAY_SIZE_Y;k++)
		{
			for(int i=0;i<ARRAY_SIZE_X;i++)
			{
				for(int j = 0;j<ARRAY_SIZE_Z;j++)
				{
					const btVector3 START( btScalar(5.0 - ARRAY_SIZE_X/2), btScalar(5.0), btScalar(3.0 - ARRAY_SIZE_Z/2) );
					btVector3 position = START + btVector3( btScalar(i), btScalar(k), btScalar(j) ) * btScalar(2.0);

					btTransform transform( btQuaternion::getIdentity(), position );
					
					const btScalar MASS(1.f);
					btRigidBody* body = createRigidBody(transform, MASS, boxShape);
					world->addRigidBody(body);
				}
			}
		}
	}
}

void Run_MovingPour(btFluidHfRigidDynamicsWorld* world)
{
	static btScalar dtSinceLastDrop = btScalar(0.0);
	static btScalar x = 4;
	static btScalar z = 4;
	static btScalar dx = btScalar(20.0);
	static btScalar dz = btScalar(30.0);
	btScalar dt = btScalar(1.0/60.);

	btAlignedObjectArray<btFluidHf*>& fluids = world->getFluidHfArray();
	btFluidHf* fluid = fluids[0];

	int minX = 2;
	int minZ = 2;
	int maxX = fluid->getNumNodesX() - 2;
	int maxZ = fluid->getNumNodesZ() - 2;

	x += dx * dt;
	if (x <= minX)
	{
		dx *= btScalar(-1.0);
		x = static_cast<btScalar>(minX);
	} 
	else if (x >= maxX) 
	{
		dx *= btScalar(-1.0);
		x = static_cast<btScalar>(maxX);
	}
	
	z += dz * dt;
	if (z <= minZ)
	{
		dz *= btScalar(-1.0);
		z = static_cast<btScalar>(minZ);
	} 
	else if (z >= maxZ) 
	{
		dz *= btScalar(-1.0);
		z = static_cast<btScalar>(maxZ);
	}

	const btScalar DROP_HEIGHT = btScalar(3.0);
	fluid->addFluidHeight( fluid->arrayIndex( static_cast<int>(x), static_cast<int>(z) ), DROP_HEIGHT );
}