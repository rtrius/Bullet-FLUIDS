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

#ifndef HF_DEMOS_H
#define HF_DEMOS_H

#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btScalar.h"

class btCollisionShape;
class btFluidHfRigidDynamicsWorld;

void Init_Floatyness(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
void Init_Bowl(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
void Init_FillPool(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
void Init_Drops(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
void Init_Wave(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
void Init_RandomDrops(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
void Init_Fill(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
void Init_Fill2(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
void Init_BlockWave(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
void Init_Ground(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
void Init_Ground2(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
void Init_MovingPour(btFluidHfRigidDynamicsWorld* world, btAlignedObjectArray<btCollisionShape*>& collisionShapes);
	
void Run_FillPool(btFluidHfRigidDynamicsWorld* world);
void Run_RandomDrops(btFluidHfRigidDynamicsWorld* world);
void Run_Fill(btFluidHfRigidDynamicsWorld* world);
void Run_Fill2(btFluidHfRigidDynamicsWorld* world);
void Run_MovingPour(btFluidHfRigidDynamicsWorld* world);

#define NUM_DEMOS 12
static void (*demo_init_functions[NUM_DEMOS])(btFluidHfRigidDynamicsWorld*, btAlignedObjectArray<btCollisionShape*>&) =
{
	Init_Floatyness,
	Init_Bowl,
	Init_FillPool,
	Init_Drops,
	Init_Wave,
	Init_RandomDrops,
	Init_Fill,
	Init_Fill2,
	Init_BlockWave,
	Init_Ground,
	Init_Ground2,
	Init_MovingPour,
};
static void (*demo_run_functions[NUM_DEMOS])(btFluidHfRigidDynamicsWorld*) =
{
	NULL, 			// Run_Floatyness
	NULL,			// Run_Bowl
	Run_FillPool, 	// Run_FillPool
	NULL, 			// Run_Drops
	NULL, 			// Run_Wave
	Run_RandomDrops,
	Run_Fill,
	Run_Fill2,
	NULL, 			// Run_BlockWave
	NULL, 			// Run_Ground
	NULL, 			// Run_Ground2
	Run_MovingPour,
};

//Camera direction
const btScalar g_ele_array[NUM_DEMOS] = 
{
	btScalar(10),
	btScalar(45),
	btScalar(35),
	btScalar(35),
	btScalar(10),
	btScalar(10),
	btScalar(35),
	btScalar(45),
	btScalar(35),
	btScalar(20),
	btScalar(20),
};

const btScalar g_azi_array[NUM_DEMOS] = 
{
	btScalar(0),
	btScalar(55),
	btScalar(245),
	btScalar(270),
	btScalar(55),
	btScalar(55),
	btScalar(180),
	btScalar(205),
	btScalar(255),
	btScalar(305),
	btScalar(305),
};
const btScalar g_cameraDistance_array[NUM_DEMOS] = 
{
	btScalar(20),
	btScalar(29),
	btScalar(43),
	btScalar(26),
	btScalar(77),
	btScalar(77),
	btScalar(77),
	btScalar(32),
	btScalar(62),
	btScalar(70),
	btScalar(70),
};


#endif