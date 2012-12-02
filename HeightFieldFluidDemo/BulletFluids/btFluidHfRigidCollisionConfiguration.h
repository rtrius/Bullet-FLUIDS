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

#ifndef BT_FLUID_HF_RIGID_COLLISION_CONFIGURATION
#define BT_FLUID_HF_RIGID_COLLISION_CONFIGURATION

#include "BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h"

///btFluidHfRigidCollisionConfiguration add btFluidHf interaction on top of btDefaultCollisionConfiguration
class	btFluidHfRigidCollisionConfiguration : public btDefaultCollisionConfiguration
{
	//default CreationFunctions, filling the m_doubleDispatch table
	btCollisionAlgorithmCreateFunc*	m_fluidHfRigidConvexCreateFunc;
	btCollisionAlgorithmCreateFunc*	m_fluidHfRigidConvexCreateFuncSwapped;
	btCollisionAlgorithmCreateFunc*	m_fluidHfBuoyantShapeCollisionCreateFunc;

public:

	btFluidHfRigidCollisionConfiguration(const btDefaultCollisionConstructionInfo& constructionInfo = btDefaultCollisionConstructionInfo());

	virtual ~btFluidHfRigidCollisionConfiguration();

	///creation of soft-soft and soft-rigid, and otherwise fallback to base class implementation
	virtual btCollisionAlgorithmCreateFunc* getCollisionAlgorithmCreateFunc(int proxyType0,int proxyType1);
};

#endif //BT_FLUID_HF_RIGID_COLLISION_CONFIGURATION

