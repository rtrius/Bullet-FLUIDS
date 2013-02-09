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
#ifndef BT_FLUID_SOFT_RIGID_COLLISION_CONFIGURATION_H
#define BT_FLUID_SOFT_RIGID_COLLISION_CONFIGURATION_H

#include "BulletSoftBody/btSoftBodyRigidBodyCollisionConfiguration.h"


//Sphere-Convex collision is processed by btConvexConvexAlgorithm, 
//which is not optimal for large numbers of SPH particle collisions.
//Enable BT_OVERRIDE_SPHERE_COLLISION to use specialized sphere collision algorithms.
//
//Warning: these algorithms have not been extensively tested(experimental).
//
//#define BT_OVERRIDE_SPHERE_COLLISION


///Preliminary soft body - SPH fluid interaction demo; this class is not supported.
class btFluidSoftRigidCollisionConfiguration : public btSoftBodyRigidBodyCollisionConfiguration
{
	btCollisionAlgorithmCreateFunc*	m_fluidRigidCreateFunc;
	btCollisionAlgorithmCreateFunc*	m_fluidRigidCreateFuncSwapped;

#ifdef BT_OVERRIDE_SPHERE_COLLISION
	btCollisionAlgorithmCreateFunc* m_sphereBoxCF;
	btCollisionAlgorithmCreateFunc* m_boxSphereCF;
	
	btCollisionAlgorithmCreateFunc* m_sphereCapsuleCF;
	btCollisionAlgorithmCreateFunc* m_capsuleSphereCF;
	
	btCollisionAlgorithmCreateFunc* m_sphereConeCF;
	btCollisionAlgorithmCreateFunc* m_coneSphereCF;
	
	btCollisionAlgorithmCreateFunc* m_sphereCylinderCF;
	btCollisionAlgorithmCreateFunc* m_cylinderSphereCF;
	
	btCollisionAlgorithmCreateFunc* m_sphereHeightfieldCF;
	btCollisionAlgorithmCreateFunc* m_heightfieldSphereCF;
#endif //BT_OVERRIDE_SPHERE_COLLISION
	
public:
	btFluidSoftRigidCollisionConfiguration( const btDefaultCollisionConstructionInfo& constructionInfo = btDefaultCollisionConstructionInfo() );
	virtual ~btFluidSoftRigidCollisionConfiguration();
	
	virtual btCollisionAlgorithmCreateFunc* getCollisionAlgorithmCreateFunc(int proxyType0, int proxyType1);
};

#endif
