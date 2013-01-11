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
#include "btFluidRigidCollisionConfiguration.h"

#include "LinearMath/btPoolAllocator.h"
#include "BulletCollision/BroadphaseCollision/btBroadphaseProxy.h"

#ifdef BT_OVERRIDE_SPHERE_COLLISION
	#include "BulletCollision/CollisionDispatch/btSphereBoxCollisionAlgorithm.h"
	#include "BulletFluids/CollisionDispatch_temp/btSphereCapsuleCollisionAlgorithm.h"	//Move to BulletCollision/CollisionDispatch
#endif //BT_OVERRIDE_SPHERE_COLLISION

#include "Sph/btFluidSphRigidCollisionAlgorithm.h"


btFluidRigidCollisionConfiguration::btFluidRigidCollisionConfiguration(const btDefaultCollisionConstructionInfo& constructionInfo)
: btDefaultCollisionConfiguration(constructionInfo)
{
	void* ptr;

#ifdef BT_OVERRIDE_SPHERE_COLLISION
	ptr = btAlignedAlloc( sizeof(btSphereBoxCollisionAlgorithm::CreateFunc), 16 );
	m_sphereBoxCF = new(ptr) btSphereBoxCollisionAlgorithm::CreateFunc;
	
	ptr = btAlignedAlloc( sizeof(btSphereBoxCollisionAlgorithm::CreateFunc), 16 );
	m_boxSphereCF = new(ptr) btSphereBoxCollisionAlgorithm::CreateFunc;
	m_boxSphereCF->m_swapped = true;
	
	
	ptr = btAlignedAlloc( sizeof(btSphereCapsuleCollisionAlgorithm::CreateFunc), 16 );
	m_sphereCapsuleCF = new(ptr) btSphereCapsuleCollisionAlgorithm::CreateFunc;
	
	ptr = btAlignedAlloc( sizeof(btSphereCapsuleCollisionAlgorithm::CreateFunc), 16 );
	m_capsuleSphereCF = new(ptr) btSphereCapsuleCollisionAlgorithm::CreateFunc;
	m_capsuleSphereCF->m_swapped = true;
#endif //BT_OVERRIDE_SPHERE_COLLISION
	
	
	ptr = btAlignedAlloc( sizeof(btFluidSphRigidCollisionAlgorithm::CreateFunc), 16 );
	m_fluidRigidCreateFunc = new(ptr) btFluidSphRigidCollisionAlgorithm::CreateFunc;
	
	ptr = btAlignedAlloc( sizeof(btFluidSphRigidCollisionAlgorithm::CreateFunc), 16 );
	m_fluidRigidCreateFuncSwapped = new(ptr) btFluidSphRigidCollisionAlgorithm::CreateFunc;
	m_fluidRigidCreateFuncSwapped->m_swapped = true;
	
	//Collision algorithms introducted by btFluidRigidCollisionConfiguration may be
	//larger than m_collisionAlgorithmPool's element size. Resize if it is not large enough.
	int maxAlgorithmSize = sizeof(btFluidSphRigidCollisionAlgorithm);
	if( m_ownsCollisionAlgorithmPool && m_collisionAlgorithmPool && maxAlgorithmSize > m_collisionAlgorithmPool->getElementSize() )
	{
		m_collisionAlgorithmPool->~btPoolAllocator();
		btAlignedFree(m_collisionAlgorithmPool);
		
		ptr = btAlignedAlloc( sizeof(btPoolAllocator), 16 );
		m_collisionAlgorithmPool = new(ptr) btPoolAllocator(maxAlgorithmSize, constructionInfo.m_defaultMaxCollisionAlgorithmPoolSize);
	}
}

btFluidRigidCollisionConfiguration::~btFluidRigidCollisionConfiguration()
{
#ifdef BT_OVERRIDE_SPHERE_COLLISION
	m_sphereBoxCF->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_sphereBoxCF);
	m_boxSphereCF->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_boxSphereCF);
	
	m_sphereCapsuleCF->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_sphereCapsuleCF);
	m_capsuleSphereCF->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_capsuleSphereCF);
#endif //BT_OVERRIDE_SPHERE_COLLISION

	m_fluidRigidCreateFunc->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_fluidRigidCreateFunc);

	m_fluidRigidCreateFuncSwapped->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_fluidRigidCreateFuncSwapped);
}

btCollisionAlgorithmCreateFunc* btFluidRigidCollisionConfiguration::getCollisionAlgorithmCreateFunc(int proxyType0, int proxyType1)
{
#ifdef BT_OVERRIDE_SPHERE_COLLISION
	if(proxyType0 == SPHERE_SHAPE_PROXYTYPE && proxyType1 == BOX_SHAPE_PROXYTYPE) return m_sphereBoxCF;
	if(proxyType0 == BOX_SHAPE_PROXYTYPE && proxyType1 == SPHERE_SHAPE_PROXYTYPE) return m_boxSphereCF;
	
	if(proxyType0 == SPHERE_SHAPE_PROXYTYPE && proxyType1 == CAPSULE_SHAPE_PROXYTYPE) return m_sphereCapsuleCF;
	if(proxyType0 == CAPSULE_SHAPE_PROXYTYPE && proxyType1 == SPHERE_SHAPE_PROXYTYPE) return m_capsuleSphereCF;
#endif //BT_OVERRIDE_SPHERE_COLLISION

	//	btFluidSph-btSoftBody interaction is not implemented
	//	temporarily use SOFTBODY_SHAPE_PROXYTYPE (replace later with FLUID_SPH_SHAPE_PROXYTYPE)

	bool collideProxyType1 = ( btBroadphaseProxy::isConvex(proxyType1)
								|| btBroadphaseProxy::isConcave(proxyType1)
								|| btBroadphaseProxy::isCompound(proxyType1) );

	if(proxyType0 == SOFTBODY_SHAPE_PROXYTYPE  && collideProxyType1) return m_fluidRigidCreateFunc;

	bool collideProxyType0 = ( btBroadphaseProxy::isConvex(proxyType0)
								|| btBroadphaseProxy::isConcave(proxyType0)
								|| btBroadphaseProxy::isCompound(proxyType0) );
	
	if(collideProxyType0 && proxyType1 == SOFTBODY_SHAPE_PROXYTYPE ) return m_fluidRigidCreateFuncSwapped;

	return btDefaultCollisionConfiguration::getCollisionAlgorithmCreateFunc(proxyType0, proxyType1);
}