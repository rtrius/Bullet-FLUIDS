/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. 
   If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_FLUID_RIGID_COLLISION_CONFIGURATION_H
#define BT_FLUID_RIGID_COLLISION_CONFIGURATION_H

#include "LinearMath/btPoolAllocator.h"
#include "BulletCollision/CollisionDispatch/btDefaultCollisionConfiguration.h"

#include "btFluidSphRigidCollisionAlgorithm.h"


///Includes btFluidSph collision support on top of btDefaultCollisionConfiguration.
class btFluidRigidCollisionConfiguration : public btDefaultCollisionConfiguration
{
	btCollisionAlgorithmCreateFunc*	m_fluidRigidConvexCreateFunc;
	btCollisionAlgorithmCreateFunc*	m_fluidRigidConvexCreateFuncSwapped;

public:

	btFluidRigidCollisionConfiguration( const btDefaultCollisionConstructionInfo& constructionInfo = btDefaultCollisionConstructionInfo() )
	: btDefaultCollisionConfiguration(constructionInfo)
	{
		void* ptr;

		ptr = btAlignedAlloc( sizeof(btFluidSphRigidCollisionAlgorithm::CreateFunc), 16 );
		m_fluidRigidConvexCreateFunc = new(ptr) btFluidSphRigidCollisionAlgorithm::CreateFunc;
		
		ptr = btAlignedAlloc( sizeof(btFluidSphRigidCollisionAlgorithm::CreateFunc), 16 );
		m_fluidRigidConvexCreateFuncSwapped = new(ptr) btFluidSphRigidCollisionAlgorithm::CreateFunc;
		m_fluidRigidConvexCreateFuncSwapped->m_swapped = true;
		
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

	virtual ~btFluidRigidCollisionConfiguration()
	{
		m_fluidRigidConvexCreateFunc->~btCollisionAlgorithmCreateFunc();
		btAlignedFree(m_fluidRigidConvexCreateFunc);

		m_fluidRigidConvexCreateFuncSwapped->~btCollisionAlgorithmCreateFunc();
		btAlignedFree(m_fluidRigidConvexCreateFuncSwapped);
	}

	virtual btCollisionAlgorithmCreateFunc* getCollisionAlgorithmCreateFunc(int proxyType0, int proxyType1)
	{
		if( proxyType0 == CUSTOM_CONCAVE_SHAPE_TYPE && btBroadphaseProxy::isConvex(proxyType1) ) return	m_fluidRigidConvexCreateFunc;

		if( btBroadphaseProxy::isConvex(proxyType0) && proxyType1 == CUSTOM_CONCAVE_SHAPE_TYPE ) return	m_fluidRigidConvexCreateFuncSwapped;

		return btDefaultCollisionConfiguration::getCollisionAlgorithmCreateFunc(proxyType0, proxyType1);
	}

};
#endif