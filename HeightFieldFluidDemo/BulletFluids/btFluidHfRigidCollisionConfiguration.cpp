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

#include "btFluidHfRigidCollisionConfiguration.h"

#include "LinearMath/btPoolAllocator.h"

#include "Hf/btFluidHfRigidCollisionAlgorithm.h"
#include "Hf/btFluidHfBuoyantShapeCollisionAlgorithm.h"

btFluidHfRigidCollisionConfiguration::btFluidHfRigidCollisionConfiguration(const btDefaultCollisionConstructionInfo& constructionInfo)
:btDefaultCollisionConfiguration(constructionInfo)
{
	void* mem;

	mem = btAlignedAlloc(sizeof(btFluidHfRigidCollisionAlgorithm::CreateFunc),16);
	m_fluidHfRigidConvexCreateFunc = new(mem) btFluidHfRigidCollisionAlgorithm::CreateFunc;

	mem = btAlignedAlloc(sizeof(btFluidHfRigidCollisionAlgorithm::CreateFunc),16);
	m_fluidHfRigidConvexCreateFuncSwapped = new(mem) btFluidHfRigidCollisionAlgorithm::CreateFunc;
	m_fluidHfRigidConvexCreateFuncSwapped->m_swapped = true;

	mem = btAlignedAlloc(sizeof(btFluidHfBuoyantShapeCollisionAlgorithm::CreateFunc),16);
	m_fluidHfBuoyantShapeCollisionCreateFunc = new(mem) btFluidHfBuoyantShapeCollisionAlgorithm::CreateFunc(m_simplexSolver, m_pdSolver);

	if (m_ownsCollisionAlgorithmPool && m_collisionAlgorithmPool)
	{
		int curElemSize = m_collisionAlgorithmPool->getElementSize();
		///calculate maximum element size, big enough to fit any collision algorithm in the memory pool
		
		int maxSize0 = sizeof(btFluidHfRigidCollisionAlgorithm);
		int maxSize1 = sizeof(btFluidHfBuoyantShapeCollisionAlgorithm);
		int maxSize2 = 0;

		int	collisionAlgorithmMaxElementSize = btMax(maxSize0,maxSize1);
		collisionAlgorithmMaxElementSize = btMax(collisionAlgorithmMaxElementSize,maxSize2);
		
		if (collisionAlgorithmMaxElementSize > curElemSize)
		{
			m_collisionAlgorithmPool->~btPoolAllocator();
			btAlignedFree(m_collisionAlgorithmPool);
			void* mem = btAlignedAlloc(sizeof(btPoolAllocator),16);
			m_collisionAlgorithmPool = new(mem) btPoolAllocator(collisionAlgorithmMaxElementSize,constructionInfo.m_defaultMaxCollisionAlgorithmPoolSize);
		}
	}
}

btFluidHfRigidCollisionConfiguration::~btFluidHfRigidCollisionConfiguration()
{
	m_fluidHfRigidConvexCreateFunc->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_fluidHfRigidConvexCreateFunc);
	
	m_fluidHfRigidConvexCreateFuncSwapped->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_fluidHfRigidConvexCreateFuncSwapped);
	
	m_fluidHfBuoyantShapeCollisionCreateFunc->~btCollisionAlgorithmCreateFunc();
	btAlignedFree(m_fluidHfBuoyantShapeCollisionCreateFunc);
}

btCollisionAlgorithmCreateFunc* btFluidHfRigidCollisionConfiguration::getCollisionAlgorithmCreateFunc(int proxyType0,int proxyType1)
{
	if(proxyType0 == HFFLUID_SHAPE_PROXYTYPE && proxyType1 == HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE)
	{
		return m_fluidHfRigidConvexCreateFunc;
	}

	if(proxyType0 == HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE && proxyType1 == HFFLUID_SHAPE_PROXYTYPE)
	{
		return m_fluidHfRigidConvexCreateFuncSwapped;
	}
	
#ifdef EXTEND_BT_FLUID_HF_BUOYANT_SHAPE_COLLISION 	//#define in btFluidHfBuoyantShapeCollisionAlgorithm.h

	bool collideProxyType1 = ( btBroadphaseProxy::isConvex(proxyType1)
								|| btBroadphaseProxy::isConcave(proxyType1)
								|| btBroadphaseProxy::isCompound(proxyType1) );
	if(proxyType0 == HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE && collideProxyType1)
	{
		return m_fluidHfBuoyantShapeCollisionCreateFunc;
	}	
	
	bool collideProxyType0 = ( btBroadphaseProxy::isConvex(proxyType0)
								|| btBroadphaseProxy::isConcave(proxyType0)
								|| btBroadphaseProxy::isCompound(proxyType0) );
	if(collideProxyType0 && proxyType1 == HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE)
	{
		return m_fluidHfBuoyantShapeCollisionCreateFunc;
	}
	
	if(proxyType0 == HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE && proxyType1 == HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE)
	{
		return m_fluidHfBuoyantShapeCollisionCreateFunc;
	}
#else
	if(proxyType0 == HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE && proxyType1 == HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE)
	{
		return m_fluidHfBuoyantShapeCollisionCreateFunc;
	}
#endif

	///fallback to the regular rigid collision shape
	return btDefaultCollisionConfiguration::getCollisionAlgorithmCreateFunc(proxyType0,proxyType1);
}

