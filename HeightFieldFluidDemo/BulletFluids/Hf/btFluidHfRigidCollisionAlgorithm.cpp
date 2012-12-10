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

#include "btFluidHfRigidCollisionAlgorithm.h"

#include "BulletCollision/CollisionDispatch/btCollisionDispatcher.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletCollision/CollisionShapes/btBoxShape.h"
#include "BulletCollision/CollisionDispatch/btCollisionObject.h"
#include "BulletDynamics/Dynamics/btRigidBody.h"

#include "btFluidHfBuoyantConvexShape.h"
#include "btFluidHf.h"


btFluidHfRigidCollisionAlgorithm::btFluidHfRigidCollisionAlgorithm(const btCollisionAlgorithmConstructionInfo& ci, 
																	const btCollisionObjectWrapper* body0Wrap,
																	const btCollisionObjectWrapper* body1Wrap, bool isSwapped)
: btCollisionAlgorithm(ci), m_isSwapped(isSwapped), 
	m_convexTrianglecallback(ci.m_dispatcher1, body0Wrap, body1Wrap, !isSwapped) // we flip the isSwapped because we are hf fluid vs. convex and callback expects convex vs. concave
{
}

void btFluidHfRigidCollisionAlgorithm::processGround(const btCollisionObjectWrapper* hfFluidWrap, const btCollisionObjectWrapper* rigidWrap,
													const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	// to perform the convex shape vs. ground terrain:
	// we pull the convex shape out of the btFluidHfBuoyantConvexShape and replace it temporarily
	const btFluidHfBuoyantConvexShape* tmpRigidShape = static_cast<const btFluidHfBuoyantConvexShape*>( rigidWrap->getCollisionShape() );
	const btConvexShape* convexShape = tmpRigidShape->getConvexShape();
	
	btCollisionObjectWrapper tempRigidWrap( rigidWrap, convexShape, rigidWrap->getCollisionObject(), rigidWrap->getWorldTransform() );
	
	btScalar triangleMargin = m_rigidCollisionObject->getCollisionShape()->getMargin();
	resultOut->setPersistentManifold(m_convexTrianglecallback.m_manifoldPtr);
	
	m_convexTrianglecallback.setTimeStepAndCounters(triangleMargin, dispatchInfo, &tempRigidWrap, hfFluidWrap, resultOut);
	m_hfFluid->forEachGroundTriangle(&m_convexTrianglecallback, m_convexTrianglecallback.getAabbMin(), m_convexTrianglecallback.getAabbMax());
	
	resultOut->refreshContactPoints();
	m_convexTrianglecallback.clearWrapperData();
}


inline btScalar linearInterpolate(btScalar a, btScalar b, btScalar x) { return (btScalar(1.0) - x) * a + x * b; }
inline btScalar bilinearInterpolate(btScalar x1, btScalar x2, btScalar y1, btScalar y2, btScalar x, btScalar y) 
{
	btScalar a = linearInterpolate(x1, x2, x);
	btScalar b = linearInterpolate(y1, y2, x);

	return linearInterpolate(a, b, y);
}

class btIDebugDraw;
class btFluidHfColumnRigidBodyCallback : public btFluidHf::btFluidHfColumnCallback
{
protected:
	btRigidBody* m_rigidBody;
	btFluidHfBuoyantConvexShape* m_buoyantShape;
	btIDebugDraw* m_debugDraw;
	
	int m_numVoxels;
	btAlignedObjectArray<btVector3> m_transformedVoxelPositions;
	btAlignedObjectArray<bool> m_voxelSubmerged;
	
	btScalar m_submergedVolume;
	btScalar m_density;
	btScalar m_floatyness;
	btScalar m_timeStep;
	
public:
	btFluidHfColumnRigidBodyCallback(btRigidBody* rigidBody, btIDebugDraw* debugDraw, btScalar density, btScalar floatyness, btScalar timeStep)
	{
		m_rigidBody = rigidBody;
		m_buoyantShape = (btFluidHfBuoyantConvexShape*)rigidBody->getCollisionShape();
		m_debugDraw = debugDraw;
		
		m_numVoxels = m_buoyantShape->getNumVoxels();
		{
			m_transformedVoxelPositions.resize(m_numVoxels);
			m_voxelSubmerged.resize(m_numVoxels);
			
			for(int i = 0; i < m_numVoxels; i++)
			{
				m_transformedVoxelPositions[i] = m_rigidBody->getWorldTransform() * m_buoyantShape->getVoxelPosition(i);
				m_voxelSubmerged[i] = false;
			}
		}
		
		m_submergedVolume = btScalar(0.0);
		m_density = density;
		m_floatyness = floatyness;
		m_timeStep = timeStep;
	}
	~btFluidHfColumnRigidBodyCallback() {}
	
	static bool sphereVsAABB (const btVector3& aabbMin, const btVector3& aabbMax, const btVector3& sphereCenter, const btScalar sphereRadius)
	{
		btScalar totalDistance = 0;

		// Accumulate the distance of the sphere's center on each axis
		for(int i = 0; i < 3; ++i) 
		{
			// If the sphere's center is outside the aabb, we've got distance on this axis
			if(sphereCenter[i] < aabbMin[i])
			{
				btScalar borderDistance = aabbMin[i] - sphereCenter[i];
				totalDistance += borderDistance * borderDistance;
			} 
			else if(sphereCenter[i] > aabbMax[i])
			{
				btScalar borderDistance = sphereCenter[i] - aabbMax[i];
				totalDistance += borderDistance * borderDistance;
			}
			// Otherwise the sphere's center is within the box on this axis, so the
			// distance will be 0 and we do not need to accumulate anything at all
		}

		// If the distance to the box is lower than the sphere's radius, both are overlapping
		return (totalDistance <= (sphereRadius * sphereRadius));
	}
	bool processColumn(btFluidHf* fluid, int w, int l)
	{
		btScalar columnVolume = btScalar(0.0);

		const btScalar voxelDiameter = m_buoyantShape->getVoxelRadius() * btScalar(2.0);
		for(int i = 0; i < m_numVoxels; i++)
		{
			if(m_voxelSubmerged[i]) continue;
			
			if(0)	//Use more gradual interaction
			{
				const btVector3& fluidPosition = fluid->getWorldTransform().getOrigin();
				const btVector3& voxelPosition = m_transformedVoxelPositions[i];
			
				const btScalar cellWidth = fluid->getGridCellWidth();
				int x = static_cast<int>( (voxelPosition.x() - fluidPosition.x()) / cellWidth );
				int z = static_cast<int>( (voxelPosition.z() - fluidPosition.z()) / cellWidth );
				
				//if( 1 <= x && x < fluid->getNumNodesX()-1 && 1 <= z && z < fluid->getNumNodesZ()-1 )
				if( x == w && z == l )
				{
					btScalar voxelMinY = voxelPosition.y() - m_buoyantShape->getVoxelRadius();
					
					btScalar xContinuous = ( (voxelPosition.x() - fluidPosition.x()) / cellWidth ) - static_cast<btScalar>(x);
					btScalar zContinuous = ( (voxelPosition.z() - fluidPosition.z()) / cellWidth ) - static_cast<btScalar>(z);
					
					//btScalar height1 = fluid->getCombinedHeight( fluid->arrayIndex(x,z) );
					//btScalar height2 = fluid->getCombinedHeight( fluid->arrayIndex(x+1,z) );
					//btScalar height3 = fluid->getCombinedHeight( fluid->arrayIndex(x,z+1) );
					//btScalar height4 = fluid->getCombinedHeight( fluid->arrayIndex(x+1,z+1) );
					//btScalar fluidHeight = bilinearInterpolate(height1, height2, height3, height4, xContinuous, zContinuous) + fluidPosition.y();
					btScalar fluidHeight = fluid->getCombinedHeight( fluid->arrayIndex(x,z) ) + fluidPosition.y();
					
					btScalar submergedDepth = fluidHeight - voxelMinY;
					if(submergedDepth > SIMD_EPSILON)
					{
						m_voxelSubmerged[i] = true;
					
						btScalar submergedFraction = btMin( submergedDepth / voxelDiameter, btScalar(1.0) );
						btScalar submergedVolume = m_buoyantShape->getVolumePerVoxel() * submergedFraction;
						
						columnVolume += submergedVolume;

						const bool APPLY_BUOYANCY_IMPULSE = true;
						if(APPLY_BUOYANCY_IMPULSE)
						{
							btScalar massDisplacedWater = submergedVolume * m_density * m_floatyness;
							btScalar force = massDisplacedWater * -fluid->getParameters().m_gravity;
							btScalar impulseMag = force * m_timeStep;
							btVector3 impulseVec = btVector3( btScalar(0.0), btScalar(1.0), btScalar(0.0) ) * impulseMag;
							
							btVector3 application_point = voxelPosition;
							btVector3 relative_position = application_point - m_rigidBody->getCenterOfMassPosition();
							m_rigidBody->applyImpulse(impulseVec, relative_position);
						}
					}
				}
				
				continue;	//Ignore rest of loop
			}
			
			btVector3 columnAabbMin, columnAabbMax;
			fluid->getAabbForColumn(w, l, columnAabbMin, columnAabbMax);
			if( sphereVsAABB(columnAabbMin, columnAabbMax, m_transformedVoxelPositions[i], m_buoyantShape->getVoxelRadius()) )
			{
				m_voxelSubmerged[i] = true;
				btScalar voxelVolume = m_buoyantShape->getVolumePerVoxel();
				columnVolume += voxelVolume;

				const bool APPLY_BUOYANCY_IMPULSE = true;
				if(APPLY_BUOYANCY_IMPULSE)
				{
					btScalar massDisplacedWater = voxelVolume * m_density * m_floatyness;
					btScalar force = massDisplacedWater * -fluid->getParameters().m_gravity;
					btScalar impulseMag = force * m_timeStep;
					btVector3 impulseVec = btVector3( btScalar(0.0), btScalar(1.0), btScalar(0.0) ) * impulseMag;
					
					btVector3 application_point = m_transformedVoxelPositions[i];
					btVector3 relative_position = application_point - m_rigidBody->getCenterOfMassPosition();
					m_rigidBody->applyImpulse(impulseVec, relative_position);
				}
			}
		}

		if( columnVolume > btScalar(0.0) )
		{
			m_submergedVolume += columnVolume;

			const btScalar MIN_DISPLACE_HEIGHT(0.1);
			const bool DISPLACE_FLUID = true;
			if(DISPLACE_FLUID && fluid->getFluidHeight( fluid->arrayIndex(w,l) ) > MIN_DISPLACE_HEIGHT) fluid->addDisplaced(w, l, columnVolume);
			
			const bool APPLY_FLUID_VELOCITY_IMPULSE = true;
			if(APPLY_FLUID_VELOCITY_IMPULSE)
			{
				int index = fluid->arrayIndex (w,l);
				btVector3 velDelta = btVector3( fluid->getVelocityX(index), btScalar(0.0), fluid->getVelocityZ(index) );
				btVector3 impulse = velDelta * m_timeStep * fluid->getParameters().m_horizontalVelocityScale;
				
				m_rigidBody->applyCentralImpulse(impulse);
			}
		}

		return true;
	}
	
	btScalar getSubmergedVolume() const { return m_submergedVolume; }
};

void btFluidHfRigidCollisionAlgorithm::processCollision(const btCollisionObjectWrapper* body0Wrap, const btCollisionObjectWrapper* body1Wrap,
														const btDispatcherInfo& dispatchInfo, btManifoldResult* resultOut)
{
	const btCollisionObjectWrapper* hfFluidWrap = (m_isSwapped) ? body1Wrap : body0Wrap;
	const btCollisionObjectWrapper* rigidWrap = (m_isSwapped) ? body0Wrap : body1Wrap;	
	
	m_hfFluid = static_cast<btFluidHf*>( const_cast<btCollisionObject*>(hfFluidWrap->getCollisionObject()) );
	m_rigidCollisionObject = const_cast<btCollisionObject*>( rigidWrap->getCollisionObject() );
	
	//Collide rigid body with ground
	processGround(hfFluidWrap, rigidWrap, dispatchInfo, resultOut);
	
	//Collide rigid body with fluid
	btRigidBody* rb = btRigidBody::upcast(m_rigidCollisionObject);
	if( rb && rb->getInvMass() != btScalar(0.0) )
	{
		btFluidHfBuoyantConvexShape* buoyantShape = (btFluidHfBuoyantConvexShape*)m_rigidCollisionObject->getCollisionShape();
		btScalar volume = buoyantShape->getTotalVolume();
		
		btScalar mass = btScalar(1.0) / rb->getInvMass();
		btScalar density = mass / volume;
		
		btScalar floatyness = buoyantShape->getFloatyness();
		
		//Collide btRigidBody voxels with btHfFluid columns
		btFluidHfColumnRigidBodyCallback columnCallback(rb, dispatchInfo.m_debugDraw, density, floatyness, dispatchInfo.m_timeStep);
		m_hfFluid->forEachFluidColumn(&columnCallback, m_convexTrianglecallback.getAabbMin(), m_convexTrianglecallback.getAabbMax());
		
		//Apply fluid friction
		btScalar submerged_volume = columnCallback.getSubmergedVolume();
		if( submerged_volume > btScalar(0.0) )
		{
			btScalar submerged_percentage = submerged_volume / volume;
			
			const btScalar mu = btScalar(1.5);
			btScalar scaled_mu = mu * submerged_percentage * dispatchInfo.m_timeStep;
			rb->applyCentralImpulse( scaled_mu * -rb->getLinearVelocity() );
			rb->applyTorqueImpulse( scaled_mu * -rb->getAngularVelocity() );
		}
		
		rb->activate(true);
	}
}
