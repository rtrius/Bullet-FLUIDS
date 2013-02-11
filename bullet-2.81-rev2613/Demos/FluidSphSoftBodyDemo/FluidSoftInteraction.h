#ifndef FLUID_SOFT_INTERACTION_H
#define FLUID_SOFT_INTERACTION_H

#include "LinearMath/btQuickProf.h"
#include "LinearMath/btAabbUtil2.h"
#include "BulletCollision/CollisionShapes/btSphereShape.h"
#include "BulletSoftBody/btSoftBodyInternals.h"
#include "BulletFluids/Sph/btFluidSph.h"

static void resolveFluidSoftContactImpulse(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, int sphParticleIndex,
									btSoftBody::Cluster* softCluster, const btVector3& normalOnSoftBody, 
									const btVector3&pointOnSoftBody, btScalar distance)
{
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	
	if( distance < btScalar(0.0) )
	{
		int i = sphParticleIndex;
		btFluidParticles& particles = fluid->internalGetParticles();
		
		btSoftBody::Body clusterBody(softCluster);
		btVector3 softLocalHitPoint = pointOnSoftBody - clusterBody.xform().getOrigin();
		
		btVector3 fluidVelocity = fluid->getVelocity(i);
		btVector3 softVelocity = clusterBody.velocity(softLocalHitPoint);
		softVelocity *= FG.m_simulationScale;
	
		btVector3 relativeVelocity = fluidVelocity - softVelocity;
		btScalar penetratingMagnitude = relativeVelocity.dot(-normalOnSoftBody);
		if( penetratingMagnitude < btScalar(0.0) ) penetratingMagnitude = btScalar(0.0);
		
		btVector3 penetratingVelocity = -normalOnSoftBody * penetratingMagnitude;
		btVector3 tangentialVelocity = relativeVelocity - penetratingVelocity;
		
		penetratingVelocity *= btScalar(1.0) + FL.m_boundaryRestitution;
		
		btScalar positionError = (-distance) * (FG.m_simulationScale/FG.m_timeStep) * FL.m_boundaryErp;
		btVector3 particleImpulse = -(penetratingVelocity + (-normalOnSoftBody*positionError) + tangentialVelocity*FL.m_boundaryFriction);
		
		const bool APPLY_IMPULSE_TO_SOFT_BODY = false;
		if(APPLY_IMPULSE_TO_SOFT_BODY)
		{
			btScalar inertiaParticle = btScalar(1.0) / FL.m_particleMass;
			
			btVector3 relPosCrossNormal = softLocalHitPoint.cross(normalOnSoftBody);
			btScalar inertiaSoft = clusterBody.invMass() + ( relPosCrossNormal * clusterBody.invWorldInertia() ).dot(relPosCrossNormal);
			
			particleImpulse *= btScalar(1.0) / (inertiaParticle + inertiaSoft);
			
			btVector3 worldScaleImpulse = -particleImpulse / FG.m_simulationScale;
			
			//	Apply the impulse to all nodes(vertices) in the soft body cluster
			//	this is incorrect, but sufficient for demonstration purposes
			btVector3 perNodeImpulse = worldScaleImpulse / static_cast<btScalar>( softCluster->m_nodes.size() );
			perNodeImpulse /= FG.m_timeStep;		//Impulse is accumulated as force
			
			for(int i = 0; i < softCluster->m_nodes.size(); ++i)
			{
				btSoftBody::Node* node = softCluster->m_nodes[i];
				node->m_f += perNodeImpulse;
			}
			
			
			particleImpulse *= inertiaParticle;
		}
		
		
		btVector3& vel = particles.m_vel[i];
		btVector3& vel_eval = particles.m_vel_eval[i];
		
		//Leapfrog integration
		btVector3 velNext = vel + particleImpulse;
		vel_eval = (vel + velNext) * btScalar(0.5);
		vel = velNext;
	}
}

///Preliminary soft body - SPH fluid interaction demo; this class is not supported.
struct ParticleSoftBodyCollisionCallback : public btDbvt::ICollide
{
	//All members must be set
	const btFluidSphParametersGlobal* m_globalParameters;
	
	btFluidSph* m_fluidSph;
	btSoftBody* m_softBody;
	
	int m_sphParticleIndex;
	btCollisionObject* m_particleObject;
	
	
	virtual void Process(const btDbvtNode* leaf)
	{
		BT_PROFILE("FluidSoft Process softBodyClusterLeaf");
		
		const btFluidSphParametersLocal& FL = m_fluidSph->getLocalParameters();
		
		btSoftBody::Cluster* softCluster = static_cast<btSoftBody::Cluster*>(leaf->data);
		btSoftClusterCollisionShape clusterShape(softCluster);
		
		btTransform& particleTransform = m_particleObject->getWorldTransform();
		btConvexShape* particleShape = static_cast<btConvexShape*>( m_particleObject->getCollisionShape() );
		
		btGjkEpaSolver2::sResults contactResult;            
		if( btGjkEpaSolver2::SignedDistance(&clusterShape, btTransform::getIdentity(), 
											particleShape, particleTransform, btVector3(1,0,0), contactResult) )
		{
			btVector3 normalOnCluster = -contactResult.normal;
			btVector3 pointOnCluster = contactResult.witnesses[0];
			btScalar distance = contactResult.distance - FL.m_particleRadius;
			
			resolveFluidSoftContactImpulse(*m_globalParameters, m_fluidSph, m_sphParticleIndex, 
											softCluster, normalOnCluster, pointOnCluster, distance);
		}
		
	}
};

///Preliminary soft body - SPH fluid interaction demo; this class is not supported.
struct FluidSphSoftBodyCollisionCallback : public btFluidSortingGrid::AabbCallback
{
	//All members must be set
	const btFluidSphParametersGlobal* m_globalParameters;
	
	btFluidSph* m_fluidSph;
	btSoftBody* m_softBody;
	
	btCollisionObject* m_particleObject;
	
	FluidSphSoftBodyCollisionCallback() {}
	
	
	virtual bool processParticles(const btFluidGridIterator FI, const btVector3& aabbMin, const btVector3& aabbMax)
	{
		BT_PROFILE("FluidSoft processParticles()");
		btTransform& particleTransform = m_particleObject->getWorldTransform();
		
		//Collide each SPH particle as a rigid body sphere against the soft body's clusters 
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
		{
			particleTransform.setOrigin( m_fluidSph->getPosition(n) );
		
			btVector3 particleMin, particleMax;
			m_particleObject->getCollisionShape()->getAabb(particleTransform, particleMin, particleMax);
			
			ParticleSoftBodyCollisionCallback dbvtCallback;
			dbvtCallback.m_globalParameters = m_globalParameters;
			dbvtCallback.m_fluidSph = m_fluidSph;
			dbvtCallback.m_softBody = m_softBody;
			dbvtCallback.m_sphParticleIndex = n;
			dbvtCallback.m_particleObject = m_particleObject;
		
			const btDbvt& softBodyClusterDbvt = m_softBody->m_cdbvt;
			btDbvtVolume particleAabb = btDbvtVolume::FromMM(particleMin, particleMax);
			
			softBodyClusterDbvt.collideTV(softBodyClusterDbvt.m_root, particleAabb, dbvtCallback);
		}
		
		return true;
	}

};

///Preliminary soft body - SPH fluid interaction demo; this class is not supported.
struct FluidSoftInteractor
{
	static void collideFluidSphWithSoftBody(const btFluidSphParametersGlobal& FG, 
											btFluidSph* fluidSph, btSoftBody* softBody, const btVector3& softAabbMin, const btVector3& softAabbMax)
	{
		BT_PROFILE("collideFluidSphWithSoftBody()");
	
		const btFluidSphParametersLocal& FL = fluidSph->getLocalParameters();
		const btFluidSortingGrid& grid = fluidSph->getGrid();
		
		btSphereShape particleShape(FL.m_particleRadius);
		
		btCollisionObject particleObject;
		particleObject.setCollisionShape(&particleShape);
		particleObject.getWorldTransform().setIdentity();
		
		//btFluidSortingGrid uses only the particle centers(without radius)
		//add the radius to the soft body AABB to avoid missing collisions
		const btVector3 fluidRadius(FL.m_particleRadius, FL.m_particleRadius, FL.m_particleRadius);
		btVector3 expandedSoftMin = softAabbMin - fluidRadius;
		btVector3 expandedSoftMax = softAabbMax + fluidRadius;
		
		FluidSphSoftBodyCollisionCallback fluidSoftCollider;
		fluidSoftCollider.m_globalParameters = &FG;
		fluidSoftCollider.m_fluidSph = fluidSph;
		fluidSoftCollider.m_softBody = softBody;
		fluidSoftCollider.m_particleObject = &particleObject;
		
		//Call FluidSphSoftBodyCollisionCallback::processParticles() for
		//each SPH fluid grid cell intersecting with the soft body's AABB
		grid.forEachGridCell(expandedSoftMin, expandedSoftMax, fluidSoftCollider);
	}
	
	static void performCollisionDetectionAndResponse(const btFluidSphParametersGlobal& FG, btAlignedObjectArray<btFluidSph*>& fluids,
													btAlignedObjectArray<btSoftBody*>& softBodies, btScalar timeStep)
	{
		BT_PROFILE("btFluidSph - btSoftBody interaction");
	
		for(int i = 0; i < fluids.size(); ++i)
		{
			btFluidSph* fluidSph = fluids[i];
			if( !fluidSph->numParticles() ) continue;
			
			btVector3 fluidAabbMin, fluidAabbMax;
			fluidSph->getAabb(fluidAabbMin, fluidAabbMax);
			
			for(int n = 0; n < softBodies.size(); ++n)
			{
				btSoftBody* softBody = softBodies[n];
				btVector3 softAabbMin, softAabbMax;
				softBody->getAabb(softAabbMin, softAabbMax);
			
				if( TestAabbAgainstAabb2(fluidAabbMin, fluidAabbMax, softAabbMin, softAabbMax) )
					collideFluidSphWithSoftBody(FG, fluidSph, softBody, softAabbMin, softAabbMax);
			}
		}
	}
};


#endif