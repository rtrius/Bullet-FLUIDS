#ifndef FLUID_SPH_HF_INTERACTION_H
#define FLUID_SPH_HF_INTERACTION_H

#include "LinearMath/btQuickProf.h"
#include "LinearMath/btAlignedObjectArray.h"

#include "BulletFluids/Sph/btFluidSph.h"
#include "BulletFluids/Hf/btFluidHf.h"

class FluidSphHfInteractor
{
	
	struct SphHfContact
	{
		int m_particleIndex;
		int m_hfColumnIndex;
		btScalar m_distance;	//World scale distance from the particle to the heightfield fluid surface
		
		SphHfContact() {}
		SphHfContact(int particleIndex, int hfColumnIndex, btScalar distance) 
		: m_particleIndex(particleIndex), m_hfColumnIndex(hfColumnIndex), m_distance(distance) {}
	};
	
	struct DescendingDistanceSortPredicate
	{
		inline bool operator() (const SphHfContact& a, const SphHfContact& b) const 
		{
			return (a.m_distance > b.m_distance);
		}
	};


	btAlignedObjectArray<SphHfContact> m_contacts;
	
public:
	void interact(const btFluidSphParametersGlobal& FG, btFluidSph* fluidSph, btFluidHf* fluidHf)
	{
		detectCollisions(FG, fluidSph, fluidHf);
		//resolveCollisionsAbsorb(FG, fluidSph, fluidHf);
		resolveCollisionsPartialAbsorb(FG, fluidSph, fluidHf);
	}

	void detectCollisions(const btFluidSphParametersGlobal& FG, btFluidSph* fluidSph, btFluidHf* fluidHf)
	{
		BT_PROFILE("FluidSphHfInteractor::detectCollisions()");
	
		m_contacts.resize(0);
	
		const btFluidSphParametersLocal& FL = fluidSph->getLocalParameters();
		
		for(int n = 0; n < fluidSph->numParticles(); ++n)
		{
			const btVector3& transformedParticlePos = fluidHf->getWorldTransform().invXform( fluidSph->getPosition(n) );
			
			int index_x = fluidHf->getCellPosX( transformedParticlePos.x() );
			int index_z = fluidHf->getCellPosZ( transformedParticlePos.z() );
			if( 1 <= index_x && index_x < fluidHf->getNumNodesX() - 1
			 && 1 <= index_z && index_z < fluidHf->getNumNodesZ() - 1 )	//AABB test
			{
				const btScalar MARGIN(0.05);
				
				btScalar particleHeight = transformedParticlePos.y() - FL.m_particleRadius - MARGIN;
				int columnIndex = index_x + index_z*fluidHf->getNumNodesX();
				if( particleHeight < fluidHf->getCombinedHeight(columnIndex) )
				{
					btScalar distance = particleHeight - fluidHf->getCombinedHeight(columnIndex);
				
					if( distance < btScalar(0.0) ) m_contacts.push_back( SphHfContact(n, columnIndex, distance) );
				}
			}
		}
	}
	
	///btHfFluid absorbs all particles in contact
	void resolveCollisionsAbsorb(const btFluidSphParametersGlobal& FG, btFluidSph* fluidSph, btFluidHf* fluidHf)
	{
		BT_PROFILE("FluidSphHfInteractor::resolveCollisionsAbsorb()");
		
		const btScalar PARTICLE_HEIGHT_CONTRIBUTION(3.0);
			
		for(int i = 0; i < m_contacts.size(); ++i)
		{
			fluidSph->markParticleForRemoval(m_contacts[i].m_particleIndex);
			fluidHf->addFluidHeight(m_contacts[i].m_hfColumnIndex, PARTICLE_HEIGHT_CONTRIBUTION);
		}
	}
	
	///btHfFluid absorbs deepest penetrating particles when the 
	///number of particles exceeds a threshold, and pushes other particles up
	void resolveCollisionsPartialAbsorb(const btFluidSphParametersGlobal& FG, btFluidSph* fluidSph, btFluidHf* fluidHf)
	{
		BT_PROFILE("FluidSphHfInteractor::resolveCollisionsPartialAbsorb()");
		
		const btFluidSphParametersLocal& FL = fluidSph->getLocalParameters();
		
		//Absorb deepest penetrating particles
		const int MAX_PARTICLES = 4096;
		const btScalar PARTICLE_HEIGHT_CONTRIBUTION(1.0);
		if( fluidSph->numParticles() > MAX_PARTICLES )
		{
			//Sort contacts, from highest to lowest distance(least to most penetrating)
			m_contacts.quickSort( FluidSphHfInteractor::DescendingDistanceSortPredicate() );
			
			//Absorb most penetrating particles, and discard those contacts
			int numPotentiallyAbsorbed = fluidSph->numParticles() - MAX_PARTICLES;
			int numAbsorbed = ( m_contacts.size() > numPotentiallyAbsorbed ) ? numPotentiallyAbsorbed : m_contacts.size();
			int numRemainingAfterAbsorb = m_contacts.size() - numAbsorbed;
			
			for(int i = numRemainingAfterAbsorb; i < m_contacts.size(); ++i)	//	verify
			{
				fluidSph->markParticleForRemoval(m_contacts[i].m_particleIndex);
				fluidHf->addFluidHeight(m_contacts[i].m_hfColumnIndex, PARTICLE_HEIGHT_CONTRIBUTION);
			}
			
			//Discard absorbed contacts
			m_contacts.resize(numRemainingAfterAbsorb);
		}
		
		//Push other particles out of the heightfield fluid
		const btScalar PARTICLE_HEIGHT_DISPLACEMENT(0.5);
		const btScalar STIFFNESS(20000.0);
		const btScalar DAMPING(256.0);
		const btVector3 NORMAL(0, 1, 0);
		for(int i = 0; i < m_contacts.size(); ++i)
		{
			int particleIndex = m_contacts[i].m_particleIndex;
			//int hfColumnIndex = m_contacts[i].m_hfColumnIndex;
		
			const btVector3& velocity = fluidSph->getVelocity(particleIndex);
			btScalar penetration = (-m_contacts[i].m_distance) * FG.m_simulationScale;
			btScalar magnitude = STIFFNESS*penetration - DAMPING*velocity.dot(NORMAL);	//Simulation scaled acceleration
		
			//fluidHf->addDisplaced(index_x, index_z, PARTICLE_HEIGHT_DISPLACEMENT);
			fluidSph->applyForce(particleIndex, NORMAL*magnitude*FL.m_particleMass);
		}
	}
};

#endif