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
#include "btFluidSolverMultiphase.h"

#include "LinearMath/btAabbUtil2.h"		//TestAabbAgainstAabb2()


void btFluidSolverMultiphase::updateGridAndCalculateSphForces(const btFluidParametersGlobal& FG, btFluidSph** fluids, int numFluids)
{
	BT_PROFILE("btFluidSolverMultiphase::updateGridAndCalculateSphForces()");
	
	m_sphData.resize(numFluids);
	for(int i = 0; i < numFluids; ++i) m_sphData[i].resize( fluids[i]->numParticles() );
			
	//
	for(int i = 0; i < numFluids; ++i) fluids[i]->insertParticlesIntoGrid();
	
	//Determine intersecting btFluidSph AABBs
	btAlignedObjectArray< btAlignedObjectArray<btFluidSph*> > interactingFluids;
	btAlignedObjectArray< btAlignedObjectArray<btFluidSolverSph::SphParticles*> > interactingSphData;
	interactingFluids.resize(numFluids);
	interactingSphData.resize(numFluids);
	
	for(int i = 0; i < numFluids; ++i)
	{
		interactingFluids[i].resize(0);
		interactingSphData[i].resize(0);
	
		btVector3 min_i, max_i;
		fluids[i]->getAabb(min_i, max_i);
	
		for(int n = 0; n < numFluids; ++n)
		{
			if(i == n) continue;
			
			btVector3 min_n, max_n;
			fluids[n]->getAabb(min_n, max_n);
			
			if( TestAabbAgainstAabb2(min_i, max_i, min_n, max_n) )
			{
				interactingFluids[i].push_back( fluids[n] );
				interactingSphData[i].push_back( &m_sphData[n] );
			}
		}
	}
	
	//
	for(int i = 0; i < numFluids; ++i) 
		sphComputePressureMultiphase( FG, fluids[i], m_sphData[i], interactingFluids[i], interactingSphData[i] );
		
	for(int i = 0; i < numFluids; ++i) 
		sphComputeForceMultiphase( FG, fluids[i], m_sphData[i], interactingFluids[i], interactingSphData[i] );
		
	for(int i = 0; i < numFluids; ++i)
	{
		applySphForce(FG, fluids[i], m_sphData[i].m_sphForce);
	}
}

void btFluidSolverMultiphase::sphComputePressureMultiphase(const btFluidParametersGlobal& FG, btFluidSph* fluid, 
															btFluidSolverSph::SphParticles& sphData,
															btAlignedObjectArray<btFluidSph*>& interactingFluids, 
															btAlignedObjectArray<btFluidSolverSph::SphParticles*>& interactingSphData)
{
	BT_PROFILE("sphComputePressureMultiphase()");
	
	const btFluidParametersLocal& FL = fluid->getLocalParameters();
	const btFluidSortingGrid& grid = fluid->getGrid();
	btFluidParticles& particles = fluid->internalGetParticles();

	for(int i = 0; i < fluid->numParticles(); ++i)
	{
//#define DENSITY_CONTRAST
#ifndef DENSITY_CONTRAST
		btScalar sum = 0.0f;
#else
		btScalar sum = FG.m_sphRadiusSquared*FG.m_sphRadiusSquared*FG.m_sphRadiusSquared;	//Self contributed density
#endif
		sphData.m_neighborTable[i].clear();

		btFluidSortingGrid::FoundCells foundCells;
		grid.findCells(particles.m_pos[i], foundCells);
		
		for(int cell = 0; cell < btFluidSortingGrid::NUM_FOUND_CELLS; cell++) 
		{
			btFluidGridIterator& FI = foundCells.m_iterators[cell];
			
			for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
			{
				if(i == n) continue;
				
				btVector3 distance = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
				btScalar distanceSquared = distance.length2();
				
				if(FG.m_sphRadiusSquared > distanceSquared) 
				{
					btScalar c = FG.m_sphRadiusSquared - distanceSquared;
					sum += c * c * c;
					
					if( !sphData.m_neighborTable[i].isFilled() ) sphData.m_neighborTable[i].addNeighbor( n, btSqrt(distanceSquared) );
					else break;
				}
			}
		}
#ifndef DENSITY_CONTRAST		
		btScalar density = sum * FL.m_particleMass * FG.m_poly6KernCoeff;
#else
		btScalar particleDensity = sum * FG.m_poly6KernCoeff;
#endif
		
		//EXTERNAL_FLUID_INTERACTION
		for(int j = 0; j < interactingFluids.size(); ++j)
		{
			btFluidSph* externalFluid = interactingFluids[j];
		
			const btFluidParametersLocal& externalFL = externalFluid->getLocalParameters();
			const btFluidSortingGrid& externalGrid = externalFluid->getGrid();
			btFluidParticles& externalParticles = externalFluid->internalGetParticles();
			
			btScalar externalSum = 0.0;	
			btFluidSortingGrid::FoundCells externalFoundCells;
			externalGrid.findCells(particles.m_pos[i], externalFoundCells);
			for(int cell = 0; cell < btFluidSortingGrid::NUM_FOUND_CELLS; cell++) 
			{
				btFluidGridIterator& FI = externalFoundCells.m_iterators[cell];
				for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
				{
					btVector3 distance = (particles.m_pos[i] - externalParticles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
					btScalar distanceSquared = distance.length2();
				
					if(FG.m_sphRadiusSquared > distanceSquared) 
					{
						btScalar c = FG.m_sphRadiusSquared - distanceSquared;
						externalSum += c * c * c;
					}
				}
			}
#ifndef DENSITY_CONTRAST				
			density += externalSum * externalFL.m_particleMass * FG.m_poly6KernCoeff;
#else
			particleDensity += externalSum * FG.m_poly6KernCoeff;
#endif
		}
		//EXTERNAL_FLUID_INTERACTION

#ifndef DENSITY_CONTRAST
		sphData.m_pressure[i] = (density - FL.m_restDensity) * FL.m_stiffness;
		sphData.m_invDensity[i] = 1.0f / density;
#else
		btScalar density = particleDensity * FL.m_particleMass;
		sphData.m_pressure[i] = (density - FL.m_restDensity) * FL.m_stiffness;
		sphData.m_invDensity[i] = 1.0f / density;
#endif
	}
	
}

void computeForceNeighborTable_Multiphase(const btFluidParametersGlobal& FG, const btScalar vterm, int particleIndex, 
											btFluidParticles& particles, btFluidSolverSph::SphParticles& sphData)
{
	int i = particleIndex;

	btVector3 force(0, 0, 0);
	for(int j = 0; j < sphData.m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = sphData.m_neighborTable[i].getNeighborIndex(j);
		
		btVector3 distance = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
		
		btScalar c = FG.m_sphSmoothRadius - sphData.m_neighborTable[i].getDistance(j);
		btScalar pterm = -0.5f * c * FG.m_spikyKernGradCoeff 
					 * ( sphData.m_pressure[i] + sphData.m_pressure[n]) / sphData.m_neighborTable[i].getDistance(j);
		btScalar dterm = c * sphData.m_invDensity[i] * sphData.m_invDensity[n];

		btVector3 forceAdded( (pterm * distance.x() + vterm * (particles.m_vel_eval[n].x() - particles.m_vel_eval[i].x())) * dterm,
							  (pterm * distance.y() + vterm * (particles.m_vel_eval[n].y() - particles.m_vel_eval[i].y())) * dterm,
							  (pterm * distance.z() + vterm * (particles.m_vel_eval[n].z() - particles.m_vel_eval[i].z())) * dterm );
		force += forceAdded;
	}
	
	sphData.m_sphForce[i] = force;
}

void btFluidSolverMultiphase::sphComputeForceMultiphase(const btFluidParametersGlobal& FG, btFluidSph* fluid, 
															btFluidSolverSph::SphParticles& sphData,
															btAlignedObjectArray<btFluidSph*>& interactingFluids, 
															btAlignedObjectArray<btFluidSolverSph::SphParticles*>& interactingSphData)
{
	BT_PROFILE("sphComputeForceMultiphase()");
	
	const btFluidParametersLocal& FL = fluid->getLocalParameters();
	btFluidParticles& particles = fluid->internalGetParticles();
	btScalar vterm = FG.m_viscosityKernLapCoeff * FL.m_viscosity;
	
	for(int i = 0; i < particles.size(); ++i)
	{
		computeForceNeighborTable_Multiphase(FG, vterm, i, particles, sphData);
		sphData.m_sphForce[i] *= FL.m_particleMass;
	}
	
	//EXTERNAL_FLUID_INTERACTION
	btScalar radius = FG.m_sphSmoothRadius / FG.m_simulationScale;
	
	for(int j = 0; j < interactingFluids.size(); ++j)
	{
		btFluidSph* externalFluid = interactingFluids[j];
		btFluidSolverSph::SphParticles& externalSphData = *interactingSphData[j];
		
		const btFluidParametersLocal& externalFL = externalFluid->getLocalParameters();
		const btFluidSortingGrid& externalGrid = externalFluid->getGrid();
		btFluidParticles& externalParticles = externalFluid->internalGetParticles();
		
		btScalar averagedViscosity = (FL.m_viscosity + externalFL.m_viscosity) * 0.5f;
		btScalar vterm2 = FG.m_viscosityKernLapCoeff * averagedViscosity;
		
		for(int i = 0; i < particles.size(); ++i)
		{
			btVector3 externalForce(0, 0, 0);
			btFluidSortingGrid::FoundCells externalFoundCells;
			externalGrid.findCells(particles.m_pos[i], externalFoundCells);
			for(int cell = 0; cell < btFluidSortingGrid::NUM_FOUND_CELLS; cell++) 
			{
				btFluidGridIterator& FI = externalFoundCells.m_iterators[cell];
				for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
				{
					btVector3 distance = (particles.m_pos[i] - externalParticles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
					btScalar distanceSquared = distance.length2();
					
					if(FG.m_sphRadiusSquared > distanceSquared) 
					{
						const bool USE_FLUIDS_V2_EQUATIONS = true;
						if(USE_FLUIDS_V2_EQUATIONS)
						{
							btScalar r = btSqrt(distanceSquared);
							
							btScalar c = FG.m_sphSmoothRadius - r;
							btScalar pterm = -0.5f * c * FG.m_spikyKernGradCoeff * ( sphData.m_pressure[i] + externalSphData.m_pressure[n]) / r;
							btScalar dterm = c * sphData.m_invDensity[i] * externalSphData.m_invDensity[n];
							
							btVector3 forceAdded( (pterm * distance.x() + vterm2 * (externalParticles.m_vel_eval[n].x() - particles.m_vel_eval[i].x())) * dterm,
												  (pterm * distance.y() + vterm2 * (externalParticles.m_vel_eval[n].y() - particles.m_vel_eval[i].y())) * dterm,
												  (pterm * distance.z() + vterm2 * (externalParticles.m_vel_eval[n].z() - particles.m_vel_eval[i].z())) * dterm );
							externalForce += forceAdded;
						}
						else
						{
							btScalar r = btSqrt(distanceSquared);
							btScalar c = FG.m_sphSmoothRadius - r;
						
								//Non-kernel
							btScalar pressureScalar = -externalFL.m_particleMass;
							pressureScalar *= btScalar(0.5) * (sphData.m_pressure[i] + externalSphData.m_pressure[n]) * externalSphData.m_invDensity[n];
								//Kernel
							pressureScalar *= FG.m_spikyKernGradCoeff;
							pressureScalar *= (c * c) / r;

							btVector3 pressureForce = distance * pressureScalar;
							
								//Non-kernel
							btScalar viscosityScalar = averagedViscosity * externalFL.m_particleMass;
							viscosityScalar *= externalSphData.m_invDensity[n];
								//Kernel
							viscosityScalar *= FG.m_viscosityKernLapCoeff;
							viscosityScalar *= c;
							
							btVector3 viscosityForce = (externalParticles.m_vel_eval[n] - particles.m_vel_eval[i]) * viscosityScalar;
							
							
							btVector3 forceAdded2 = pressureForce + viscosityForce;
							externalForce += forceAdded2;
						}
					}
				}
			}
			
			sphData.m_sphForce[i] += externalForce *= externalFL.m_particleMass;
		}
	}
	//EXTERNAL_FLUID_INTERACTION
}
