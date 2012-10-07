/* FluidSolverMultiphase.cpp

	ZLib license
	This software is provided 'as-is', without any express or implied
	warranty. In no event will the authors be held liable for any damages
	arising from the use of this software.
	
	Permission is granted to anyone to use this software for any purpose,
	including commercial applications, and to alter it and redistribute it
	freely, subject to the following restrictions:
	
	1. The origin of this software must not be misrepresented; you must not
	   claim that you wrote the original software. If you use this software
	   in a product, an acknowledgment in the product documentation would be
	   appreciated but is not required.
	2. Altered source versions must be plainly marked as such, and must not be
	   misrepresented as being the original software.
	3. This notice may not be removed or altered from any source distribution.
*/
#include "FluidSolverMultiphase.h"

#include "LinearMath/btAabbUtil2.h"		//TestAabbAgainstAabb2()

void FluidSolverMultiphase::stepSimulation(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids)
{
	BT_PROFILE("FluidSolverMultiphase::stepSimulation()");

	//
	for(int i = 0; i < fluids->size(); ++i) (*fluids)[i]->insertParticlesIntoGrid();
	
	//Determine intersecting FluidSph AABBs
	btAlignedObjectArray< btAlignedObjectArray<FluidSph*> > interactingFluids;
	interactingFluids.resize( fluids->size() );
	for(int i = 0; i < fluids->size(); ++i)
	{
		interactingFluids[i].resize(0);
	
		btVector3 min_i, max_i;
		(*fluids)[i]->getCurrentAabb(FG, &min_i, &max_i);
	
		for(int n = 0; n < fluids->size(); ++n)
		{
			if(i == n) continue;
			
			btVector3 min_n, max_n;
			(*fluids)[n]->getCurrentAabb(FG, &min_n, &max_n);
			
			if( TestAabbAgainstAabb2(min_i, max_i, min_n, max_n) )interactingFluids[i].push_back( (*fluids)[n] );
		}
	}
	
	//
	for(int i = 0; i < fluids->size(); ++i) 
		sphComputePressure( FG, (*fluids)[i], &interactingFluids[i] );
		
	for(int i = 0; i < fluids->size(); ++i) 
		sphComputeForce( FG, (*fluids)[i], &interactingFluids[i] );
		
	for(int i = 0; i < fluids->size(); ++i)
		integrate( FG, (*fluids)[i]->getLocalParameters(), &(*fluids)[i]->internalGetFluidParticles() );
}

void FluidSolverMultiphase::sphComputePressure(const FluidParametersGlobal &FG, FluidSph *fluid, btAlignedObjectArray<FluidSph*> *interactingFluids)
{
	BT_PROFILE("FluidSolverMultiphase::sphComputePressure()");
	
	const FluidParametersLocal &FL = fluid->getLocalParameters();
	FluidParticles &particles = fluid->internalGetFluidParticles();
	FluidSortingGrid &grid = fluid->internalGetGrid();
	
	btScalar radius = FG.m_sphSmoothRadius / FG.m_simulationScale;

	for(int i = 0; i < fluid->numParticles(); ++i)
	{
		btScalar sum = 0.0;	
		particles.m_neighborTable[i].clear();

		FluidSortingGrid::FoundCells foundCells;
		grid.findCells(particles.m_pos[i], radius, &foundCells);
		
		for(int cell = 0; cell < FluidSortingGrid::NUM_FOUND_CELLS; cell++) 
		{
			FluidGridIterator &FI = foundCells.m_iterators[cell];
			
			for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
			{
				if(i == n) continue;
				
				btVector3 distance = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
				btScalar distanceSquared = distance.length2();
				
				if(FG.m_sphRadiusSquared > distanceSquared) 
				{
					btScalar c = FG.m_sphRadiusSquared - distanceSquared;
					sum += c * c * c;
					
					if( !particles.m_neighborTable[i].isFilled() ) particles.m_neighborTable[i].addNeighbor( n, btSqrt(distanceSquared) );
					else break;
				}
			}
		}
		
		btScalar density = sum * FL.m_particleMass * FG.m_poly6KernCoeff;
		
		//EXTERNAL_FLUID_INTERACTION
		for(int j = 0; j < interactingFluids->size(); ++j)
		{
			FluidSph *externalFluid = (*interactingFluids)[j];
		
			const FluidParametersLocal &externalFL = externalFluid->getLocalParameters();
			FluidParticles &externalParticles = externalFluid->internalGetFluidParticles();
			FluidSortingGrid &externalGrid = externalFluid->internalGetGrid();
			
			btScalar externalSum = 0.0;	
			FluidSortingGrid::FoundCells externalFoundCells;
			externalGrid.findCells(particles.m_pos[i], radius, &externalFoundCells);
			for(int cell = 0; cell < FluidSortingGrid::NUM_FOUND_CELLS; cell++) 
			{
				FluidGridIterator &FI = externalFoundCells.m_iterators[cell];
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
			
			density += externalSum * externalFL.m_particleMass * FG.m_poly6KernCoeff;
		}
		//EXTERNAL_FLUID_INTERACTION
		
		particles.m_pressure[i] = (density - FL.m_restDensity) * FL.m_stiffness;
		particles.m_invDensity[i] = 1.0f / density;
	}
}

void computeForceNeighborTable_Multiphase(const FluidParametersGlobal &FG, const btScalar vterm, int particleIndex, FluidParticles *fluids)
{
	int i = particleIndex;

	btVector3 force(0, 0, 0);
	for(int j = 0; j < fluids->m_neighborTable[i].numNeighbors(); j++ ) 
	{
		int n = fluids->m_neighborTable[i].getNeighborIndex(j);
		
		btVector3 distance = (fluids->m_pos[i] - fluids->m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
		
		btScalar c = FG.m_sphSmoothRadius - fluids->m_neighborTable[i].getDistance(j);
		btScalar pterm = -0.5f * c * FG.m_spikyKernGradCoeff 
					 * ( fluids->m_pressure[i] + fluids->m_pressure[n]) / fluids->m_neighborTable[i].getDistance(j);
		btScalar dterm = c * fluids->m_invDensity[i] * fluids->m_invDensity[n];

		btVector3 forceAdded( (pterm * distance.x() + vterm * (fluids->m_vel_eval[n].x() - fluids->m_vel_eval[i].x())) * dterm,
							  (pterm * distance.y() + vterm * (fluids->m_vel_eval[n].y() - fluids->m_vel_eval[i].y())) * dterm,
							  (pterm * distance.z() + vterm * (fluids->m_vel_eval[n].z() - fluids->m_vel_eval[i].z())) * dterm );
		force += forceAdded;
	}
	
	fluids->m_sph_force[i] = force;
}

void FluidSolverMultiphase::sphComputeForce(const FluidParametersGlobal &FG, FluidSph *fluid, btAlignedObjectArray<FluidSph*> *interactingFluids)
{
	BT_PROFILE("sphComputeForce()");
	
	const FluidParametersLocal &FL = fluid->getLocalParameters();
	FluidParticles &fluids = fluid->internalGetFluidParticles();
	btScalar vterm = FG.m_viscosityKernLapCoeff * FL.m_viscosity;
	
	for(int i = 0; i < fluids.size(); ++i)
		computeForceNeighborTable_Multiphase(FG, vterm, i, &fluids);
	
	//EXTERNAL_FLUID_INTERACTION
	btScalar radius = FG.m_sphSmoothRadius / FG.m_simulationScale;
	
	for(int j = 0; j < interactingFluids->size(); ++j)
	{
		FluidSph *externalFluid = (*interactingFluids)[j];
			
		const FluidParametersLocal &externalFL = externalFluid->getLocalParameters();
		FluidParticles &externalParticles = externalFluid->internalGetFluidParticles();
		FluidSortingGrid &externalGrid = externalFluid->internalGetGrid();
		
		btScalar averagedViscosity = (FL.m_viscosity + externalFL.m_viscosity) * 0.5f;
		btScalar vterm2 = FG.m_viscosityKernLapCoeff * averagedViscosity;
		
		for(int i = 0; i < fluids.size(); ++i)
		{
			btVector3 externalForce(0, 0, 0);
			FluidSortingGrid::FoundCells externalFoundCells;
			externalGrid.findCells(fluids.m_pos[i], radius, &externalFoundCells);
			for(int cell = 0; cell < FluidSortingGrid::NUM_FOUND_CELLS; cell++) 
			{
				FluidGridIterator &FI = externalFoundCells.m_iterators[cell];
				for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
				{
					btVector3 distance = (fluids.m_pos[i] - externalParticles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
					btScalar distanceSquared = distance.length2();
					
					if(FG.m_sphRadiusSquared > distanceSquared) 
					{
						const bool USE_FLUIDS_V2_EQUATIONS = true;
						if(USE_FLUIDS_V2_EQUATIONS)
						{
							btScalar r = btSqrt(distanceSquared);
							
							btScalar c = FG.m_sphSmoothRadius - r;
							btScalar pterm = -0.5f * c * FG.m_spikyKernGradCoeff * ( fluids.m_pressure[i] + externalParticles.m_pressure[n]) / r;
							btScalar dterm = c * fluids.m_invDensity[i] * externalParticles.m_invDensity[n];
							
							btVector3 forceAdded( (pterm * distance.x() + vterm2 * (externalParticles.m_vel_eval[n].x() - fluids.m_vel_eval[i].x())) * dterm,
												  (pterm * distance.y() + vterm2 * (externalParticles.m_vel_eval[n].y() - fluids.m_vel_eval[i].y())) * dterm,
												  (pterm * distance.z() + vterm2 * (externalParticles.m_vel_eval[n].z() - fluids.m_vel_eval[i].z())) * dterm );
							externalForce += forceAdded;
						}
						else
						{
							btScalar r = btSqrt(distanceSquared);
							btScalar c = FG.m_sphSmoothRadius - r;
						
								//Non-kernel
							btScalar pressureScalar = -externalFL.m_particleMass;
							pressureScalar *= btScalar(0.5) * (fluids.m_pressure[i] + externalParticles.m_pressure[n]) * externalParticles.m_invDensity[n];
								//Kernel
							pressureScalar *= FG.m_spikyKernGradCoeff;
							pressureScalar *= (c * c) / r;

							btVector3 pressureForce = distance * pressureScalar;
							
								//Non-kernel
							btScalar viscosityScalar = averagedViscosity * externalFL.m_particleMass;
							viscosityScalar *= externalParticles.m_invDensity[n];
								//Kernel
							viscosityScalar *= FG.m_viscosityKernLapCoeff;
							viscosityScalar *= c;
							
							btVector3 viscosityForce = (externalParticles.m_vel_eval[n] - fluids.m_vel_eval[i]) * viscosityScalar;
							
							
							btVector3 forceAdded2 = pressureForce + viscosityForce;
							externalForce += forceAdded2;
						}
					}
				}
			}
			
			fluids.m_sph_force[i] += externalForce;
		}
	}
	//EXTERNAL_FLUID_INTERACTION
}
