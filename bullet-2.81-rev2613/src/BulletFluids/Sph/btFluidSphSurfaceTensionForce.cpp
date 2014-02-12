/*
Bullet-FLUIDS 
Copyright (c) 2012-2014 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/

#include "btFluidSphSurfaceTensionForce.h"

#include "btFluidSphSolver.h"

#include "btFluidSortingGrid.h"
#include "btFluidSph.h"

//Issues with the surface tension force:
//- Not stable for strong surface tension values(e.g. 1.0). (may be related to low stiffness/compressibility or time step too high)
//- Particles inside the fluid should have a colorFieldGradient value near 0, but this is not the case currently.
//- Paper does not specify which kernel is used for the colorFieldGradient; should not be an issue as long
//as the gradient balances out inside the fluid, but may be related to the colorFieldGradient issue on the above line.
#define USE_SPIKY_KERNEL_GRADIENT

void computeColorFieldGradientNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, int particleIndex, 
												const btFluidParticles& particles,
												const btAlignedObjectArray<btFluidSphNeighbors>& neighborTables, 
												const btAlignedObjectArray<btScalar>& invDensity,
												btAlignedObjectArray<btVector3>& colorFieldGradient)
{
	int i = particleIndex;

	for(int j = 0; j < neighborTables[i].numNeighbors(); j++) 
	{
		int n = neighborTables[i].getNeighborIndex(j);
		
		btVector3 n_to_i = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
		btScalar distance = neighborTables[i].getDistance(j);
		//distance = (distance < SIMD_EPSILON) ? SIMD_EPSILON : distance;
		
#ifdef USE_SPIKY_KERNEL_GRADIENT
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		btScalar spikyKernelGradientScalar = closeness * closeness;
		btVector3 spikyKernelGradient_in = n_to_i * (spikyKernelGradientScalar / distance);
		
		colorFieldGradient[i] += spikyKernelGradient_in * invDensity[n];
		colorFieldGradient[n] += spikyKernelGradient_in * invDensity[i];
#else
		btScalar squaredDistance = distance * distance;
		btScalar squaredCloseness = FG.m_sphRadiusSquared - squaredDistance;
		
		btScalar poly6KernelGradientScalar = squaredCloseness * squaredCloseness;
		btVector3 poly6KernelGradient_in = n_to_i * (poly6KernelGradientScalar);
		colorFieldGradient[i] += poly6KernelGradient_in * invDensity[n];
		colorFieldGradient[n] += poly6KernelGradient_in * invDensity[i];
#endif
	}
}
void computeColorFieldGradientInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
											const btFluidSortingGrid& grid, const btFluidParticles& particles,
											const btAlignedObjectArray<btFluidSphNeighbors>& neighborTables, 
											const btAlignedObjectArray<btScalar>& invDensity,
												btAlignedObjectArray<btVector3>& colorFieldGradient)
{
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computeColorFieldGradientNeighborTableSymmetric(FG, particleIndex, particles, neighborTables, invDensity, colorFieldGradient);
	}
}
void btFluidSphSurfaceTensionForce::computeColorFieldGradient(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, 
														const btAlignedObjectArray<btFluidSphNeighbors>& neighborTables, 
														const btAlignedObjectArray<btScalar>& invDensity)
{
	BT_PROFILE("btFluidSphSolverDefault::computeColorFieldGradient()");
	
	int numParticles = fluid->numParticles();
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	const btFluidSortingGrid& grid = fluid->getGrid();
	btFluidParticles& particles = fluid->internalGetParticles();
	
	for(int i = 0; i < numParticles; ++i) m_colorFieldGradient[i].setValue(0, 0, 0);
	
	for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
	{
		const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
		if( !multithreadingGroup.size() ) continue;
		
		//If using multiple threads, this needs to be moved to a separate function
		{
			for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
				computeColorFieldGradientInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, 
														neighborTables, invDensity, m_colorFieldGradient);
		}
	}
	
	
#ifdef USE_SPIKY_KERNEL_GRADIENT
	for(int i = 0; i < numParticles; ++i) m_colorFieldGradient[i] *= FG.m_sphSmoothRadius * FL.m_sphParticleMass * FG.m_spikyKernGradCoeff;
#else
	const btScalar poly6KernGradCoeff = btScalar(-945.0) / ( btScalar(32.0) * SIMD_PI * btPow(FG.m_sphSmoothRadius, 9) );
	for(int i = 0; i < numParticles; ++i) m_colorFieldGradient[i] *= FG.m_sphSmoothRadius * FL.m_sphParticleMass * poly6KernGradCoeff;
#endif
}


void computeSurfaceTensionForceNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, 
													const btFluidSphSurfaceTensionForce::Coefficients& ST,
													int particleIndex, const btFluidParticles& particles,
													const btAlignedObjectArray<btFluidSphNeighbors>& neighborTables, 
													const btAlignedObjectArray<btScalar>& density,
													const btAlignedObjectArray<btVector3>& colorFieldGradient,
													btAlignedObjectArray<btVector3>& out_surfaceTensionForce)
{
	int i = particleIndex;

	for(int j = 0; j < neighborTables[i].numNeighbors(); j++) 
	{
		int n = neighborTables[i].getNeighborIndex(j);
		
		btScalar distance = neighborTables[i].getDistance(j);
		//distance = (distance < SIMD_EPSILON) ? SIMD_EPSILON : distance;
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		
		btScalar cubedCloseness = closeness*closeness*closeness;
		btScalar cubedDistance = distance*distance*distance;

		btScalar C_partialResult = cubedCloseness * cubedDistance;
		
		if(distance < ST.m_sphSmoothRadiusHalved)
		{
			C_partialResult *= btScalar(2.0);
			C_partialResult -= ST.m_C_add_coeff;
		}

		btScalar C = ST.m_C_multiply_coeff * C_partialResult;
		
		btVector3 n_to_i = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
		btVector3 normal_n_to_i = n_to_i / distance;
		
		btVector3 cohesionForce = normal_n_to_i * (C * FL.m_sphParticleMass);
		btVector3 curvatureForce = colorFieldGradient[i] - colorFieldGradient[n];
		
		btScalar K = -btScalar(1.0) / (density[i] + density[n]);

		btVector3 surfaceTensionForce_in = (cohesionForce + curvatureForce) * K;
		
		out_surfaceTensionForce[i] += surfaceTensionForce_in;
		out_surfaceTensionForce[n] += -surfaceTensionForce_in;
	}
}
void computeSurfaceTensionForceInCellSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, 
											const btFluidSphSurfaceTensionForce::Coefficients& ST, int gridCellIndex, 
											const btFluidSortingGrid& grid,  const btFluidParticles& particles,
											const btAlignedObjectArray<btFluidSphNeighbors>& neighborTables, 
											const btAlignedObjectArray<btScalar>& density,
											const btAlignedObjectArray<btVector3>& colorFieldGradient,
											btAlignedObjectArray<btVector3>& out_surfaceTensionForce)
{
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computeSurfaceTensionForceNeighborTableSymmetric(FG, FL, ST, particleIndex, particles, neighborTables, density, 
														colorFieldGradient, out_surfaceTensionForce);
	}
}

void btFluidSphSurfaceTensionForce::computeSurfaceTensionForce(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, 
															const btFluidSphSurfaceTensionForce::Coefficients& ST,
															const btAlignedObjectArray<btFluidSphNeighbors>& neighborTables, 
															const btAlignedObjectArray<btScalar>& density)
{
	BT_PROFILE("btFluidSphSolverDefault::computeSurfaceTensionForce()");
	
	int numParticles = fluid->numParticles();
	const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	const btFluidSortingGrid& grid = fluid->getGrid();
	btFluidParticles& particles = fluid->internalGetParticles();
	
	for(int i = 0; i < numParticles; ++i) m_surfaceTensionForce[i] = btVector3(0, 0, 0);

	for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
	{
		const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
		if( !multithreadingGroup.size() ) continue;
		
		//If using multiple threads, this needs to be moved to a separate function
		{
			for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
				computeSurfaceTensionForceInCellSymmetric(FG, FL, ST, multithreadingGroup[cell], grid, particles, neighborTables, density,
															m_colorFieldGradient, m_surfaceTensionForce);
		}
	}
	
	for(int i = 0; i < numParticles; ++i) 
		m_surfaceTensionForce[i] *= FL.m_surfaceTension * FL.m_sphParticleMass * (btScalar(2.0) * FL.m_restDensity);
}
