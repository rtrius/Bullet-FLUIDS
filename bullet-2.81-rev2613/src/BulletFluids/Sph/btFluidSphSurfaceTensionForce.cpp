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
		distance = (distance < SIMD_EPSILON) ? SIMD_EPSILON : distance;
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;
		btScalar spikyKernelGradientScalar = closeness * closeness;
		btVector3 spikyKernelGradient_in = n_to_i * (spikyKernelGradientScalar / distance);
		
		colorFieldGradient[i] += spikyKernelGradient_in * invDensity[n];
		colorFieldGradient[n] += spikyKernelGradient_in * invDensity[i];
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
	
	for(int i = 0; i < numParticles; ++i) m_colorFieldGradient[i] *= FG.m_sphSmoothRadius * FL.m_sphParticleMass * FG.m_spikyKernGradCoeff;  
}


void computeSurfaceTensionForceNeighborTableSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, 
													int particleIndex, const btFluidParticles& particles,
													const btAlignedObjectArray<btFluidSphNeighbors>& neighborTables, 
													const btAlignedObjectArray<btScalar>& invDensity,
													const btAlignedObjectArray<btVector3>& colorFieldGradient,
													btAlignedObjectArray<btVector3>& out_surfaceTensionForce)
{
	int i = particleIndex;

	for(int j = 0; j < neighborTables[i].numNeighbors(); j++) 
	{
		int n = neighborTables[i].getNeighborIndex(j);
		
		btScalar distance = neighborTables[i].getDistance(j);
		distance = (distance < SIMD_EPSILON) ? SIMD_EPSILON : distance;
		
		btScalar closeness = FG.m_sphSmoothRadius - distance;

		const btScalar C_coefficient = btScalar(32.0) / ( SIMD_PI * btPow(FG.m_sphSmoothRadius, btScalar(9.0)) );
		btScalar cubedCloseness = closeness*closeness*closeness;
		btScalar cubedDistance = distance*distance*distance;

		btScalar C_partialResult = cubedCloseness * cubedDistance;
		if(distance < btScalar(0.5) * FG.m_sphSmoothRadius)
		{
			C_partialResult *= btScalar(2.0);
			C_partialResult -= btPow(FG.m_sphSmoothRadius, btScalar(6.0)) / btScalar(64.0);
		}

		btScalar C = C_coefficient * C_partialResult;
		
		btVector3 n_to_i = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		//Simulation-scale distance
		btVector3 normal_n_to_i = n_to_i / distance;
		
		btVector3 cohesionForce = -normal_n_to_i * (C * FL.m_sphParticleMass * FL.m_sphParticleMass);
		btVector3 curvatureForce = -FL.m_sphParticleMass * (colorFieldGradient[i] - colorFieldGradient[n]);

		btScalar density_i = btScalar(1.0) / invDensity[i];
		btScalar density_n = btScalar(1.0) / invDensity[n];
		btScalar K = (btScalar(2.0) * FL.m_restDensity) / (density_i + density_n);

		btVector3 surfaceTensionForce_in = (cohesionForce + curvatureForce) * K;
		
		out_surfaceTensionForce[i] += surfaceTensionForce_in;
		out_surfaceTensionForce[n] += -surfaceTensionForce_in;
	}
}
void computeSurfaceTensionForceInCellSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, 
											int gridCellIndex, const btFluidSortingGrid& grid,  const btFluidParticles& particles,
											const btAlignedObjectArray<btFluidSphNeighbors>& neighborTables, 
											const btAlignedObjectArray<btScalar>& invDensity,
											const btAlignedObjectArray<btVector3>& colorFieldGradient,
											btAlignedObjectArray<btVector3>& out_surfaceTensionForce)
{
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int particleIndex = currentCell.m_firstIndex; particleIndex <= currentCell.m_lastIndex; ++particleIndex)
	{
		computeSurfaceTensionForceNeighborTableSymmetric(FG, FL, particleIndex, particles, neighborTables, invDensity, 
														colorFieldGradient, out_surfaceTensionForce);
	}
}

void btFluidSphSurfaceTensionForce::computeSurfaceTensionForce(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, 
															const btAlignedObjectArray<btFluidSphNeighbors>& neighborTables, 
															const btAlignedObjectArray<btScalar>& invDensity)
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
				computeSurfaceTensionForceInCellSymmetric(FG, FL, multithreadingGroup[cell], grid, particles, neighborTables, invDensity,
															m_colorFieldGradient, m_surfaceTensionForce);
		}
	}
	
	const btScalar SURFACE_TENSION_STRENGTH(300.0);
	for(int i = 0; i < numParticles; ++i) m_surfaceTensionForce[i] *= SURFACE_TENSION_STRENGTH;
}
