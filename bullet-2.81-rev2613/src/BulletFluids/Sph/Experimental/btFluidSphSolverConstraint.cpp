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
#include "btFluidSphSolverConstraint.h"

#include "BulletFluids/Sph/btFluidSortingGrid.h"

#include "LinearMath/btQuickprof.h"		//BT_PROFILE(name) macro

///SPH Kernel (1 - (d/h)^2)^3, used for the calculation of density and Jacobian entries
btScalar computeKernel(btScalar simulationScaleDistance, btScalar sphSmoothRadius)
{
	btScalar d = simulationScaleDistance / sphSmoothRadius;
	btScalar kernel = btScalar(1.0) - d*d;
	return (kernel * kernel * kernel);
}

void calculateSumsInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
									const btFluidSortingGrid& grid, btFluidParticles& particles,
									btFluidSphSolverConstraint::CfParticles& cfData)
{
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	if(currentCell.m_firstIndex <= currentCell.m_lastIndex)	//if cell is not empty
	{
		btFluidSortingGrid::FoundCells foundCells;
		grid.findCellsSymmetric(particles.m_pos[currentCell.m_firstIndex], foundCells);
		
		for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
		{
			//Remove particle, with index i, from grid cell to prevent self-particle interactions
			++foundCells.m_iterators[0].m_firstIndex;	//Local cell; currentCell == foundCells.m_iterators[0]
			
			for(int cell = 0; cell < btFluidSortingGrid::NUM_FOUND_CELLS_SYMMETRIC; cell++) 
			{
				btFluidGridIterator& FI = foundCells.m_iterators[cell];
				
				for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
				{
					//Simulation-scale distance
					btVector3 difference = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;		
					btScalar distanceSquared = difference.length2();
					
					if(FG.m_sphRadiusSquared > distanceSquared)
					{
						btScalar distance = btSqrt(distanceSquared);
						
						btScalar weight = computeKernel(distance, FG.m_sphSmoothRadius);
						
						cfData.m_density[i] += weight;
						cfData.m_density[n] += weight;
						
						if( !cfData.m_neighborTable[i].isFilled() ) cfData.m_neighborTable[i].addNeighbor(n, distance);
						else if( !cfData.m_neighborTable[n].isFilled() ) cfData.m_neighborTable[n].addNeighbor(i, distance);
						else 
						{
							cell = btFluidSortingGrid::NUM_FOUND_CELLS_SYMMETRIC;	//Break out of outer loop 
							break;
						}
					}
				}
			}
		}
	}
}

void setupConstraintsInCellSymmetric(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL,
									int gridCellIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
									btFluidSphSolverConstraint::CfParticles& cfData)
{
	BT_PROFILE("setupConstraintsInCellSymmetric()");
	
	btFluidGridIterator currentCell = grid.getGridCell(gridCellIndex);
	for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
	{
		for(int j = 0; j < cfData.m_neighborTable[i].numNeighbors(); j++) 
		{
			int n = cfData.m_neighborTable[i].getNeighborIndex(j);
			btScalar distance = cfData.m_neighborTable[i].getDistance(j);
			if(distance >= FG.m_sphSmoothRadius) continue;
			
			btScalar weight = computeKernel(distance, FG.m_sphSmoothRadius);
			
			btScalar matrixDiagonal = btScalar(2.0) * weight * weight;
			
			//m_A is the matrix diagonal; as stated by Constraint Fluids, the matrix element A_ii
			//is the result of the dot product of row i of the Jacobian with itself - hence the 'weight * weight' term
			//	
			//Regarding the Jacobian matrix, J:
			//	-Each row contains multiple blocks, each block being a (1 x 3) row vector
			//	-Non-diagonal block of a row i: (mass / restDensity) * kernel_in * direction_in
			//		kernel_in is the kernel result between particles i and n
			//		direction_in is the normal vector from n to i(not verified)
			//			-Since a normal vector dot itself == 1, it does not need to be included in the calculation
			//	-Diagonal block is the sum of all other elements on that row(this is the reason for the 'btScalar(2.0)' term)
			//
			//	This differs from the constraint fluid paper in that it uses a different kernel,
			//	and does not include the constant mass / restDensity^2 term (squared since the linear system is (J * J^T) * lambda = b)
			cfData.m_A[i] += matrixDiagonal;
			cfData.m_A[n] += matrixDiagonal;
		}
	}
}

btScalar calculateDpSum(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL,
						int particleIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
						btFluidSphSolverConstraint::CfParticles& cfData)
{
	btFluidSortingGrid::FoundCells foundCells;
	grid.findCellsSymmetric(particles.m_pos[particleIndex], foundCells);
	
	int i = particleIndex;
	
	btScalar dpSum(0.0);
	
	for(int cell = 0; cell < btFluidSortingGrid::NUM_FOUND_CELLS_SYMMETRIC; cell++) 
	{
		const btFluidGridIterator& FI = foundCells.m_iterators[cell];
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
		{
			btVector3 difference = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;	//Simulation-scale distance
			btScalar distanceSquared = difference.length2();
			
			if(FG.m_sphRadiusSquared > distanceSquared && i != n)
			{
				btScalar distance = btSqrt(distanceSquared);
				btVector3 direction = difference / distance;
				
				btScalar weight = computeKernel(distance, FG.m_sphSmoothRadius);
				
				//Compute the J*v term of the column vector b(right side of the equation)
				//Each row of J (1 x 3N) consists of several row vectors(1 x 3), each row vector being:
				//	(mass / restDensity) * kernel_in * direction_in
				//	--Mass and rest density are assumed to be 1, and as a result is omitted here
				//
				//If following the 'Constraint Fluids' paper, this should be velocity[i] - velocity[n]
				btScalar dv = direction.dot(particles.m_vel[n] - particles.m_vel[i]);
				
				dpSum += weight * dv;
			}
		}
	}
	
	return dpSum;
}

void applyImpulses(const btFluidSphParametersGlobal& FG, const btFluidSphParametersLocal& FL, btScalar magnitude,
					int particleIndex, const btFluidSortingGrid& grid, btFluidParticles& particles,
					btFluidSphSolverConstraint::CfParticles& cfData)
{
	btFluidSortingGrid::FoundCells foundCells;
	grid.findCellsSymmetric(particles.m_pos[particleIndex], foundCells);
	
	int i = particleIndex;
	
	for(int cell = 0; cell < btFluidSortingGrid::NUM_FOUND_CELLS_SYMMETRIC; cell++) 
	{
		const btFluidGridIterator& FI = foundCells.m_iterators[cell];
		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
		{
			btVector3 difference = (particles.m_pos[i] - particles.m_pos[n]) * FG.m_simulationScale;	//Simulation-scale distance
			btScalar distanceSquared = difference.length2();
			
			if(FG.m_sphRadiusSquared > distanceSquared && i != n)
			{
				btScalar distance = btSqrt(distanceSquared);
				btVector3 direction = difference / distance;
				
				btScalar weight = computeKernel(distance, FG.m_sphSmoothRadius);
				
				//Directly apply impulses, as opposed to solving for lambda(v_next = v_prev + J^T * lambda)
				//Per-neighbor entry in J is: (mass / restDensity) * kernel_in * direction_in
				btVector3 impulse_in = direction * magnitude * weight;
				
				particles.m_vel[i] += impulse_in;
				particles.m_vel[n] += -impulse_in;
			}
		}
	}
}


void btFluidSphSolverConstraint::updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids)
{
	BT_PROFILE("btFluidSphSolverConstraint::updateGridAndCalculateSphForces()");
	
	//SPH data is discarded/recalculated every frame, so only 1
	//set of arrays are needed if there is no fluid-fluid interaction.
	if( m_cfData.size() != 1 ) m_cfData.resize(1);
	
	for(int fluidIndex = 0; fluidIndex < numFluids; ++fluidIndex)
	{
		btFluidSph* fluid = fluids[fluidIndex];
		int numParticles = fluid->numParticles();
		if(!numParticles) continue;
		
		const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
		const btFluidSortingGrid& grid = fluid->getGrid();
		btFluidParticles& particles = fluid->internalGetParticles();
		btFluidSphSolverConstraint::CfParticles& cfData = m_cfData[0];
		
		if( numParticles > cfData.size() ) cfData.resize(numParticles);
		
		fluid->insertParticlesIntoGrid();
		
		//Compute density, and load neighbor tables
		{
			for(int i = 0; i < numParticles; ++i) cfData.m_neighborTable[i].clear();
			
			for(int i = 0; i < numParticles; ++i) cfData.m_density[i] = btScalar(0.0);
			
			for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
			{
				const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
				
				for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
					calculateSumsInCellSymmetric(FG, multithreadingGroup[cell], grid, particles, cfData);
			}
		}
		
		//Setup constraints for each particle
		{
			//Current implementation assumes that rest density is always 1.0; must modify Jacobian to use different rest density
			const btScalar REST_DENSITY(1.0);
			const btScalar BAUMGARTE(0.1);
			
			//This differs from the GDC2014 presentation, which uses (FL.m_restDensity - cfData.m_density[i])
			for(int i = 0; i < numParticles; ++i) cfData.m_bias[i] = (cfData.m_density[i] - REST_DENSITY) * BAUMGARTE;
			
			for(int i = 0; i < numParticles; ++i) cfData.m_A[i] = btScalar(0.0);
			
			for(int group = 0; group < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++group)
			{
				const btAlignedObjectArray<int>& multithreadingGroup = grid.internalGetMultithreadingGroup(group);
				
				for(int cell = 0; cell < multithreadingGroup.size(); ++cell)
					setupConstraintsInCellSymmetric(FG, FL, multithreadingGroup[cell], grid, particles, cfData);
			}
		}
		
		//Solve density constraints, using Projected Gauss Seidel
		const int NUM_ITERATIONS = 5;
		for(int iteration = 0; iteration < NUM_ITERATIONS; ++iteration)
		{
			for(int i = 0; i < numParticles; ++i)
			{
				//dpSum is positive if other particles are moving towards particle i(density is increasing)
				//(if following the 'Constraint Fluids' paper, should be negative in case above -- here, it is negated in calculateDpSum())
				btScalar dpSum = calculateDpSum(FG, FL, i, grid, particles, cfData);
				
				btScalar target = dpSum + cfData.m_bias[i];		//This is an element of the 'b' vector of the equation Ax = b
				
				//Calculate lambda_i using Gauss Seidel / Jacobi: 
				//lambda_i = (b_i - A_i * lambda) / A_ii
				//
				//Since each iteration assumes that the lambda/magnitude vector is all 0s, we only need to compute (b_i / A_ii)
				//The velocities(on the right side of the equation) are updated directly, so lambda remains 0
				btScalar magnitude = target / cfData.m_A[i];
				
				//Allow only positive(repulsive) impulses; this is not necessary 
				//if using the exact same method as in the paper, 
				//but is needed for stability in the current implementation
				magnitude = btMax( btScalar(0.0), magnitude );
				
				applyImpulses(FG, FL, magnitude, i, grid, particles, cfData);
			}
		}
	}
}
