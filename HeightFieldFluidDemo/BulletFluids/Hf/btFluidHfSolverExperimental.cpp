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
#include "btFluidHfSolverExperimental.h"

#include "LinearMath/btMinMax.h"

btScalar bilinearInterpolateWithBounds(const btFluidColumns& columns, const btAlignedObjectArray<btScalar>& array, 
											btScalar iPos, btScalar jPos, btScalar& out_min, btScalar& out_max)
{
	int i = (int)iPos;
	int j = (int)jPos;

	btScalar iParam1 = iPos - btScalar(i);
	btScalar iParam0 = btScalar(1.0) - iParam1;

	btScalar jParam1 = jPos - btScalar(j);
	btScalar jParam0 = btScalar(1.0) - jParam1;

	btScalar SW = array[columns.getIndex(i, j)];
	btScalar SE = array[columns.getIndex(i+1, j)];
	btScalar NW = array[columns.getIndex(i, j+1)];
	btScalar NE = array[columns.getIndex(i+1, j+1)];

	out_min = btMin( btMin(SW, NW), btMin(SE, NE) );
	out_max = btMax( btMax(SW, NW), btMax(SE, NE) );
	
	btScalar a = jParam0 * SW + jParam1 * NW;
	btScalar b = jParam0 * SE + jParam1 * NE;
	return iParam0 * a + iParam1 * b;
}

btScalar advectForwardWithBounds(const btFluidColumns& columns, const btAlignedObjectArray<btScalar>& array, 
						btScalar i, btScalar j, btScalar vel_x, btScalar vel_z, btScalar dt, btScalar& out_min, btScalar& out_max)
{
	// trace particle backwards in time
	btScalar srcI = i - vel_x * dt * columns.m_gridCellWidthInv;
	btScalar srcJ = j - vel_z * dt * columns.m_gridCellWidthInv;

	// srcI and srcJ are indices into the array,
	// we need to clamp them to within the domain
	srcI = btMax (srcI, btScalar(0.0));
	srcI = btMin (srcI, btScalar(columns.m_numNodesX-1));
	srcJ = btMax (srcJ, btScalar(0.0));
	srcJ = btMin (srcJ, btScalar(columns.m_numNodesZ-1));

	return bilinearInterpolateWithBounds(columns, array, srcI, srcJ, out_min, out_max);
}
btScalar btFluidHfSolverExperimental::advectReverse(const btFluidColumns& columns, const btAlignedObjectArray<btScalar>& array, 
													btScalar i, btScalar j, btScalar vel_x, btScalar vel_z, btScalar dt)
{
	// trace particle forward in time
	btScalar srcI = i + vel_x * dt * columns.m_gridCellWidthInv;
	btScalar srcJ = j + vel_z * dt * columns.m_gridCellWidthInv;

	// srcI and srcJ are indices into the array,
	// we need to clamp them to within the domain
	srcI = btMax (srcI, btScalar(0.0));
	srcI = btMin (srcI, btScalar(columns.m_numNodesX-1));
	srcJ = btMax (srcJ, btScalar(0.0));
	srcJ = btMin (srcJ, btScalar(columns.m_numNodesZ-1));

	return bilinearInterpolate(columns, array, srcI, srcJ);
}

void btFluidHfSolverExperimental::advectFluid_macCormack(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	//Forward advect
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index]) continue;

			btScalar vel_x = hfParameters.m_globalVelocityX;
			btScalar vel_z = hfParameters.m_globalVelocityZ;

			vel_x += (columns.m_vel_x[index]+columns.m_vel_x[index+1]) * btScalar(0.5);
			vel_z += (columns.m_vel_z[index]+columns.m_vel_z[index+columns.m_numNodesX]) * btScalar(0.5);

			m_tempForwardAdvect[index] = advectForwardWithBounds(columns, columns.m_fluidDepth, btScalar(i), btScalar(j), vel_x, vel_z, dt,
															m_tempForwardAdvectMin[index], m_tempForwardAdvectMax[index]);
		}
	}
	
	//Reverse advect to get error
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index]) continue;

			btScalar vel_x = hfParameters.m_globalVelocityX;
			btScalar vel_z = hfParameters.m_globalVelocityZ;

			vel_x += (columns.m_vel_x[index]+columns.m_vel_x[index+1]) * btScalar(0.5);
			vel_z += (columns.m_vel_z[index]+columns.m_vel_z[index+columns.m_numNodesX]) * btScalar(0.5);

			m_tempReverseAdvect[index] = advectReverse(columns, m_tempForwardAdvect, btScalar(i), btScalar(j), vel_x, vel_z, dt);
		}
	}
	
	//Combine forward advect and reverse advect to get result
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index]) continue;
			
			//Error is calculated by comparing the original fluidDepth with the forward-reverse advected fluidDepth
			//As the fluid is advected twice, the error is also doubled
			btScalar error = (columns.m_fluidDepth[index] - m_tempReverseAdvect[index]) * btScalar(0.5);
			btScalar advectedValue = m_tempForwardAdvect[index] + error;
			
			//Revert to semi-Lagrangian if result is not within bounds of values used for bilinear interpolation
			btScalar min = m_tempForwardAdvectMin[index];
			btScalar max = m_tempForwardAdvectMax[index];
			columns.m_fluidDepth[index] = (min <= advectedValue && advectedValue <= max) ? advectedValue : m_tempForwardAdvect[index];
		}
	}
}
void btFluidHfSolverExperimental::advectVelocityX_macCormack(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	//Forward advect
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index]) continue;

			btScalar vel_x = hfParameters.m_globalVelocityX;
			btScalar vel_z = hfParameters.m_globalVelocityZ;

			vel_x += columns.m_vel_x[index];
			vel_z += (columns.m_vel_z[index]+columns.m_vel_z[index-1]+columns.m_vel_z[index+columns.m_numNodesX]+columns.m_vel_z[index+columns.m_numNodesX-1]) * btScalar(0.25);

			m_tempForwardAdvect[index] = advectForwardWithBounds(columns, columns.m_vel_x, btScalar(i), btScalar(j), vel_x, vel_z, dt,
															m_tempForwardAdvectMin[index], m_tempForwardAdvectMax[index]);
		}
	}
	
	//Reverse advect to get error
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index]) continue;

			btScalar vel_x = hfParameters.m_globalVelocityX;
			btScalar vel_z = hfParameters.m_globalVelocityZ;

			vel_x += columns.m_vel_x[index];
			vel_z += (columns.m_vel_z[index]+columns.m_vel_z[index-1]+columns.m_vel_z[index+columns.m_numNodesX]+columns.m_vel_z[index+columns.m_numNodesX-1]) * btScalar(0.25);

			m_tempReverseAdvect[index] = advectReverse(columns, m_tempForwardAdvect, btScalar(i), btScalar(j), vel_x, vel_z, dt);
		}
	}
	
	//Combine forward advect and reverse advect to get result
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index]) continue;
			
			//As the fluid is advected twice, the error is also doubled
			btScalar error = (columns.m_vel_x[index] - m_tempReverseAdvect[index]) * btScalar(0.5);
			btScalar advectedValue = m_tempForwardAdvect[index] + error;
			
			//Revert to semi-Lagrangian if result is not within bounds of values used for bilinear interpolation
			btScalar min = m_tempForwardAdvectMin[index];
			btScalar max = m_tempForwardAdvectMax[index];
			columns.m_vel_x[index] = (min <= advectedValue && advectedValue <= max) ? advectedValue : m_tempForwardAdvect[index];
		}
	}
}
void btFluidHfSolverExperimental::advectVelocityZ_macCormack(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{

	//Forward advect
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index]) continue;

			btScalar vel_x = hfParameters.m_globalVelocityX;
			btScalar vel_z = hfParameters.m_globalVelocityZ;

			vel_x += (columns.m_vel_x[index]+columns.m_vel_x[index+1]+columns.m_vel_x[index-columns.m_numNodesX]+columns.m_vel_x[index-columns.m_numNodesX+1]) * btScalar(0.25);
			vel_z += columns.m_vel_z[index];

			m_tempForwardAdvect[index] = advectForwardWithBounds(columns, columns.m_vel_z, btScalar(i), btScalar(j), vel_x, vel_z, dt,
															m_tempForwardAdvectMin[index], m_tempForwardAdvectMax[index]);
		}
	}
	
	//Reverse advect to get error
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index]) continue;

			btScalar vel_x = hfParameters.m_globalVelocityX;
			btScalar vel_z = hfParameters.m_globalVelocityZ;

			vel_x += (columns.m_vel_x[index]+columns.m_vel_x[index+1]+columns.m_vel_x[index-columns.m_numNodesX]+columns.m_vel_x[index-columns.m_numNodesX+1]) * btScalar(0.25);
			vel_z += columns.m_vel_z[index];

			m_tempReverseAdvect[index] = advectReverse(columns, m_tempForwardAdvect, btScalar(i), btScalar(j), vel_x, vel_z, dt);
		}
	}
	
	//Combine forward advect and reverse advect to get result
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index]) continue;
			
			//As the fluid is advected twice, the error is also doubled
			btScalar error = (columns.m_vel_z[index] - m_tempReverseAdvect[index]) * btScalar(0.5);
			btScalar advectedValue = m_tempForwardAdvect[index] + error;
			
			//Revert to semi-Lagrangian if result is not within bounds of values used for bilinear interpolation
			btScalar min = m_tempForwardAdvectMin[index];
			btScalar max = m_tempForwardAdvectMax[index];
			columns.m_vel_z[index] = (min <= advectedValue && advectedValue <= max) ? advectedValue : m_tempForwardAdvect[index];
		}
	}
}


void btFluidHfSolverExperimental::advectVelocityX_new(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex (i, j);
			if (!columns.m_active[index]) continue;
			
			btScalar vel_x = hfParameters.m_globalVelocityX;
			btScalar vel_z = hfParameters.m_globalVelocityZ;

			vel_x += columns.m_vel_x[index];
			vel_z += (columns.m_vel_z[index]+columns.m_vel_z[index-1]+columns.m_vel_z[index+columns.m_numNodesX]+columns.m_vel_z[index+columns.m_numNodesX-1]) * btScalar(0.25);

			columns.m_temp[index] = advect(columns, columns.m_vel_x, btScalar(i), btScalar(j), vel_x, vel_z, dt);
		}
	}

	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex (i, j);
			columns.m_vel_x[index] = columns.m_temp[index];
		}
	}
}
void btFluidHfSolverExperimental::advectVelocityZ_new(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex (i, j);
			if (!columns.m_active[index]) continue;
			
			btScalar vel_x = hfParameters.m_globalVelocityX;
			btScalar vel_z = hfParameters.m_globalVelocityZ;

			vel_x += (columns.m_vel_x[index]+columns.m_vel_x[index+1]+columns.m_vel_x[index-columns.m_numNodesX]+columns.m_vel_x[index-columns.m_numNodesX+1]) * btScalar(0.25);
			vel_z += columns.m_vel_z[index];

			columns.m_temp[index] = advect(columns, columns.m_vel_z, btScalar(i), btScalar(j), vel_x, vel_z, dt);
		}
	}
	
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex (i, j);
			columns.m_vel_z[index] = columns.m_temp[index];
		}
	}
}

void btFluidHfSolverExperimental::updateHeight_new(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	const btScalar dampingAndTimeStep = hfParameters.m_heightDamping * dt;

	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index])
			{
				columns.m_combinedHeight[index] = columns.m_ground[index] + columns.m_fluidDepth[index];
				continue;
			}
			
			btScalar u1, u2, v1, v2;
			if(1)
			{
				u1 = columns.m_vel_x[index+1] * (columns.m_fluidDepth[index] + columns.m_fluidDepth[index+1]) * btScalar(0.5);
				u2 = columns.m_vel_x[index] * (columns.m_fluidDepth[index] + columns.m_fluidDepth[index-1]) * btScalar(0.5);
				v1 = columns.m_vel_z[index+columns.m_numNodesX] * (columns.m_fluidDepth[index] + columns.m_fluidDepth[index+columns.m_numNodesX]) * btScalar(0.5);
				v2 = columns.m_vel_z[index] * (columns.m_fluidDepth[index] + columns.m_fluidDepth[index-columns.m_numNodesX]) * btScalar(0.5);
			}
			else
			{
				//Evaluate height in upwind direction(unstable)
				u1 = columns.m_vel_x[index+1] * ( columns.m_vel_x[index+1] < btScalar(0.0) ) ? columns.m_fluidDepth[index+1] : columns.m_fluidDepth[index];
				u2 = columns.m_vel_x[index] * ( columns.m_vel_x[index] < btScalar(0.0) ) ? columns.m_fluidDepth[index] : columns.m_fluidDepth[index-1];
				v1 = columns.m_vel_z[index+columns.m_numNodesX] * ( columns.m_vel_z[index+columns.m_numNodesX] < btScalar(0.0) ) ? columns.m_fluidDepth[index+columns.m_numNodesX] : columns.m_fluidDepth[index];
				v2 = columns.m_vel_z[index] * ( columns.m_vel_z[index] < btScalar(0.0) ) ? columns.m_fluidDepth[index] : columns.m_fluidDepth[index-columns.m_numNodesX];
			}
			
			btScalar fluidDelta = ((u1 - u2) + (v1 - v2)) * (-columns.m_gridCellWidthInv);
			columns.m_fluidDepth[index] += fluidDelta * dampingAndTimeStep;
			columns.m_fluidDepth[index] = btMax( columns.m_fluidDepth[index], btScalar(0.0) );
			
			columns.m_combinedHeight[index] = columns.m_ground[index] + columns.m_fluidDepth[index];
		}
	}
}

void btFluidHfSolverExperimental::applyBoundaryConditions(btFluidColumns& columns)
{
	const btScalar ACTIVE_CELL_EPSILON = btScalar(0.0001) * columns.m_gridCellWidth;
	
	//
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			
			if( (columns.m_fluidDepth[index-1] <= ACTIVE_CELL_EPSILON && columns.m_ground[index-1] > columns.m_combinedHeight[index])
			 || (columns.m_fluidDepth[index] <= ACTIVE_CELL_EPSILON && columns.m_ground[index] > columns.m_combinedHeight[index-1]) )
			{
				columns.m_vel_x[index] = btScalar(0.0);
			}
			
			if( (columns.m_fluidDepth[index-columns.m_numNodesX] <= ACTIVE_CELL_EPSILON && columns.m_ground[index-columns.m_numNodesX] > columns.m_combinedHeight[index])
			 || (columns.m_fluidDepth[index] <= ACTIVE_CELL_EPSILON && columns.m_ground[index] > columns.m_combinedHeight[index-columns.m_numNodesX]) )
			{
				columns.m_vel_z[index] = btScalar(0.0);
			}
		}
	
	//Edges
	const btScalar BOUNDARY_FLUID_DEPTH = btScalar(0.0);
	for(int j = 0; j < columns.m_numNodesZ; j++)
	{
		int indexL = columns.getIndex(0, j);						//Left(X-)
		columns.m_fluidDepth[indexL] = BOUNDARY_FLUID_DEPTH;
		columns.m_combinedHeight[indexL] = columns.m_ground[indexL];
		
		columns.m_vel_x[indexL+1] = btScalar(0.0);
	}
	for(int j = 0; j < columns.m_numNodesZ; j++)
	{
		int indexR = columns.getIndex(columns.m_numNodesX-1, j);	//Right(X+)
		columns.m_fluidDepth[indexR] = BOUNDARY_FLUID_DEPTH;
		columns.m_combinedHeight[indexR] = columns.m_ground[indexR];
		
		columns.m_vel_x[indexR] = btScalar(0.0);
	}
	for(int i = 0; i < columns.m_numNodesX; i++)
	{
		int indexB = columns.getIndex(i, 0);						//Bottom(Z-)
		columns.m_fluidDepth[indexB] = BOUNDARY_FLUID_DEPTH;
		columns.m_combinedHeight[indexB] = columns.m_ground[indexB];
		
		columns.m_vel_z[indexB+columns.m_numNodesX] = btScalar(0.0);
	}
	for(int i = 0; i < columns.m_numNodesX; i++)
	{
		int indexT = columns.getIndex(i, columns.m_numNodesZ-1);	//Top(Z+)
		columns.m_fluidDepth[indexT] = BOUNDARY_FLUID_DEPTH;
		columns.m_combinedHeight[indexT] = columns.m_ground[indexT];
		
		columns.m_vel_z[indexT] = btScalar(0.0);
	}
}

void btFluidHfSolverExperimental::transferDisplaced_new(const btScalar volumeDisplacementScale, btFluidColumns& columns)
{
	int currentDisplacedIndex = columns.m_displacedIndex;
	int nextFrameDisplacedIndex = (columns.m_displacedIndex + 1) % 2;
	
	for (int i = 2; i < columns.m_numNodesX - 2; i++)
	{
		for (int j = 2; j < columns.m_numNodesZ - 2; j++)
		{
			int center = columns.getIndex(i, j);
			int left = columns.getIndex(i-1, j);
			int right = columns.getIndex(i+1, j);
			int down = columns.getIndex(i, j-1);
			int up = columns.getIndex(i, j+1);
		
			btScalar displacement = columns.m_displaced[currentDisplacedIndex][center] - columns.m_displaced[nextFrameDisplacedIndex][center];
			displacement /= columns.m_gridCellWidth * columns.m_gridCellWidth;	//displacement is in volume, but we want to change the height
			displacement *= volumeDisplacementScale;
			
			btScalar quarterDisplacement = displacement * btScalar(0.25);
			columns.m_fluidDepth[left] += quarterDisplacement;
			columns.m_fluidDepth[right] += quarterDisplacement;
			columns.m_fluidDepth[down] += quarterDisplacement;
			columns.m_fluidDepth[up] += quarterDisplacement;
			
			//Displacement becomes negative when an object is removed from the fluid(fluid is removed).
			//Since the fluid height is clamped to above 0, a negative fluid depth here
			//indicates that the volume has incorrectly increased.
			if( columns.m_fluidDepth[left] < -SIMD_EPSILON )printf("negative displaced volume 1\n");
			if( columns.m_fluidDepth[right] < -SIMD_EPSILON )printf("negative displaced volume 2\n");
			if( columns.m_fluidDepth[down] < -SIMD_EPSILON )printf("negative displaced volume 3\n");
			if( columns.m_fluidDepth[up] < -SIMD_EPSILON )printf("negative displaced volume 4\n");
			
			columns.m_fluidDepth[left] = btMax( columns.m_fluidDepth[left], btScalar(0.0) );
			columns.m_fluidDepth[right] = btMax( columns.m_fluidDepth[right], btScalar(0.0) );
			columns.m_fluidDepth[down] = btMax( columns.m_fluidDepth[down], btScalar(0.0) );
			columns.m_fluidDepth[up] = btMax( columns.m_fluidDepth[up], btScalar(0.0) );
			
			//OPTIMIZATION: zero out next frames r value
			columns.m_displaced[nextFrameDisplacedIndex][center] = btScalar(0.0);
		}
	}
	
	columns.m_displacedIndex = nextFrameDisplacedIndex; // flip frame
}
