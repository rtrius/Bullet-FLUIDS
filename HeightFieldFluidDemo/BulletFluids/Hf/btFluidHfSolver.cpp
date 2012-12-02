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

#include "btFluidHfSolver.h"

#include "LinearMath/btMinMax.h"

btScalar btFluidHfSolverDefault::bilinearInterpolate(const btFluidColumns& columns, const btAlignedObjectArray<btScalar>& array, btScalar iPos, btScalar jPos)
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

	btScalar a = jParam0 * SW + jParam1 * NW;
	btScalar b = jParam0 * SE + jParam1 * NE;
	return iParam0 * a + iParam1 * b;
}

btScalar btFluidHfSolverDefault::advect(const btFluidColumns& columns, const btAlignedObjectArray<btScalar>& array, 
						btScalar i, btScalar j, btScalar vel_x, btScalar vel_z,btScalar dt)
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

	return bilinearInterpolate(columns, array, srcI, srcJ);
}

void btFluidHfSolverDefault::advectEta(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex (i, j);

			if (!columns.m_active[index])
				continue;

			btScalar vel_x = hfParameters.m_globalVelocityX;
			btScalar vel_z = hfParameters.m_globalVelocityZ;

			vel_x += (columns.m_vel_x[index]+columns.m_vel_x[index+1]) * btScalar(0.5);
			vel_z += (columns.m_vel_z[index]+columns.m_vel_z[index+columns.m_numNodesX]) * btScalar(0.5);

			columns.m_temp[index] = advect(columns, columns.m_fluidDepth, btScalar(i), btScalar(j), vel_x, vel_z, dt);
		}
	}
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex (i, j);
			columns.m_fluidDepth[index] = columns.m_temp[index];
		}
	}
}
void btFluidHfSolverDefault::advectVelocityU(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex (i, j);
			if (!columns.m_active[index])
			{
				continue;
			}
			btScalar vel_x = hfParameters.m_globalVelocityX;
			btScalar vel_z = hfParameters.m_globalVelocityZ;

			vel_x += columns.m_vel_x[index];
			vel_z += (columns.m_vel_z[index]+columns.m_vel_z[index+1]+columns.m_vel_z[index+columns.m_numNodesX]+columns.m_vel_z[index+columns.m_numNodesX+1]) * btScalar(0.25);

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
void btFluidHfSolverDefault::advectVelocityV(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex (i, j);
			if (!columns.m_active[index])
			{
				continue;
			}
			btScalar u = hfParameters.m_globalVelocityX;
			btScalar v = hfParameters.m_globalVelocityZ;

			u += (columns.m_vel_x[index]+columns.m_vel_x[index+1]+columns.m_vel_x[index+columns.m_numNodesX]+columns.m_vel_x[index+columns.m_numNodesX+1]) * btScalar(0.25);
			v += columns.m_vel_z[index];
			

			columns.m_temp[index] = advect(columns, columns.m_vel_z, btScalar(i), btScalar(j), u, v, dt);
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

void btFluidHfSolverDefault::updateHeight(const btScalar dt, btFluidColumns& columns)
{
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
			btScalar deta = -btScalar(0.5f) * columns.m_fluidDepth[index] * dt * columns.m_gridCellWidthInv * ( (columns.m_vel_x[index+1] - columns.m_vel_x[index]) + (columns.m_vel_z[index+columns.m_numNodesX] - columns.m_vel_z[index]));
			columns.m_fluidDepth[index] += deta;
			columns.m_combinedHeight[index] = columns.m_ground[index] + btMax(columns.m_fluidDepth[index],btScalar(0.0f));
		}
	}
}
void btFluidHfSolverDefault::updateVelocity(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 2; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index])
			{
				continue;
			}
			columns.m_vel_x[index] += hfParameters.m_gravity * dt * columns.m_gridCellWidthInv * (columns.m_combinedHeight[index]-columns.m_combinedHeight[index-1]);
		}
	}

	for (int j = 2; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_active[index])
			{
				continue;
			}
			columns.m_vel_z[index] += hfParameters.m_gravity * dt * columns.m_gridCellWidthInv * (columns.m_combinedHeight[index]-columns.m_combinedHeight[index-columns.m_numNodesX]);
		}
	}
}

void btFluidHfSolverDefault::computeFlagsAndFillRatio(const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			btScalar hMin = computeHmin(columns, i, j);
			btScalar hMax = computeHmax(hfParameters.m_heightEpsilon, columns, i, j);
			btScalar etaMax = computeEtaMax(columns, i, j);
			
			btScalar h = columns.m_combinedHeight[columns.getIndex(i,j)];
			if (h <= hMin && etaMax < hfParameters.m_dryCellEpsilon)
			{
				columns.m_active[columns.getIndex(i,j)] = false;
				columns.m_fillRatio[columns.getIndex(i,j)] = btScalar(0.0f);
			} 
			else if (h > hMax)
			{
				columns.m_active[columns.getIndex(i,j)] = true;
				columns.m_fillRatio[columns.getIndex(i,j)] = btScalar(1.0f);
			} 
			else
			{
				columns.m_active[columns.getIndex(i,j)] = true;
				columns.m_fillRatio[columns.getIndex(i,j)] = (h - hMin) / (hMax - hMin);
			}
			
		}
	}
}

btScalar btFluidHfSolverDefault::computeHmin(btFluidColumns& columns, int i, int j)
{
	btScalar h1 = columns.m_combinedHeight[columns.getIndex(i-1,j-1)];
	btScalar h2 = columns.m_combinedHeight[columns.getIndex(i-1,j+1)];
	btScalar h3 = columns.m_combinedHeight[columns.getIndex(i+1,j-1)];
	btScalar h4 = columns.m_combinedHeight[columns.getIndex(i+1,j+1)];
	btScalar minh = btMin( btMin(h1, h2), btMin(h3, h4) );
	
	btScalar h = columns.m_combinedHeight[columns.getIndex(i,j)];
	return (minh + h) * btScalar(0.5f);
}

btScalar btFluidHfSolverDefault::computeHmax(const btScalar epsHeight, btFluidColumns& columns, int i, int j)
{
	btScalar h1 = columns.m_combinedHeight[columns.getIndex(i-1,j-1)];
	btScalar h2 = columns.m_combinedHeight[columns.getIndex(i-1,j+1)];
	btScalar h3 = columns.m_combinedHeight[columns.getIndex(i+1,j-1)];
	btScalar h4 = columns.m_combinedHeight[columns.getIndex(i+1,j+1)];
	btScalar maxh = btMax( btMax(h1, h2), btMax(h3, h4) );
	
	btScalar h = columns.m_combinedHeight[columns.getIndex(i,j)];
	return (maxh + h) * btScalar(0.5f) + epsHeight;
}

btScalar btFluidHfSolverDefault::computeEtaMax(btFluidColumns& columns, int i, int j)
{
	btScalar eta1 = columns.m_fluidDepth[columns.getIndex(i-1,j-1)];
	btScalar eta2 = columns.m_fluidDepth[columns.getIndex(i-1,j+1)];
	btScalar eta3 = columns.m_fluidDepth[columns.getIndex(i+1,j-1)];
	btScalar eta4 = columns.m_fluidDepth[columns.getIndex(i+1,j+1)];
	btScalar maxeta = btMax( btMax(eta1, eta2), btMax(eta3,eta4) );

	btScalar eta = columns.m_fluidDepth[columns.getIndex(i,j)];
	return (maxeta + eta) * btScalar(0.5f);
}
void btFluidHfSolverDefault::transferDisplaced(const btScalar volumeDisplacementScale, btFluidColumns& columns)
{
	int currentDisplacedIndex = columns.m_displacedIndex;
	int nextFrameDisplacedIndex = (columns.m_displacedIndex + 1) % 2;

	for (int i = 2; i < columns.m_numNodesX - 2; i++)
	{
		for (int j = 2; j < columns.m_numNodesZ - 2; j++)
		{
			int centerIndex = columns.getIndex(i, j);
		
			btScalar deltaR = columns.m_displaced[currentDisplacedIndex][centerIndex] - columns.m_displaced[nextFrameDisplacedIndex][centerIndex];
			deltaR /= columns.m_gridCellWidth * columns.m_gridCellWidth;	//deltaR is in volume, but we want to change the height
			deltaR *= volumeDisplacementScale;
			
			btScalar qdeltaR = deltaR / btScalar(4.0f);
			columns.m_fluidDepth[columns.getIndex(i-1,j-1)] += qdeltaR;
			columns.m_fluidDepth[columns.getIndex(i-1,j+1)] += qdeltaR;
			columns.m_fluidDepth[columns.getIndex(i+1,j-1)] += qdeltaR;
			columns.m_fluidDepth[columns.getIndex(i+1,j+1)] += qdeltaR;
			columns.m_fluidDepth[centerIndex] -= deltaR;
			
			// OPTIMIZATION: zero out next frames r value
			columns.m_displaced[nextFrameDisplacedIndex][centerIndex] = btScalar(0.0);
		}
	}
	
	columns.m_displacedIndex = nextFrameDisplacedIndex; // flip frame
}

void btFluidHfSolverDefault::setReflectBoundaryLeft(btFluidColumns& columns)
{
	for (int j = 0; j < columns.m_numNodesZ; j++)
	{
		int indexL = columns.getIndex(0, j);

		columns.m_combinedHeight[indexL] = columns.m_combinedHeight[indexL+1];
		columns.m_vel_x[indexL+1] = btScalar(0.0);
		columns.m_vel_z[indexL] = btScalar(0.0);
	}
}

void btFluidHfSolverDefault::setReflectBoundaryRight(btFluidColumns& columns)
{
	for (int j = 0; j < columns.m_numNodesZ; j++)
	{
		int indexR = columns.getIndex(columns.m_numNodesX-1, j);

		columns.m_combinedHeight[indexR] = columns.m_combinedHeight[indexR-1];
		columns.m_vel_x[indexR-1] = btScalar(0.0);
		columns.m_vel_z[indexR] = btScalar(0.0);
	}
}

void btFluidHfSolverDefault::setReflectBoundaryBottom(btFluidColumns& columns)
{
	for (int i = 0; i < columns.m_numNodesX; i++)
	{
		int indexT = columns.getIndex(i, 0);
		
		columns.m_combinedHeight[indexT] = columns.m_combinedHeight[indexT+columns.m_numNodesX];
		columns.m_vel_z[indexT+columns.m_numNodesX] = btScalar(0.0);
		columns.m_vel_x[indexT] = btScalar(0.0);
	}
}

void btFluidHfSolverDefault::setReflectBoundaryTop(btFluidColumns& columns)
{
	for (int i = 0; i < columns.m_numNodesX; i++)
	{
		int indexB = columns.getIndex(i, columns.m_numNodesZ-1);

		columns.m_combinedHeight[indexB] = columns.m_combinedHeight[indexB-columns.m_numNodesX];
		columns.m_vel_z[indexB-columns.m_numNodesX] = btScalar(0.0);
		columns.m_vel_x[indexB] = btScalar(0.0);
	}
}