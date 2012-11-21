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
						btScalar i, btScalar j, btScalar di, btScalar dj,btScalar dt)
{
	// trace particle backwards in time
	btScalar srcI = i - di * dt * columns.m_gridCellWidthInv;
	btScalar srcJ = j - dj * dt * columns.m_gridCellWidthInv;

	// srcI and srcJ are indices into the array,
	// we need to clamp them to within the domain
	srcI = btMax (srcI, btScalar(0.0));
	srcI = btMin (srcI, btScalar(columns.m_numNodesWidth-1));
	srcJ = btMax (srcJ, btScalar(0.0));
	srcJ = btMin (srcJ, btScalar(columns.m_numNodesLength-1));

	return bilinearInterpolate(columns, array, srcI, srcJ);
}

void btFluidHfSolverDefault::advectEta(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesLength-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesWidth-1; i++)
		{
			int index = columns.getIndex (i, j);

			if (!columns.m_flags[index])
				continue;

			btScalar u = hfParameters.m_globalVelocityU;
			btScalar v = hfParameters.m_globalVelocityV;

			u += (columns.m_u[index]+columns.m_u[index+1]) * btScalar(0.5);
			v += (columns.m_v[index]+columns.m_v[index+columns.m_numNodesWidth]) * btScalar(0.5);

			columns.m_temp[index] = advect(columns, columns.m_eta, btScalar(i), btScalar(j), u, v, dt);
		}
	}
	for (int j = 1; j < columns.m_numNodesLength-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesWidth-1; i++)
		{
			int index = columns.getIndex (i, j);
			columns.m_eta[index] = columns.m_temp[index];
		}
	}
}

void btFluidHfSolverDefault::advectVelocityU(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesLength-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesWidth-1; i++)
		{
			int index = columns.getIndex (i, j);
			if (!columns.m_flags[index])
			{
				continue;
			}
			btScalar u = hfParameters.m_globalVelocityU;
			btScalar v = hfParameters.m_globalVelocityV;

			u += columns.m_u[index];
			v += (columns.m_v[index]+columns.m_v[index+1]+columns.m_v[index+columns.m_numNodesWidth]+columns.m_v[index+columns.m_numNodesWidth+1]) * btScalar(0.25);

			columns.m_temp[index] = advect(columns, columns.m_u, btScalar(i), btScalar(j), u, v, dt);
		}
	}

	for (int j = 1; j < columns.m_numNodesLength-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesWidth-1; i++)
		{
			int index = columns.getIndex (i, j);
			columns.m_u[index] = columns.m_temp[index];
		}
	}
}

void btFluidHfSolverDefault::advectVelocityV(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesLength-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesWidth-1; i++)
		{
			int index = columns.getIndex (i, j);
			if (!columns.m_flags[index])
			{
				continue;
			}
			btScalar u = hfParameters.m_globalVelocityU;
			btScalar v = hfParameters.m_globalVelocityV;

			u += (columns.m_u[index]+columns.m_u[index+1]+columns.m_u[index+columns.m_numNodesWidth]+columns.m_u[index+columns.m_numNodesWidth+1]) * btScalar(0.25);
			v += columns.m_v[index];
			

			columns.m_temp[index] = advect(columns, columns.m_v, btScalar(i), btScalar(j), u, v, dt);
		}
	}
	
	for (int j = 1; j < columns.m_numNodesLength-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesWidth-1; i++)
		{
			int index = columns.getIndex (i, j);
			columns.m_v[index] = columns.m_temp[index];
		}
	}
}

void btFluidHfSolverDefault::updateHeight(const btScalar dt, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesLength-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesWidth-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_flags[index])
			{
				columns.m_height[index] = columns.m_ground[index] + columns.m_eta[index];
				continue;
			}
			btScalar deta = -btScalar(0.5f) * columns.m_eta[index] * dt * columns.m_gridCellWidthInv * ( (columns.m_u[index+1] - columns.m_u[index]) + (columns.m_v[index+columns.m_numNodesWidth] - columns.m_v[index]));
			columns.m_eta[index] += deta;
			columns.m_height[index] = columns.m_ground[index] + btMax(columns.m_eta[index],btScalar(0.0f));
		}
	}
}
void btFluidHfSolverDefault::updateVelocity(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesLength-1; j++)
	{
		for (int i = 2; i < columns.m_numNodesWidth-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_flags[index])
			{
				continue;
			}
			columns.m_u[index] += hfParameters.m_gravity * dt * columns.m_gridCellWidthInv * (columns.m_height[index]-columns.m_height[index-1]);
		}
	}

	for (int j = 2; j < columns.m_numNodesLength-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesWidth-1; i++)
		{
			int index = columns.getIndex(i, j);
			if (!columns.m_flags[index])
			{
				continue;
			}
			columns.m_v[index] += hfParameters.m_gravity * dt * columns.m_gridCellWidthInv * (columns.m_height[index]-columns.m_height[index-columns.m_numNodesWidth]);
		}
	}
}

void btFluidHfSolverDefault::computeFlagsAndFillRatio(const btFluidHfParameters& hfParameters, btFluidColumns& columns)
{
	for (int j = 1; j < columns.m_numNodesLength-1; j++)
	{
		for (int i = 1; i < columns.m_numNodesWidth-1; i++)
		{
			btScalar hMin = computeHmin(columns, i, j);
			btScalar hMax = computeHmax(hfParameters.m_epsHeight, columns, i, j);
			btScalar etaMax = computeEtaMax(columns, i, j);
			
			btScalar h = columns.m_height[columns.getIndex(i,j)];
			if (h <= hMin && etaMax < hfParameters.m_epsEta)
			{
				columns.m_flags[columns.getIndex(i,j)] = false;
				columns.m_fillRatio[columns.getIndex(i,j)] = btScalar(0.0f);
			} 
			else if (h > hMax)
			{
				columns.m_flags[columns.getIndex(i,j)] = true;
				columns.m_fillRatio[columns.getIndex(i,j)] = btScalar(1.0f);
			} 
			else
			{
				columns.m_flags[columns.getIndex(i,j)] = true;
				columns.m_fillRatio[columns.getIndex(i,j)] = (h - hMin) / (hMax - hMin);
			}
			
		}
	}
}

btScalar btFluidHfSolverDefault::computeHmin(btFluidColumns& columns, int i, int j)
{
	btScalar h1 = columns.m_height[columns.getIndex(i-1,j-1)];
	btScalar h2 = columns.m_height[columns.getIndex(i-1,j+1)];
	btScalar h3 = columns.m_height[columns.getIndex(i+1,j-1)];
	btScalar h4 = columns.m_height[columns.getIndex(i+1,j+1)];
	btScalar minh = btMin( btMin(h1, h2), btMin(h3, h4) );
	
	btScalar h = columns.m_height[columns.getIndex(i,j)];
	return (minh + h) * btScalar(0.5f);
}

btScalar btFluidHfSolverDefault::computeHmax(const btScalar epsHeight, btFluidColumns& columns, int i, int j)
{
	btScalar h1 = columns.m_height[columns.getIndex(i-1,j-1)];
	btScalar h2 = columns.m_height[columns.getIndex(i-1,j+1)];
	btScalar h3 = columns.m_height[columns.getIndex(i+1,j-1)];
	btScalar h4 = columns.m_height[columns.getIndex(i+1,j+1)];
	btScalar maxh = btMax( btMax(h1, h2), btMax(h3, h4) );
	
	btScalar h = columns.m_height[columns.getIndex(i,j)];
	return (maxh + h) * btScalar(0.5f) + epsHeight;
}

btScalar btFluidHfSolverDefault::computeEtaMax(btFluidColumns& columns, int i, int j)
{
	btScalar eta1 = columns.m_eta[columns.getIndex(i-1,j-1)];
	btScalar eta2 = columns.m_eta[columns.getIndex(i-1,j+1)];
	btScalar eta3 = columns.m_eta[columns.getIndex(i+1,j-1)];
	btScalar eta4 = columns.m_eta[columns.getIndex(i+1,j+1)];
	btScalar maxeta = btMax( btMax(eta1, eta2), btMax(eta3,eta4) );

	btScalar eta = columns.m_eta[columns.getIndex(i,j)];
	return (maxeta + eta) * btScalar(0.5f);
}

void btFluidHfSolverDefault::transferDisplaced(const btScalar volumeDisplacementScale, btFluidColumns& columns)
{
	for (int i = 2; i < columns.m_numNodesWidth - 2; i++)
	{
		for (int j = 2; j < columns.m_numNodesLength - 2; j++)
		{
			btScalar deltaR = columns.m_r[columns.m_rIndex][columns.getIndex(i,j)] - columns.m_r[(columns.m_rIndex+1)%2][columns.getIndex(i,j)];
			deltaR /= columns.m_gridCellWidth * columns.m_gridCellWidth;	//deltaR is in volume, but we want to change the height
			deltaR *= volumeDisplacementScale;
			
			btScalar qdeltaR = deltaR / btScalar(4.0f);
			columns.m_eta[columns.getIndex(i-1,j-1)] += qdeltaR;
			columns.m_eta[columns.getIndex(i-1,j+1)] += qdeltaR;
			columns.m_eta[columns.getIndex(i+1,j-1)] += qdeltaR;
			columns.m_eta[columns.getIndex(i+1,j+1)] += qdeltaR;
			columns.m_eta[columns.getIndex(i,j)] -= deltaR;
			
			// OPTIMIZATION: zero out next frames r value
			columns.m_r[(columns.m_rIndex+1)%2][columns.getIndex(i,j)] = btScalar(0.0);
		}
	}
	
	columns.m_rIndex = (columns.m_rIndex + 1) % 2; // flip frame
}

void btFluidHfSolverDefault::setReflectBoundaryLeft(btFluidColumns& columns)
{
	for (int j = 0; j < columns.m_numNodesLength; j++)
	{
		int indexL = columns.getIndex(0, j);

		columns.m_height[indexL] = columns.m_height[indexL+1];
		columns.m_u[indexL+1] = btScalar(0.0);
		columns.m_v[indexL] = btScalar(0.0);
	}
}

void btFluidHfSolverDefault::setReflectBoundaryRight(btFluidColumns& columns)
{
	for (int j = 0; j < columns.m_numNodesLength; j++)
	{
		int indexR = columns.getIndex(columns.m_numNodesWidth-1, j);

		columns.m_height[indexR] = columns.m_height[indexR-1];
		columns.m_u[indexR-1] = btScalar(0.0);
		columns.m_v[indexR] = btScalar(0.0);
	}
}

void btFluidHfSolverDefault::setReflectBoundaryBottom(btFluidColumns& columns)
{
	for (int i = 0; i < columns.m_numNodesWidth; i++)
	{
		int indexT = columns.getIndex(i, 0);
		
		columns.m_height[indexT] = columns.m_height[indexT+columns.m_numNodesWidth];
		columns.m_v[indexT+columns.m_numNodesWidth] = btScalar(0.0);
		columns.m_u[indexT] = btScalar(0.0);
	}
}

void btFluidHfSolverDefault::setReflectBoundaryTop(btFluidColumns& columns)
{
	for (int i = 0; i < columns.m_numNodesWidth; i++)
	{
		int indexB = columns.getIndex(i, columns.m_numNodesLength-1);

		columns.m_height[indexB] = columns.m_height[indexB-columns.m_numNodesWidth];
		columns.m_v[indexB-columns.m_numNodesWidth] = btScalar(0.0);
		columns.m_u[indexB] = btScalar(0.0);
	}
}