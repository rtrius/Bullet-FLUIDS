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
#ifndef BT_FLUID_HF_SOLVER_EXPERIMENTAL_H
#define BT_FLUID_HF_SOLVER_EXPERIMENTAL_H

#include "btFluidColumns.h"
#include "btFluidHfSolver.h"

class btFluidHfSolverExperimental : public btFluidHfSolverDefault
{
	//For MacCormack advection
	btAlignedObjectArray<btScalar> m_tempForwardAdvect;
	btAlignedObjectArray<btScalar> m_tempForwardAdvectMin;
	btAlignedObjectArray<btScalar> m_tempForwardAdvectMax;
	btAlignedObjectArray<btScalar> m_tempReverseAdvect;
	
public:
	virtual void stepSimulation(btScalar timeStep, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
	{
		int numColumns = columns.numColumns();
		if( m_tempForwardAdvect.size() < numColumns )
		{
			m_tempForwardAdvect.resize(numColumns);
			m_tempForwardAdvectMin.resize(numColumns);
			m_tempForwardAdvectMax.resize(numColumns);
			m_tempReverseAdvect.resize(numColumns);
		}
		
	
		transferDisplaced(hfParameters.m_volumeDisplacementScale, columns);
		
		btScalar initialHeight = getTotalFluidHeight(columns);
		
		const bool USE_MACCORMACK_METHOD = false;
		if(USE_MACCORMACK_METHOD)
		{
			//advectFluid_macCormack(timeStep, hfParameters, columns);
			advectVelocityX_macCormack(timeStep, hfParameters, columns);
			advectVelocityZ_macCormack(timeStep, hfParameters, columns);
		}
		else
		{
			//advectEta(timeStep, hfParameters, columns);
			advectVelocityX_new(timeStep, hfParameters, columns);
			advectVelocityZ_new(timeStep, hfParameters, columns);
		}
		
		//Most volume loss occurs during advection
		{
			btScalar heightAfterAdvection = getTotalFluidHeight(columns);
			printf("heightAfterAdvection - initialHeight: %f \n", heightAfterAdvection - initialHeight);
			printf("initialHeight: %f \n", initialHeight);
			printf("heightAfterAdvection: %f \n", heightAfterAdvection);
		}
		
		//applyBoundaryConditions(columns);
		//updateHeight(timeStep, columns);
		updateHeight_new(timeStep, hfParameters, columns);
		
		if(1)
		{
			for (int j = 1; j < columns.m_numNodesZ-1; j++)
				for (int i = 1; i < columns.m_numNodesX-1; i++)
				{
					int index = columns.getIndex(i, j);
					
					columns.m_active[index] = true;
				}
		}
		else computeFlagsAndFillRatio(hfParameters, columns);
		
		updateVelocity(timeStep, hfParameters, columns);
		
		applyBoundaryConditions(columns);
		
		//Apply fluid flow threshold
		const btScalar ACTIVE_CELL_EPSILON = btScalar(0.0001) * columns.m_gridCellWidth;
		const btScalar minFlowHeight = hfParameters.m_fluidFlowThreshold;
		if( hfParameters.m_fluidFlowThreshold != btScalar(0.0) )
		{
			for (int j = 1; j < columns.m_numNodesZ-1; j++)
				for (int i = 1; i < columns.m_numNodesX-1; i++)
				{
					int index = columns.getIndex(i, j);
					
					if( (columns.m_combinedHeight[index-1] > columns.m_combinedHeight[index] && columns.m_fluidDepth[index-1] < minFlowHeight && columns.m_fluidDepth[index] < ACTIVE_CELL_EPSILON)
					 || (columns.m_combinedHeight[index] > columns.m_combinedHeight[index-1] && columns.m_fluidDepth[index] < minFlowHeight && columns.m_fluidDepth[index-1] < ACTIVE_CELL_EPSILON) )
					{
						columns.m_vel_x[index] = btScalar(0.0);
					}
					
					if( (columns.m_combinedHeight[index-columns.m_numNodesX] > columns.m_combinedHeight[index] && columns.m_fluidDepth[index-columns.m_numNodesX] < minFlowHeight && columns.m_fluidDepth[index] < ACTIVE_CELL_EPSILON)
					 || (columns.m_combinedHeight[index] > columns.m_combinedHeight[index-columns.m_numNodesX] && columns.m_fluidDepth[index] < minFlowHeight && columns.m_fluidDepth[index-columns.m_numNodesX] < ACTIVE_CELL_EPSILON) )
					{
						columns.m_vel_z[index] = btScalar(0.0);
					}
				}
		}
		
		btScalar finalHeight = getTotalFluidHeight(columns);
		printf("finalHeight - initialHeight: %f \n", finalHeight - initialHeight);
		printf("initialHeight: %f \n", initialHeight);
		printf("finalHeight: %f \n", finalHeight);
		printf("\n");
	}
	
protected:

	static btScalar advectReverse(const btFluidColumns& columns, const btAlignedObjectArray<btScalar>& array, 
									btScalar i, btScalar j, btScalar vel_x, btScalar vel_z, btScalar dt);

	void advectFluid_macCormack(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	void advectVelocityX_macCormack(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	void advectVelocityZ_macCormack(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	
	static void advectVelocityX_new(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	static void advectVelocityZ_new(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	
	static void updateHeight_new(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	
	static void applyBoundaryConditions(btFluidColumns& columns);
};


#endif