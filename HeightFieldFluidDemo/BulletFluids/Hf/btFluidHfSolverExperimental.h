#ifndef BT_FLUID_HF_SOLVER_EXPERIMENTAL_H
#define BT_FLUID_HF_SOLVER_EXPERIMENTAL_H

#include "btFluidColumns.h"
#include "btFluidHfSolver.h"

class btFluidHfSolverExperimental : public btFluidHfSolverDefault
{
public:
	virtual void stepSimulation(btScalar timeStep, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
	{
		btScalar initialHeight = getTotalFluidHeight(columns);
	
		transferDisplaced(hfParameters.m_volumeDisplacementScale, columns);
		
		advectEta(timeStep, hfParameters, columns);
		{
			btScalar heightAfterAdvection = getTotalFluidHeight(columns);
			printf("heightAfterAdvection - initialHeight: %f \n", heightAfterAdvection - initialHeight);
		}
		
		advectVelocityX_New(timeStep, hfParameters, columns);
		advectVelocityZ_New(timeStep, hfParameters, columns);
		
		updateHeight(timeStep, columns);
		
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
		
		btScalar finalHeight = getTotalFluidHeight(columns);
		printf("finalHeight - initialHeight: %f \n", finalHeight - initialHeight);
		printf("finalHeight: %f \n", finalHeight);
		printf("\n");
	}
	
protected:
	static void advectVelocityX_New(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	static void advectVelocityZ_New(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	
	static void applyBoundaryConditions(btFluidColumns& columns);
};


#endif