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
		updateHeight_new(timeStep, columns);
		
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
	
	static void updateHeight_new(const btScalar dt, btFluidColumns& columns);
	
	static void applyBoundaryConditions(btFluidColumns& columns);
};


#endif