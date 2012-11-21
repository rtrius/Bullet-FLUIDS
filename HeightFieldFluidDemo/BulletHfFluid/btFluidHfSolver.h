#ifndef BT_FLUID_HF_SOLVER_H
#define BT_FLUID_HF_SOLVER_H

#include "btFluidColumns.h"

class btFluidHfSolver
{
public:
	virtual void stepSimulation(btScalar timeStep, const btFluidHfParameters& hfParameters, btFluidColumns& columns) = 0;
};


class btFluidHfSolverDefault : public btFluidHfSolver
{
public:
	virtual void stepSimulation(btScalar timeStep, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
	{
		transferDisplaced(hfParameters.m_volumeDisplacementScale, columns);
		
		advectEta(timeStep, hfParameters, columns);
		advectVelocityU(timeStep, hfParameters, columns);
		advectVelocityV(timeStep, hfParameters, columns);
		
		updateHeight(timeStep, columns);
		computeFlagsAndFillRatio(hfParameters, columns);
		updateVelocity(timeStep, hfParameters, columns);	
		
		setReflectBoundaryLeft(columns);
		setReflectBoundaryRight(columns);
		setReflectBoundaryTop(columns);
		setReflectBoundaryBottom(columns);
	}
	
protected:
	static btScalar bilinearInterpolate(const btFluidColumns& columns, const btAlignedObjectArray<btScalar>& array, btScalar iPos, btScalar jPos);
	static btScalar advect(const btFluidColumns& columns, const btAlignedObjectArray<btScalar>& array, 
							btScalar i, btScalar j, btScalar di, btScalar dj,btScalar dt);
	static void advectEta(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	static void advectVelocityU(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	static void advectVelocityV(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	
	static void updateHeight(const btScalar dt, btFluidColumns& columns);
	static void updateVelocity(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	
public:
	static void computeFlagsAndFillRatio(const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	
protected:
	static btScalar computeHmin(btFluidColumns& columns, int i, int j);
	static btScalar computeHmax(const btScalar epsHeight, btFluidColumns& columns, int i, int j);
	static btScalar computeEtaMax(btFluidColumns& columns, int i, int j);
	static void transferDisplaced(const btScalar volumeDisplacementScale, btFluidColumns& columns);

	static void setReflectBoundaryLeft(btFluidColumns& columns);
	static void setReflectBoundaryRight(btFluidColumns& columns);
	static void setReflectBoundaryBottom(btFluidColumns& columns);
	static void setReflectBoundaryTop(btFluidColumns& columns);
};


class btFluidHfSolverExperimental : public btFluidHfSolverDefault
{
public:
	virtual void stepSimulation(btScalar timeStep, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
	{
		//transferDisplaced(hfParameters.m_volumeDisplacementScale, columns);
		
		advectEta(timeStep, hfParameters, columns);
		advectVelocityU(timeStep, hfParameters, columns);
		advectVelocityV(timeStep, hfParameters, columns);
		
		updateHeight(timeStep, columns);
		
		if(0)
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
		
		//Set additional reflect cells
		if(0)
		{
			const btScalar ACTIVE_CELL_EPSILON = btScalar(0.0001) * columns.m_gridCellWidth;
		
			for (int j = 1; j < columns.m_numNodesZ-1; j++)
				for (int i = 1; i < columns.m_numNodesX-1; i++)
				{
					int index = columns.getIndex(i, j);
					
					if( !columns.m_active[index] )
					{
						columns.m_vel_x[index] = btScalar(0.0);
						columns.m_vel_z[index] = btScalar(0.0);
					}
					
					if( (columns.m_fluidDepth[index] <= ACTIVE_CELL_EPSILON && columns.m_ground[index] > columns.m_combinedHeight[index+1])
					 || (columns.m_fluidDepth[index+1] <= ACTIVE_CELL_EPSILON && columns.m_ground[index+1] > columns.m_combinedHeight[index]) )
					{
						//columns.m_vel_x[index] = btScalar(0.0);
					}
					
					if( (columns.m_fluidDepth[index] <= ACTIVE_CELL_EPSILON && columns.m_ground[index] > columns.m_combinedHeight[index+columns.m_numNodesX])
					 || (columns.m_fluidDepth[index+columns.m_numNodesX] <= ACTIVE_CELL_EPSILON && columns.m_ground[index+columns.m_numNodesX] > columns.m_combinedHeight[index]) )
					{
						//columns.m_vel_z[index] = btScalar(0.0);
					}
				}
		}
		
		setReflectBoundaryLeft(columns);
		setReflectBoundaryRight(columns);
		setReflectBoundaryTop(columns);
		setReflectBoundaryBottom(columns);
	}
};


#endif