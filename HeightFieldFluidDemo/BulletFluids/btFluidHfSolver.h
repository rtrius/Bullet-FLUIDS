#ifndef BT_FLUID_HF_SOLVER_H
#define BT_FLUID_HF_SOLVER_H

#include "LinearMath/btMinMax.h"

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
	static btScalar advectReverse(const btFluidColumns& columns, const btAlignedObjectArray<btScalar>& array, 
									btScalar i, btScalar j, btScalar vel_x, btScalar vel_z,btScalar dt)
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

public:
	virtual void stepSimulation(btScalar timeStep, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
	{
		//transferDisplaced(hfParameters.m_volumeDisplacementScale, columns);
		
		advectEta(timeStep, hfParameters, columns);
		advectVelocityU(timeStep, hfParameters, columns);
		advectVelocityV(timeStep, hfParameters, columns);
		
		if(0)
		{
			for (int j = 1; j < columns.m_numNodesZ-1; j++)
			{
				for (int i = 1; i < columns.m_numNodesX-1; i++)
				{
					int index = columns.getIndex(i, j);
					if (!columns.m_active[index])
					{
						columns.m_combinedHeight[index] = columns.m_ground[index] + columns.m_fluidDepth[index];
						columns.m_combinedHeight[index] = columns.m_ground[index] + btMax(columns.m_fluidDepth[index],btScalar(0.0));
						continue;
					}
					
					btScalar u1, u2, v1, v2;
					if(0)
					{
						u1 = columns.m_vel_x[index] * (columns.m_fluidDepth[index] + columns.m_fluidDepth[index+1]) * btScalar(0.5);
						u2 = columns.m_vel_x[index-1] * (columns.m_fluidDepth[index] + columns.m_fluidDepth[index-1]) * btScalar(0.5);
						v1 = columns.m_vel_z[index] * (columns.m_fluidDepth[index] + columns.m_fluidDepth[index+columns.m_numNodesX]) * btScalar(0.5);
						v2 = columns.m_vel_z[index-columns.m_numNodesX] * (columns.m_fluidDepth[index] + columns.m_fluidDepth[index-columns.m_numNodesX]) * btScalar(0.5);
					}
					else
					{
						u1 = columns.m_vel_x[index] * ( columns.m_vel_x[index] > btScalar(0.0) ) ? columns.m_fluidDepth[index] : columns.m_fluidDepth[index+1];
						u2 = columns.m_vel_x[index-1] * ( columns.m_vel_x[index-1] > btScalar(0.0) ) ? columns.m_fluidDepth[index-1] : columns.m_fluidDepth[index];
						v1 = columns.m_vel_z[index] * ( columns.m_vel_z[index] > btScalar(0.0) ) ? columns.m_fluidDepth[index] : columns.m_fluidDepth[index+columns.m_numNodesX];
						v2 = columns.m_vel_z[index-columns.m_numNodesX] * ( columns.m_vel_z[index-columns.m_numNodesX] > btScalar(0.0) ) ? columns.m_fluidDepth[index-columns.m_numNodesX] : columns.m_fluidDepth[index];
					}
					
					btScalar deta = ((u1 - u2) + (v1 - v2)) * - columns.m_gridCellWidthInv;
					columns.m_fluidDepth[index] += deta * timeStep;
					columns.m_fluidDepth[index] = btMax( columns.m_fluidDepth[index], btScalar(0.0) );
					columns.m_combinedHeight[index] = columns.m_ground[index] + btMax(columns.m_fluidDepth[index],btScalar(0.0));
				}
			}
		}
		else updateHeight(timeStep, columns);
		
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
		
		if(0)
		{
			for (int j = 1; j < columns.m_numNodesZ-1; j++)
			{
				for (int i = 1; i < columns.m_numNodesX-2; i++)
				{
					int index = columns.getIndex(i, j);
					if (!columns.m_active[index])
					{
						continue;
					}
					columns.m_vel_x[index] += hfParameters.m_gravity * timeStep * columns.m_gridCellWidthInv * (columns.m_combinedHeight[index+1]-columns.m_combinedHeight[index]);
					columns.m_vel_x[index] = btMin(columns.m_vel_x[index], btScalar(0.5) * columns.m_gridCellWidth / timeStep);
				}
			}

			for (int j = 1; j < columns.m_numNodesZ-2; j++)
			{
				for (int i = 1; i < columns.m_numNodesX-1; i++)
				{
					int index = columns.getIndex(i, j);
					if (!columns.m_active[index])
					{
						continue;
					}
					columns.m_vel_z[index] += hfParameters.m_gravity * timeStep * columns.m_gridCellWidthInv * (columns.m_combinedHeight[index+columns.m_numNodesX]-columns.m_combinedHeight[index]);
					columns.m_vel_z[index] = btMin(columns.m_vel_z[index], btScalar(0.5) * columns.m_gridCellWidth / timeStep);
				}
			}
		}
		else updateVelocity(timeStep, hfParameters, columns);	
		
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