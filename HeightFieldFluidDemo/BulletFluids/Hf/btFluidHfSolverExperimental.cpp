#include "btFluidHfSolverExperimental.h"

#include "LinearMath/btMinMax.h"


void btFluidHfSolverExperimental::advectVelocityX_New(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
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
void btFluidHfSolverExperimental::advectVelocityZ_New(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
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

void btFluidHfSolverExperimental::applyBoundaryConditions(btFluidColumns& columns)
{
	const btScalar ACTIVE_CELL_EPSILON = btScalar(0.0001) * columns.m_gridCellWidth;

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
	}
	for(int j = 0; j < columns.m_numNodesZ; j++)					//Right(X+)
	{
		int indexR = columns.getIndex(columns.m_numNodesX-1, j);
		columns.m_fluidDepth[indexR] = BOUNDARY_FLUID_DEPTH;
		columns.m_combinedHeight[indexR] = columns.m_ground[indexR];
	}
	for(int i = 0; i < columns.m_numNodesX; i++)
	{
		int indexT = columns.getIndex(i, 0);						//Top(Z-)
		columns.m_fluidDepth[indexT] = BOUNDARY_FLUID_DEPTH;
		columns.m_combinedHeight[indexT] = columns.m_ground[indexT];
	}
	for(int i = 0; i < columns.m_numNodesX; i++)
	{
		int indexB = columns.getIndex(i, columns.m_numNodesZ-1);	//Bottom(Z+)
		columns.m_fluidDepth[indexB] = BOUNDARY_FLUID_DEPTH;
		columns.m_combinedHeight[indexB] = columns.m_ground[indexB];
	}
}


