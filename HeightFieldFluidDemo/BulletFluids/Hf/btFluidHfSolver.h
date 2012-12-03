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

#ifndef BT_FLUID_HF_SOLVER_H
#define BT_FLUID_HF_SOLVER_H


#include "btFluidColumns.h"

///Determines how the columns of a btFluidHf move.
class btFluidHfSolver
{
public:
	virtual void stepSimulation(btScalar timeStep, const btFluidHfParameters& hfParameters, btFluidColumns& columns) = 0;
};


#include <cstdio>	//	remove
inline btScalar getTotalFluidHeight(const btFluidColumns& columns)
{
	btScalar fluidHeight = btScalar(0.0);
	
	for (int j = 1; j < columns.m_numNodesZ-1; j++)
		for (int i = 1; i < columns.m_numNodesX-1; i++)
		{
			int index = columns.getIndex(i, j);
			fluidHeight += columns.m_fluidDepth[index];
		}
		
	return fluidHeight;
}

class btFluidHfSolverDefault : public btFluidHfSolver
{
public:
	virtual void stepSimulation(btScalar timeStep, const btFluidHfParameters& hfParameters, btFluidColumns& columns)
	{
		btScalar initialHeight = getTotalFluidHeight(columns);
	
	
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
		
		
		btScalar finalHeight = getTotalFluidHeight(columns);
		printf("finalHeight - initialHeight: %f \n", finalHeight - initialHeight);
		printf("finalHeight: %f \n", finalHeight);
		printf("\n");
	}
	
	static void computeFlagsAndFillRatio(const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	
protected:
	static btScalar bilinearInterpolate(const btFluidColumns& columns, const btAlignedObjectArray<btScalar>& array, btScalar iPos, btScalar jPos);
	static btScalar advect(const btFluidColumns& columns, const btAlignedObjectArray<btScalar>& array, 
							btScalar i, btScalar j, btScalar di, btScalar dj,btScalar dt);
	static void advectEta(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	static void advectVelocityU(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	static void advectVelocityV(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	
	static void updateHeight(const btScalar dt, btFluidColumns& columns);
	static void updateVelocity(const btScalar dt, const btFluidHfParameters& hfParameters, btFluidColumns& columns);
	
	static btScalar computeHmin(btFluidColumns& columns, int i, int j);
	static btScalar computeHmax(const btScalar epsHeight, btFluidColumns& columns, int i, int j);
	static btScalar computeEtaMax(btFluidColumns& columns, int i, int j);
	static void transferDisplaced(const btScalar volumeDisplacementScale, btFluidColumns& columns);

	static void setReflectBoundaryLeft(btFluidColumns& columns);
	static void setReflectBoundaryRight(btFluidColumns& columns);
	static void setReflectBoundaryBottom(btFluidColumns& columns);
	static void setReflectBoundaryTop(btFluidColumns& columns);
};


#endif
