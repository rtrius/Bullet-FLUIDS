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

#endif