/*
Bullet-FLUIDS 
Copyright (c) 2012-2013 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_FLUID_SPH_SOLVER_PBF_H
#define BT_FLUID_SPH_SOLVER_PBF_H

#include "BulletFluids/Sph/btFluidSphSolver.h"

///Work in progress; currently not functional(do not use)
///Rigid body interaction not implemented
///@remarks
///This solver is based on the method described in: \n
///"Position Based Fluids" \n
///M. Macklin and M. Muller. \n
///ACM Transactions on Graphics(TOG) - Proceedings of ACM SIGGRAPH 2013, v.32 n.4, July 2013. \n
class btFluidSphSolverPBF : public btFluidSphSolver
{
public:
	struct PbfParticles
	{
		btAlignedObjectArray<btFluidSphNeighbors> m_neighborTable;
		
		btAlignedObjectArray<btVector3> m_predictedPosition;
		btAlignedObjectArray<btVector3> m_predictedVelocity;
		
		btAlignedObjectArray<btVector3> m_deltaPosition;
		btAlignedObjectArray<btVector3> m_nextVelocity;
		btAlignedObjectArray<btVector3> m_xsphViscosity;
		
		btAlignedObjectArray<btScalar> m_density;
		btAlignedObjectArray<btScalar> m_scalingFactorDenominator;
		btAlignedObjectArray<btScalar> m_scalingFactor;
		
		int size() const { return m_neighborTable.size(); }
		void resize(int newSize)
		{
			m_neighborTable.resize(newSize);
			
			m_predictedPosition.resize(newSize);
			m_predictedVelocity.resize(newSize);
			
			m_deltaPosition.resize(newSize);
			m_nextVelocity.resize(newSize);
			m_xsphViscosity.resize(newSize);
			
			m_density.resize(newSize);
			m_scalingFactorDenominator.resize(newSize);
			m_scalingFactor.resize(newSize);
		}
	};

protected:
	btAlignedObjectArray<btFluidSphSolverPBF::PbfParticles> m_pbfData;

public:
	virtual void updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids);
	
	virtual bool isPositionBasedSolver() const { return true; }
	
	static void findNeighborsInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
											const btFluidSortingGrid& grid, btFluidParticles& particles,
											btFluidSphSolverPBF::PbfParticles& m_pbfData);
};

#endif


