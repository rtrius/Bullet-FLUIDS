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
#ifndef BT_FLUID_SPH_SOLVER_IISPH_H
#define BT_FLUID_SPH_SOLVER_IISPH_H

#include "BulletFluids/Sph/btFluidSphSolver.h"

///Work in progress; currently not functional(do not use)
///@remarks
///This solver implements the method described in: \n
///"Implicit Incompressible SPH"
///M. Ihmsen, J. Cornelis, B. Solenthaler, C. Horvath, and M. Teschner
class btFluidSphSolverIISPH : public btFluidSphSolver
{
public:
	struct IiSphParticles
	{
		btAlignedObjectArray<btFluidSphNeighbors> m_neighborTable;
		
		btAlignedObjectArray<btVector3> m_predictedAcceleration;
		btAlignedObjectArray<btVector3> m_predictedVelocity;
		btAlignedObjectArray<btVector3> m_d_ii;
		btAlignedObjectArray<btVector3> m_d_ij_pj_sum;
		btAlignedObjectArray<btVector3> m_pressureForce;
		
		btAlignedObjectArray<btScalar> m_density;
		btAlignedObjectArray<btScalar> m_density_adv;
		btAlignedObjectArray<btScalar> m_a_ii;
		btAlignedObjectArray<btScalar> m_equation13_sum;
		btAlignedObjectArray<btScalar> m_pressure;
		
		int size() const { return m_neighborTable.size(); }
		void resize(int newSize)
		{
			m_neighborTable.resize(newSize);
			
			m_predictedAcceleration.resize(newSize);
			m_predictedVelocity.resize(newSize);
			m_d_ii.resize(newSize);
			m_d_ij_pj_sum.resize(newSize);
			m_pressureForce.resize(newSize);
			
			m_density.resize(newSize);
			m_density_adv.resize(newSize);
			m_a_ii.resize(newSize);
			m_equation13_sum.resize(newSize);
			m_pressure.resize(newSize);
		}
	};

protected:
	btAlignedObjectArray<btFluidSphSolverIISPH::IiSphParticles> m_iiSphdata;

public:
	virtual void updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids);
	
	static void calculateSumsInCellSymmetric(const btFluidSphParametersGlobal& FG, int gridCellIndex, 
											const btFluidSortingGrid& grid, btFluidParticles& particles,
											btFluidSphSolverIISPH::IiSphParticles& sphData);
};

#endif


