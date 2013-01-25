/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_FLUID_SPH_SOLVER_PCISPH_H
#define BT_FLUID_SPH_SOLVER_PCISPH_H

#include "BulletFluids/Sph/btFluidSphSolver.h"

///Do not use; work in progress solver with improved incompressibility
class btFluidSphSolverPCISPH : public btFluidSphSolver
{
public:
	struct PciSphParticles
	{
		btAlignedObjectArray<btFluidSphNeighbors> m_neighborTable;
		
		btAlignedObjectArray<btVector3> m_viscosityForce;
		btAlignedObjectArray<btVector3> m_pressureForce;
		btAlignedObjectArray<btVector3> m_predictedPosition;
		btAlignedObjectArray<btVector3> m_predictedVelocity;
		
		btAlignedObjectArray<btScalar> m_density;
		btAlignedObjectArray<btScalar> m_densityError;
		btAlignedObjectArray<btScalar> m_pressure;
		
		int size() const { return m_neighborTable.size(); }
		void resize(int newSize)
		{
			m_neighborTable.resize(newSize);
			
			m_viscosityForce.resize(newSize);
			m_pressureForce.resize(newSize);
			m_predictedPosition.resize(newSize);
			m_predictedVelocity.resize(newSize);
			
			m_density.resize(newSize);
			m_densityError.resize(newSize);
			m_pressure.resize(newSize);
		}
	};

protected:
	btAlignedObjectArray<btFluidSphSolverPCISPH::PciSphParticles> m_pciSphData;

public:
	virtual void updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids);
};

#endif


