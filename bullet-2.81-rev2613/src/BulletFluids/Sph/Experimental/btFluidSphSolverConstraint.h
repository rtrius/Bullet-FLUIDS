/*
Bullet-FLUIDS
Copyright (c) 2012-2014 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_FLUID_SPH_SOLVER_CONSTRAINT_H
#define BT_FLUID_SPH_SOLVER_CONSTRAINT_H

#include "BulletFluids/Sph/btFluidSphSolver.h"


//Note: The implementation here follows D. Gustafsson's GDC 2014 presentation, "Constraint fluids in Sprinkle".
//It does *not* use the SPOOK integrator as described in the "Constraint Fluids" paper. 
//Instead, the 'Sequential Impulse' method is used. In either case, the main effort is to solve a linear system of the form Ax = b.
//
//The equation to solve is:
//		(1/mass) * (J * J^T) * lambda = (zeta - J * v)
//		(mass/restDensity^2) * (G * G^T) * lambda = zeta - ((mass/restDensity) * G*v)
//where 
//		J = (mass/restDensity) * G						//J is the Jacobian matrix, and G is the Jacobian matrix with some constant terms removed
//		zeta = (density - restDensity) * baumgarte		//Constraint error, C(x1, x2, ..., xN)
//
//		Note that: 
//			-There is no (1/dt) term on the right side as we solve for an impulse and not a force
//			-Since the mass of all particles is constant, the M^-1 matrix is replaced with (1 / mass)
//			-As this is a draft version, the terms accounting for gravity and external forces are omitted
//
//For a single iteration of PGS(1 row):
//		(mass/restDensity^2) * A_ii * lambda_i = zeta_i - ((mass/restDensity) * G_i*v)
//		lambda_i = [zeta_i - ((mass/restDensity) * G_i*v)] / [(mass/restDensity^2) * A_ii]
//
//Discrepancies between "Constraint fluids in Sprinkle" with the "Constraint Fluids" paper:
//	-Uses the Sequential Impulse equation (from "Iterative Dynamics with Temporal Coherence" by E. Catto).
//	-Mass and rest density are always 1, and so do not need to be included in the Jacobian matrix.
//	-Different kernel (1 - (d/h)^2)^3 is used for both SPH density and Jacobian.
//	-Impulses are directly applied (v_next = v_prev + J^T * lambda)
//	-lambda is clamped to the range [0, infinity], which means that only repulsive impulses are allowed.
//	Clamping improves stability, but also changes the problem from a linear system into a LCP.
//
//The assumption that (mass == 1) and (restDensity == 1), results in the simplified linear system:
//		(G * G^T) * lambda = zeta - G*v
//as well as a simplified PGS iteration:
//		lambda_i = [zeta_i - G_i * v] / A_ii		//The (G * G^T)_i * lambda term is missing since the impulses are immediately applied
//

///Sample implementation of the Constraint Fluids method, for incompressible fluid simulation
///@remarks
///This solver is based on: \n
///"Constraint Fluids". \n
///K. Bodin, C. Lacoursiere, M. Servin. \n
class btFluidSphSolverConstraint : public btFluidSphSolver
{
public:
	struct CfParticles
	{
		btAlignedObjectArray<btFluidSphNeighbors> m_neighborTable;
		
		btAlignedObjectArray<btScalar> m_density;
		btAlignedObjectArray<btScalar> m_bias;		///<Scaled constraint error
		btAlignedObjectArray<btScalar> m_A;			///<Contains diagonal elements of the (J * M^-1 * J^T) matrix
		
		int size() const { return m_neighborTable.size(); }
		void resize(int newSize)
		{
			m_neighborTable.resize(newSize);
			
			m_density.resize(newSize);
			m_bias.resize(newSize);
			m_A.resize(newSize);
		}
	};

protected:
	btAlignedObjectArray<btFluidSphSolverConstraint::CfParticles> m_cfData;

public:
	virtual void updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids);
};

#endif


