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
#ifndef BT_FLUID_SPH_SOLVER_MULTIPHASE_H
#define BT_FLUID_SPH_SOLVER_MULTIPHASE_H

#include "btFluidSphSolver.h"


///@brief Experimental, unoptimized solver supporting fluid-fluid interaction.
///@remarks
///This solver has issues when btFluidSph with differing btFluidSphParametersLocal interact:
/// - Fluid particles will stick to boundaries; for instance, some particles of lighter fluids
///will not rise even if heavier fluid particles are on top.
class btFluidSphSolverMultiphase : public btFluidSphSolverDefault
{
public:
	virtual void updateGridAndCalculateSphForces(const btFluidSphParametersGlobal& FG, btFluidSph** fluids, int numFluids);
	
protected:
	virtual void sphComputePressureMultiphase(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, btFluidSphSolverDefault::SphParticles& sphData,
												btAlignedObjectArray<btFluidSph*>& interactingFluids, 
												btAlignedObjectArray<btFluidSphSolverDefault::SphParticles*>& interactingSphData);
	virtual void sphComputeForceMultiphase(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, btFluidSphSolverDefault::SphParticles& sphData,
												btAlignedObjectArray<btFluidSph*>& interactingFluids, 
												btAlignedObjectArray<btFluidSphSolverDefault::SphParticles*>& interactingSphData);
};

#endif