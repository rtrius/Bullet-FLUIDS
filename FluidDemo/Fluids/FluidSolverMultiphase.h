/* FluidSolverMultiphase.h

	ZLib license
	This software is provided 'as-is', without any express or implied
	warranty. In no event will the authors be held liable for any damages
	arising from the use of this software.
	
	Permission is granted to anyone to use this software for any purpose,
	including commercial applications, and to alter it and redistribute it
	freely, subject to the following restrictions:
	
	1. The origin of this software must not be misrepresented; you must not
	   claim that you wrote the original software. If you use this software
	   in a product, an acknowledgment in the product documentation would be
	   appreciated but is not required.
	2. Altered source versions must be plainly marked as such, and must not be
	   misrepresented as being the original software.
	3. This notice may not be removed or altered from any source distribution.
*/
#ifndef FLUID_SOLVER_MULTIPHASE_H
#define FLUID_SOLVER_MULTIPHASE_H

#include "FluidSolver.h"


///@brief Experimental solver supporting fluid-fluid interaction.
///@remarks
///This solver has issues when FluidSph with differing FluidParametersLocal interact:
/// - Fluid particles will stick to boundaries; for instance, some particles of lighter fluids
///will not rise even if heavier fluid particles are on top.
/// - Fluid particles with greater mass / rest density will rise, 
///while 'lighter' fluids will sink. 
class FluidSolverMultiphase : public FluidSolver
{
public:
	virtual void stepSimulation(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids);
	
protected:
	virtual void sphComputePressure(const FluidParametersGlobal &FG, FluidSph *fluid, btAlignedObjectArray<FluidSph*> *interactingFluids);
	virtual void sphComputeForce(const FluidParametersGlobal &FG, FluidSph *fluid, btAlignedObjectArray<FluidSph*> *interactingFluids);
};

#endif