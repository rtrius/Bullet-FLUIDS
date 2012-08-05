/** FluidSolver.h

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
#ifndef FLUID_SOLVER_H
#define FLUID_SOLVER_H

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

#include "FluidSph.h"


class FluidSolver
{
public:
	virtual void stepSimulation(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids) = 0;
	
protected:
	virtual void integrate(const FluidParametersGlobal &FG, const FluidParametersLocal &FL, FluidParticles *fluids);
};

class FluidSolverGridNeighbor : public FluidSolver
{
public:
	virtual void stepSimulation(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids)
	{
		BT_PROFILE("FluidSolverGridNeighbor::stepSimulation()");
	
		//
		for(int i = 0; i < fluids->size(); ++i)
		{
			(*fluids)[i]->removeMarkedParticles();
			(*fluids)[i]->insertParticlesIntoGrid();
		}
		
		//
		for(int i = 0; i < fluids->size(); ++i) 
			sphComputePressure( FG, (*fluids)[i] );
			
		for(int i = 0; i < fluids->size(); ++i) 
			sphComputeForce( FG, (*fluids)[i] );
			
		for(int i = 0; i < fluids->size(); ++i)
			integrate( FG, (*fluids)[i]->getLocalParameters(), &(*fluids)[i]->internalGetFluidParticles() );
	}
	
protected:
	virtual void sphComputePressure(const FluidParametersGlobal &FG, FluidSph *fluid);
	virtual void sphComputeForce(const FluidParametersGlobal &FG, FluidSph *fluid);
	
	void sphComputeForceGrid(const FluidParametersGlobal &FG, FluidSph *fluid);
};

class FluidSolverReducedGridNeighbor : public FluidSolverGridNeighbor
{
protected:
	virtual void sphComputePressure(const FluidParametersGlobal &FG, FluidSph *fluid) { sphComputePressureGridReduce(FG, fluid); }
	virtual void sphComputeForce(const FluidParametersGlobal &FG, FluidSph *fluid) { sphComputeForceReduce(FG, fluid); } 
	
	void sphComputePressureGridReduce(const FluidParametersGlobal &FG, FluidSph *fluid);
	void sphComputeForceReduce(const FluidParametersGlobal &FG, FluidSph *fluid);
};

#endif


