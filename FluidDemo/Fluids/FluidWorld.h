/* FluidWorld.h
	Copyright (C) 2012 Jackson Lee

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
#ifndef FLUID_WORLD_H
#define FLUID_WORLD_H

#include "FluidSph.h"
#include "FluidSolver.h"


///@brief Coordinates several FluidSph and global fluid properities.
///@remarks 
///Terminology:
/// - World - a set of fluids.
/// - Fluid - a set of particles.
/// - Particle - a point with position and velocity that may be influenced by other particles using SPH.
class FluidWorld
{
	FluidParametersGlobal 			m_globalParameters;

	btAlignedObjectArray<FluidSph*> m_fluids;

	FluidSolver 					*m_fluidSolver;
	
public:
	FluidWorld(FluidSolver *fluidSolver) : m_fluidSolver(fluidSolver) {}
	
	const FluidParametersGlobal& getGlobalParameters() const { return m_globalParameters; }
	void setGlobalParameters(const FluidParametersGlobal &FG) { m_globalParameters = FG; }
	
	void addFluid(FluidSph *fluid) 
	{
		btAssert(fluid);
		btAssert( m_fluids.findLinearSearch(fluid) == m_fluids.size() );
		
		m_fluids.push_back(fluid);
	}
	void removeFluid(FluidSph *fluid) { m_fluids.remove(fluid);	}	//May swap elements in m_fluids
	int getNumFluids() const { return m_fluids.size(); }
	FluidSph* getFluid(int index) { return m_fluids[index]; }
	
	void stepSimulation() { m_fluidSolver->stepSimulation(m_globalParameters, &m_fluids); }
	
	void setFluidSolver(FluidSolver *solver) { m_fluidSolver = solver; }
	FluidSolver* getFluidSolver() { return m_fluidSolver; }
};


#endif
