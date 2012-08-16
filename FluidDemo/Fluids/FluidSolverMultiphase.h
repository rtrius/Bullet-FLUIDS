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

inline bool isIntersectingAabb(const btVector3 &minA, const btVector3 &maxA, 
							   const btVector3 &minB, const btVector3 &maxB)
{
	return ( minA.x() <= maxB.x() && minB.x() <= maxA.x() 
		  && minA.y() <= maxB.y() && minB.y() <= maxA.y()
		  && minA.z() <= maxB.z() && minB.z() <= maxA.z() );
}

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
	virtual void stepSimulation(const FluidParametersGlobal &FG, btAlignedObjectArray<FluidSph*> *fluids)
	{
		BT_PROFILE("FluidSolverMultiphase::stepSimulation()");
	
		//
		for(int i = 0; i < fluids->size(); ++i)
		{
			(*fluids)[i]->removeMarkedParticles();
			(*fluids)[i]->insertParticlesIntoGrid();
		}
		
		//Determine intersecting FluidSph AABBs
		btAlignedObjectArray< btAlignedObjectArray<FluidSph*> > interactingFluids;
		interactingFluids.resize( fluids->size() );
		for(int i = 0; i < fluids->size(); ++i)
		{
			interactingFluids[i].resize(0);
		
			btVector3 min_i, max_i;
			(*fluids)[i]->getCurrentAabb(FG, &min_i, &max_i);
		
			for(int n = 0; n < fluids->size(); ++n)
			{
				if(i == n) continue;
				
				btVector3 min_n, max_n;
				(*fluids)[n]->getCurrentAabb(FG, &min_n, &max_n);
				
				if( isIntersectingAabb(min_i, max_i, min_n, max_n) )interactingFluids[i].push_back( (*fluids)[n] );
			}
		}
		
		//
		for(int i = 0; i < fluids->size(); ++i) 
			sphComputePressure( FG, (*fluids)[i], &interactingFluids[i] );
			
		for(int i = 0; i < fluids->size(); ++i) 
			sphComputeForce( FG, (*fluids)[i], &interactingFluids[i] );
			
		for(int i = 0; i < fluids->size(); ++i)
			integrate( FG, (*fluids)[i]->getLocalParameters(), &(*fluids)[i]->internalGetFluidParticles() );
	}
	
protected:
	virtual void sphComputePressure(const FluidParametersGlobal &FG, FluidSph *fluid, btAlignedObjectArray<FluidSph*> *interactingFluids);
	virtual void sphComputeForce(const FluidParametersGlobal &FG, FluidSph *fluid, btAlignedObjectArray<FluidSph*> *interactingFluids);
};

#endif