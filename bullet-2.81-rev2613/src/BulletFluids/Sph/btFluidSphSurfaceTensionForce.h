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
#ifndef BT_FLUID_SPH_SURFACE_TENSION_FORCE
#define BT_FLUID_SPH_SURFACE_TENSION_FORCE

#include "btFluidSph.h"

class btFluidSphNeighbors;

///Computes and applies a surface tension force
///@remarks Based on:  \n
///"Versatile Surface Tension and Adhesion for SPH Fluids".  \n
///N. Akinci, G. Akinci, and M. Teschner.  \n
///ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2013), Nov. 2013.
class btFluidSphSurfaceTensionForce
{
protected:
	btAlignedObjectArray<btVector3> m_colorFieldGradient;
	btAlignedObjectArray<btVector3> m_surfaceTensionForce;
	
public:

	///Can only be called after particle neighbors and density are computed, and should only be called in the SPH solver.
	void computeAndApplySurfaceTensionForce(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, 
											const btAlignedObjectArray<btFluidSphNeighbors>& neighbors, 
											const btAlignedObjectArray<btScalar>& invDensity,
											btAlignedObjectArray<btVector3>& out_accumulatedForce)
	{
		int numParticles = fluid->numParticles();
		const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	
		m_colorFieldGradient.resize(numParticles);
		m_surfaceTensionForce.resize(numParticles);
		
		//
		computeColorFieldGradient(FG, fluid, neighbors, invDensity);
		computeSurfaceTensionForce(FG, fluid, neighbors, invDensity);
		
		//Apply surface tension force
		for(int i = 0; i < numParticles; ++i) out_accumulatedForce[i] += m_surfaceTensionForce[i];
	}

protected:
	void computeColorFieldGradient(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, 
									const btAlignedObjectArray<btFluidSphNeighbors>& neighbors, 
									const btAlignedObjectArray<btScalar>& invDensity);

	void computeSurfaceTensionForce(const btFluidSphParametersGlobal& FG, btFluidSph* fluid, 
									const btAlignedObjectArray<btFluidSphNeighbors>& neighborTables, 
									const btAlignedObjectArray<btScalar>& invDensity);
};

#endif
