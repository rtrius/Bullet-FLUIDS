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

///(Work in progress) Computes and applies a surface tension force
///@remarks Based on:  \n
///"Versatile Surface Tension and Adhesion for SPH Fluids".  \n
///N. Akinci, G. Akinci, and M. Teschner.  \n
///ACM Transactions on Graphics (Proc. SIGGRAPH Asia 2013), Nov. 2013.
class btFluidSphSurfaceTensionForce
{
protected:
	btAlignedObjectArray<btScalar> m_density;	//Default SPH solver only stores inverted density
	
	btAlignedObjectArray<btVector3> m_colorFieldGradient;
	btAlignedObjectArray<btVector3> m_surfaceTensionForce;
	
public:

	///Coefficients of the spline function C(); equation 2 in the paper
	struct Coefficients
	{
		btScalar m_C_multiply_coeff;		///< 32.0 / (pi * h^9)
		btScalar m_C_add_coeff;				///< h^6 / 64.0
		
		//
		btScalar m_sphSmoothRadiusHalved;
	};

	///Can only be called after particle neighbors and density are computed, and should only be called in the SPH solver.
	void computeAndApplySurfaceTensionForce(const btFluidSphParametersGlobal& FG, btFluidSph* fluid,
											const btAlignedObjectArray<btFluidSphNeighbors>& neighbors, 
											const btAlignedObjectArray<btScalar>& invDensity,
											btAlignedObjectArray<btVector3>& out_accumulatedAcceleration)
	{
		int numParticles = fluid->numParticles();
		const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
	
		//
		{
			m_density.resize(numParticles);
			m_colorFieldGradient.resize(numParticles);
			m_surfaceTensionForce.resize(numParticles);
		}
		
		//
		for(int i = 0; i < numParticles; ++i) m_density[i] = btScalar(1.0) / invDensity[i];
		
		//
		btFluidSphSurfaceTensionForce::Coefficients ST;
		ST.m_C_multiply_coeff = btScalar(32.0) / ( SIMD_PI * btPow(FG.m_sphSmoothRadius, btScalar(9.0)) );
		ST.m_C_add_coeff = btPow(FG.m_sphSmoothRadius, btScalar(6.0)) / btScalar(64.0);
		ST.m_sphSmoothRadiusHalved = FG.m_sphSmoothRadius * btScalar(0.5);
		
		//
		computeColorFieldGradient(FG, fluid, neighbors, invDensity);
		computeSurfaceTensionForce(FG, fluid, ST, neighbors, m_density);
		
		//Apply surface tension force
		for(int i = 0; i < numParticles; ++i) out_accumulatedAcceleration[i] += m_surfaceTensionForce[i] / FL.m_sphParticleMass;
	}

protected:
	void computeColorFieldGradient(const btFluidSphParametersGlobal& FG, btFluidSph* fluid,
									const btAlignedObjectArray<btFluidSphNeighbors>& neighbors, 
									const btAlignedObjectArray<btScalar>& invDensity);

	void computeSurfaceTensionForce(const btFluidSphParametersGlobal& FG, btFluidSph* fluid,
									const btFluidSphSurfaceTensionForce::Coefficients& ST,
									const btAlignedObjectArray<btFluidSphNeighbors>& neighborTables, 
									const btAlignedObjectArray<btScalar>& density);
};

#endif
