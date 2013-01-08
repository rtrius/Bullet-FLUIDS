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
//Portions of this file based on FLUIDS v.2 - SPH Fluid Simulator for CPU and GPU
//Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com
#ifndef BT_FLUID_SPH_PARAMETERS_H
#define BT_FLUID_SPH_PARAMETERS_H
	
#include "LinearMath/btVector3.h"


///@brief Contains characteristics shared by all btFluidSph inside a btFluidRigidDynamicsWorld.
struct btFluidSphParametersGlobal
{
	btScalar m_timeStep;				///<Seconds; simulation becomes unstable at > ~0.004s timestep( with setDefaultParameters() ).
	
	///@brief N*m_simulationScale converts N into simulation scale; N/m_simulationScale converts N into world scale.
	///@remarks As SPH fluid simulations are scale sensitive, the simulation
	///(fluid-fluid interaction) is performed at a physically-correct
	///'simulation scale', which is typically much smaller than the 
	///'world scale' at which the particles are rendered.
	///@remarks Note that the grid cell size is dependent on the simulation scale, 
	///so btFluidSph::setGridCellSize() should also be called after changing this.
	btScalar m_simulationScale;
	btScalar m_speedLimit;				///<Acceleration/force limit; simulation scale; meters/second.
	btScalar m_sphSmoothRadius;			///<SPH particle interaction radius; use setSphInteractionRadius() to set this; simulation scale; meters.
	
	///@name Kernel function coefficients; dependent on m_sphSmoothRadius; use setSphInteractionRadius() to set these.
	///@{
	btScalar m_sphRadiusSquared;		///<m_sphSmoothRadius^2.
	btScalar m_poly6KernCoeff;			///<Coefficient of the poly6 kernel; for density calculation.
	btScalar m_spikyKernGradCoeff;		///<Coefficient of the gradient of the spiky kernel; for pressure force calculation.
	btScalar m_viscosityKernLapCoeff;	///<Coefficient of the Laplacian of the viscosity kernel; for viscosity force calculation.
	btScalar m_initialSum; 				///<Self-contributed particle density; should generally be within [0.0, m_sphRadiusSquared^3] (for Wpoly6).
	///@}
	
	btFluidSphParametersGlobal() { setDefaultParameters(); }
	void setDefaultParameters()
	{
		m_timeStep = btScalar(0.003);
		
		m_simulationScale 	 = btScalar(0.004);
		m_speedLimit 		 = btScalar(200.0);	
		
		setSphInteractionRadius( btScalar(0.01) );
	}
	
	///The grid cell size is dependent on this radius, so btFluidSph::setGridCellSize() should also be called after this.
	void setSphInteractionRadius(btScalar radius)
	{
		m_sphSmoothRadius = radius;
	
		m_sphRadiusSquared = m_sphSmoothRadius * m_sphSmoothRadius;
		
		//Wpoly6 kernel (denominator part) - 2003 Muller, p.4
		m_poly6KernCoeff = btScalar(315.0) / ( btScalar(64.0) * SIMD_PI * btPow(m_sphSmoothRadius, 9) );
		
		m_spikyKernGradCoeff = btScalar(-45.0) / ( SIMD_PI * btPow(m_sphSmoothRadius, 6) );
		
		//Laplacian of viscocity (denominator): PI h^6
		m_viscosityKernLapCoeff = btScalar(45.0) / ( SIMD_PI * btPow(m_sphSmoothRadius, 6) );
		
		//
		//m_initialSum = btScalar(0.0);
		//m_initialSum = m_sphRadiusSquared*m_sphRadiusSquared*m_sphRadiusSquared;	//poly6 kernel partial result
		m_initialSum = m_sphRadiusSquared*m_sphRadiusSquared*m_sphRadiusSquared * btScalar(0.25);
	}
};

///@brief Contains the properties of a single btFluidSph.
struct btFluidSphParametersLocal
{
	btVector3 m_aabbBoundaryMin;		///<Particles cannot move below this boundary; world scale; meters.
	btVector3 m_aabbBoundaryMax;		///<Particles cannot move above this boundary; world scale; meters.
	int m_enableAabbBoundary;			///<If nonzero, the particles are confined to m_aabbBoundaryMin and m_aabbBoundaryMax.
	
	btVector3 m_gravity;				///<Simulation scale; meters / seconds^2.
	
	btScalar m_viscosity;				///<Higher values increase the fluid's resistance to flow; force calculation; pascal*seconds(Pa*s).
	btScalar m_restDensity;				///<Used for pressure calculation; kilograms/meters^3
	btScalar m_sphParticleMass;			///<Mass of a single particle when calculating SPH density and force; kilograms.
	btScalar m_stiffness;				///<Gas constant; higher values make a less compressible, more unstable fluid; pressure calculation; joules.
	
	btScalar m_particleDist;			///<Used to determine particle spacing for btFluidEmitter; simulation scale; meters. 
	btScalar m_particleRadius;			///<For collision detection and collision response; world scale; meters.
	btScalar m_particleMass;			///<Mass of a single particle when colliding with rigid bodies and applying forces; kilograms.
	
	btScalar m_boundaryStiff;			///<Spring coefficient; controls the magnitude of the boundary repulsion force.
	btScalar m_boundaryDamp;			///<Damping coefficient; controls the influence of relative velocity on the boundary repulsion force.
	btScalar m_boundaryFriction;		///<Fraction of tangential velocity removed per frame; [0.0, 1.0]; higher values more unstable.
	btScalar m_boundaryRestitution;		///<Fraction of reflected velocity(bounciness); [0.0, 1.0]; higher values more unstable.
	btScalar m_boundaryErp;				///<Controls how quickly penetration is removed(per frame impulse: penetration_depth*m_boundaryErp).
	
	btFluidSphParametersLocal() { setDefaultParameters(); }
	void setDefaultParameters()
	{
		m_aabbBoundaryMin.setValue(-BT_LARGE_FLOAT, -BT_LARGE_FLOAT, -BT_LARGE_FLOAT);
		m_aabbBoundaryMax.setValue(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);
		m_enableAabbBoundary = 0;
	
		m_gravity.setValue(0, btScalar(-9.8), 0);
	
		m_viscosity 	= btScalar(0.2);
		m_restDensity 	= btScalar(600.0);
		m_sphParticleMass  = btScalar(0.00020543);
		m_stiffness 	= btScalar(1.5);
		
		m_particleDist = btPow( m_sphParticleMass/m_restDensity, btScalar(1.0/3.0) );
		m_particleRadius = btScalar(1.0);
		m_particleMass = btScalar(0.00020543);
		
		m_boundaryStiff	= btScalar(20000.0);
		m_boundaryDamp 	= btScalar(256.0);
		m_boundaryFriction 	= btScalar(0.0);
		m_boundaryRestitution = btScalar(0.0);
		m_boundaryErp = btScalar(0.0375);
	}
};


#endif


