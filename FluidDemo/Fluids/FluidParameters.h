/*
  FLUIDS v.1 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

  ZLib license
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
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

#ifndef FLUID_PARAMETERS_H
#define FLUID_PARAMETERS_H
	
#include "LinearMath/btVector3.h"


///@brief Contains characteristics shared by all fluids inside a FluidWorld.
struct FluidParametersGlobal
{
	btVector3 m_planeGravity;			///<Simulation scale; meters / seconds^2.
	btVector3 m_pointGravityPosition;	///<World scale; meters.
	btScalar m_pointGravity;			///<Simulation scale; meters / seconds^2.
	
	btScalar m_timeStep;				///<Seconds; simulation becomes unstable at > ~0.004s timestep( with setDefaultParameters() ).
	
	///As SPH fluid simulations are scale sensitive, the simulation
	///(fluid-fluid interaction) is performed at a physically-correct
	///'simulation scale', which is typically much smaller than the 
	///'world scale' at which the particles are rendered.
	btScalar m_simulationScale;			///<N*m_simulationScale converts N into simulation scale; N/m_simulationScale converts N into world scale.
	btScalar m_particleRadius;			///<For collsion detection/integration; simulation scale; meters.
	btScalar m_speedLimit;				///<Acceleration/force limit; simulation scale; meters/second.
	btScalar m_sphSmoothRadius;			///<SPH particle interaction radius; use setSphInteractionRadius() to set this; simulation scale; meters.
	
	//Kernel function coefficients; dependent on m_sphSmoothRadius
	btScalar m_sphRadiusSquared;		///<m_sphSmoothRadius^2.
	btScalar m_poly6KernCoeff;			///<Coefficient of the poly6 kernel; for density calculation.
	btScalar m_spikyKernGradCoeff;		///<Coefficient of the gradient of the spiky kernel; for pressure force calculation.
	btScalar m_viscosityKernLapCoeff;	///<Coefficient of the Laplacian of the viscosity kernel; for viscosity force calculation.
	
	FluidParametersGlobal() { setDefaultParameters(); }
	void setDefaultParameters()
	{
		m_planeGravity.setValue(0, btScalar(-9.8), 0);
		m_pointGravityPosition.setValue(0, 0, 0);
		m_pointGravity = btScalar(0.0);

		m_timeStep = btScalar(0.003); 	//0.001 == for point grav
		
		m_simulationScale 	 = btScalar(0.004);
		m_particleRadius 	 = btScalar(0.004);
		m_speedLimit 		 = btScalar(200.0);	
		
		setSphInteractionRadius( btScalar(0.01) );
	}
	
	void setSphInteractionRadius(btScalar radius)
	{
		m_sphSmoothRadius = radius;
	
		m_sphRadiusSquared = m_sphSmoothRadius * m_sphSmoothRadius;
		
		//Wpoly6 kernel (denominator part) - 2003 Muller, p.4
		m_poly6KernCoeff = btScalar(315.0) / ( btScalar(64.0) * SIMD_PI * btPow(m_sphSmoothRadius, 9) );
		
		m_spikyKernGradCoeff = btScalar(-45.0) / ( SIMD_PI * btPow(m_sphSmoothRadius, 6) );
		
		//Laplacian of viscocity (denominator): PI h^6
		m_viscosityKernLapCoeff = btScalar(45.0) / ( SIMD_PI * btPow(m_sphSmoothRadius, 6) );
	}
};

///@brief Contains the properties of a single FluidSph.
struct FluidParametersLocal
{
	btVector3 m_volumeMin;				///<Particles cannot move below this boundary; world scale; meters.
	btVector3 m_volumeMax;				///<Particles cannot move above this boundary; world scale; meters.
	
	btScalar m_viscosity;				///<Higher values increase the fluid's resistance to flow; force calculation; pascal*seconds(Pa*s).
	btScalar m_restDensity;				///<Used for pressure calculation; kilograms/meters^3
	btScalar m_particleMass;			///<Used for density calculation and collision response; kilograms.
	btScalar m_stiffness;				///<Gas constant; higher values make a less compressible, more unstable fluid; pressure calculation; joules.
	
	btScalar m_boundaryStiff;			///<Spring coefficient; controls the magnitude of the boundary repulsion force.
	btScalar m_boundaryDamp;			///<Damping coefficient; controls the influence of relative velocity on the boundary repulsion force.
	btScalar m_boundaryFriction;		///<Fraction of tangential velocity removed per frame; [0.0, 1.0]; higher values more unstable.
	
	btScalar m_particleDist;			///<Used to determine particle spacing for FluidEmitter; simulation scale; meters. 
	
	FluidParametersLocal() { setDefaultParameters(); }
	void setDefaultParameters()
	{
		m_viscosity 	= btScalar(0.2);
		m_restDensity 	= btScalar(600.0);
		m_particleMass 	= btScalar(0.00020543);
		m_stiffness 	= btScalar(0.5);
		
		m_boundaryStiff	= btScalar(20000.0);
		m_boundaryDamp 	= btScalar(256.0);
		m_boundaryFriction 	= btScalar(0.0);
		
		m_particleDist = btPow( m_particleMass/m_restDensity, btScalar(1.0/3.0) );
	}
};


#endif


