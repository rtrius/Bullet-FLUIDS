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

//Since SPH fluid simulation is scale sensitive, the simulation
//(fluid-fluid interaction) is performed at a physically-correct
//'simulation scale', which is typically much smaller than the 
//'world scale' at which the particles are rendered.
struct FluidParametersGlobal
{
	btVector3 m_planeGravity;
	btVector3 m_pointGravityPosition;
	btScalar m_pointGravity;
	
	btScalar m_timeStep;				//Seconds; simulation becomes unstable at > ~0.004s timestep
	
	//
	btScalar sph_simscale;				//N*simscale converts N into simulation scale; N/simscale converts N into world scale
	btScalar sph_pradius;				//Meters; Collision detection/integration; sim scale
	btScalar sph_smoothradius;			//Simulation scale; Meters
	btScalar sph_limit;					//Acceleration/force limit
	
	//Kernel function coefficients
	btScalar m_R2;
	btScalar m_Poly6Kern;				//Density calculation
	btScalar m_LapKern;					//Viscosity force calculation
	btScalar m_SpikyKern;				//Pressure force calculation
	
	FluidParametersGlobal() { setDefaultParameters(); }
	void setDefaultParameters()
	{
		m_planeGravity.setValue(0, btScalar(-9.8), 0);
		m_pointGravityPosition.setValue(0, 0, 0);
		m_pointGravity = btScalar(0.0);

		m_timeStep = btScalar(0.003); 	//0.001 == for point grav
		
		//SPH Parameters
		{
			sph_simscale 	 = btScalar(0.004);		//Unit size
			sph_pradius 	 = btScalar(0.004);		//m
			sph_smoothradius = btScalar(0.01);		//m 
			sph_limit 		 = btScalar(200.0);		//m/s
		}
		
		//SPH Kernels
		{
			m_R2 = sph_smoothradius * sph_smoothradius;
			
			//Wpoly6 kernel (denominator part) - 2003 Muller, p.4
			m_Poly6Kern = btScalar(315.0) / ( btScalar(64.0) * SIMD_PI * btPow(sph_smoothradius, 9) );
			
			//Laplacian of viscocity (denominator): PI h^6
			m_SpikyKern = btScalar(-45.0) / ( SIMD_PI * btPow(sph_smoothradius, 6) );
			
			m_LapKern = btScalar(45.0) / ( SIMD_PI * btPow(sph_smoothradius, 6) );
		}
	}
};

struct FluidParametersLocal
{
	btVector3 m_volumeMin;				//World scale
	btVector3 m_volumeMax;				//World scale
	
	btScalar m_viscosity;				//Force calculation
	btScalar m_restDensity;				//Pressure/density calculation
	btScalar m_particleMass;			//Pressure/density calculation
	btScalar m_intstiff;				//Pressure/density calculation
	btScalar m_extstiff;				//Integration; sim scale
	btScalar m_extdamp;					//Integration; sim scale
	
	btScalar m_particleDist;			//Meters; used to determine particle spacing
	
	FluidParametersLocal() { setDefaultParameters(); }
	void setDefaultParameters()
	{
		m_viscosity 	= btScalar(0.2);			//Pascal-second (Pa.s) = 1 kg m^-1 s^-1  (see wikipedia page on viscosity)
		m_restDensity 	= btScalar(600.0);		//kg / m^3
		m_particleMass 	= btScalar(0.00020543);	//kg
		m_intstiff 		= btScalar(0.5);
		m_extstiff		= btScalar(20000.0);
		m_extdamp 		= btScalar(256.0);
		
		m_particleDist = btPow( m_particleMass/m_restDensity, btScalar(1.0/3.0) );
	}
};


#endif


