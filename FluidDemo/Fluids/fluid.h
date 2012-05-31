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

#ifndef DEF_FLUID
#define DEF_FLUID
	
#include "vector3df.h"
#include "LinearMath/btVector3.h"

struct Fluid 
{
	Vector3DF		pos;				//'current' position
	int				nextFluidIndex;		//Index of the next Fluid in the same grid cell(forward linked list of 'struct Fluid')
	
	Vector3DF		vel;				//'current + (1/2)*timestep' velocity; used for leapfrog integration
	Vector3DF		vel_eval;			//'current' velocity; used to compute forces, determine collisions
	
	float			pressure;			//Smoothed Particle Hydrodynamics
	float			density;	
	Vector3DF		sph_force;
	
	Vector3DF		externalAcceleration;		//This is applied in the next simulation step, then set to 0
	
	Vector3DF		prev_pos;			//CCD_TEST
};

class Neighbors
{
	static const int MAX_NEIGHBORS = 80;

	unsigned short m_count;
	unsigned short m_particleIndicies[MAX_NEIGHBORS];
	float m_distances[MAX_NEIGHBORS];
	
public:
	inline unsigned short numNeighbors() const { return m_count; }
	inline unsigned short getNeighborIndex(int index) const { return m_particleIndicies[index]; }
	inline float getDistance(int index) const { return m_distances[index]; }
	inline bool isFilled() const { return (m_count >= MAX_NEIGHBORS); }
	
	inline void clear() { m_count = 0; }
	inline void addNeighbor(unsigned short neighborIndex, float distance)
	{
		m_particleIndicies[m_count] = neighborIndex;
		m_distances[m_count] = distance;
		++m_count;
	}
};

//Since SPH fluid simulation is scale sensitive, the simulation
//(fluid-fluid interaction) is performed at a physically-correct
//'simulation scale', which is typically much smaller than the 
//'world scale' at which the particles are rendered.
struct FluidParameters
{
	double sph_simscale;				//N*simscale converts N into simulation scale; N/simscale converts N into world scale
	double sph_visc;					//Force calculation
	double sph_restdensity;				//Pressure/density calculation
	double sph_pmass;					//Pressure/density calculation
	double sph_pradius;					//Meters; Collision detection/integration; sim scale
	double sph_pdist;					//Meters;	used to determine particle spacing
	double sph_smoothradius;			//Simulation scale; Meters
	double sph_intstiff;				//Pressure/density calculation
	double sph_extstiff;				//Integration; sim scale
	double sph_extdamp;					//Integration; sim scale
	double sph_limit;					//Acceleration/force limit
	
	double m_timeStep;					//Seconds; simulation becomes unstable at > ~0.004s timestep
	
	//Kernel functions (constant?)
	double m_R2;
	double m_Poly6Kern;
	double m_LapKern;
	double m_SpikyKern;	
	
	//
	Vector3DF m_volumeMin;	//World scale
	Vector3DF m_volumeMax;	//World scale
	
	float m_pointGravity;
	Vector3DF m_pointGravityPosition;
	
	Vector3DF m_planeGravity;
};


///Since OpenCL support for doubles is optional,
///it is necessary to convert all doubles in 
///FluidParameters to float before sending it to OpenCL.
struct FluidParameters_float
{	
	float sph_simscale;
	float sph_visc;
	float sph_restdensity;
	float sph_pmass;
	float sph_pradius;
	float sph_pdist;
	float sph_smoothradius;
	float sph_intstiff;
	float sph_extstiff;
	float sph_extdamp;
	float sph_limit;
	float m_timeStep;
	float m_R2;
	float m_Poly6Kern;
	float m_LapKern;
	float m_SpikyKern;
	Vector3DF m_volumeMin;	//If btVector3 is used in place of Vector3DF, this may contain doubles
	Vector3DF m_volumeMax;
	float m_pointGravity;
	Vector3DF m_pointGravityPosition;
	Vector3DF m_planeGravity;
	
	FluidParameters_float(const FluidParameters &FP)
	{
		sph_simscale = FP.sph_simscale;
		sph_visc = FP.sph_visc;
		sph_restdensity = FP.sph_restdensity;
		sph_pmass = FP.sph_pmass;
		sph_pradius = FP.sph_pradius;
		sph_pdist = FP.sph_pdist;
		sph_smoothradius = FP.sph_smoothradius;
		sph_intstiff = FP.sph_intstiff;
		sph_extstiff = FP.sph_extstiff;
		sph_extdamp = FP.sph_extdamp;
		sph_limit = FP.sph_limit;
		m_timeStep = FP.m_timeStep;
		m_R2 = FP.m_R2;
		m_Poly6Kern = FP.m_Poly6Kern;
		m_LapKern = FP.m_LapKern;
		m_SpikyKern = FP.m_SpikyKern;
		m_volumeMin = FP.m_volumeMin;
		m_volumeMax = FP.m_volumeMax;
		m_pointGravity = FP.m_pointGravity;
		m_pointGravityPosition = FP.m_pointGravityPosition;
		m_planeGravity = FP.m_planeGravity;
	}
};

#endif













