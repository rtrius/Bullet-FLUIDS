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
	
#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"

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

const int INVALID_PARTICLE_INDEX = -1;
struct Fluids
{
	int m_maxParticles;

	//Parallel arrays
	btAlignedObjectArray<btVector3> m_pos;						//Current position
	btAlignedObjectArray<btVector3> m_vel;						//'Current + (1/2)*timestep' velocity; used for leapfrog integration
	btAlignedObjectArray<btVector3> m_vel_eval;					//Current velocity
	btAlignedObjectArray<btVector3> m_sph_force;				//SPH
	btAlignedObjectArray<btVector3> m_externalAcceleration;		//This is applied during stepSimulation(), then set to 0
	btAlignedObjectArray<btVector3> m_prev_pos;					//	CCD_TEST
	btAlignedObjectArray<float> m_pressure;						//SPH
	btAlignedObjectArray<float> m_density;						//SPH
	btAlignedObjectArray<int> m_nextFluidIndex;					//Index of the next Fluid in the same grid cell(forward linked list)
	
	btAlignedObjectArray<Neighbors>  m_neighborTable;

	
	Fluids() : m_maxParticles(0) {}
	
	int	size() const	{ return m_pos.size(); }

	int addFluid(const btVector3 &position);
	void removeFluid(int index);		//Changes indicies; invalidates grid
	
	void resize(int size);		//Does not initialize particles if(size > current_size)
	
	void setMaxParticles(int maxNumParticles);
};

//Since SPH fluid simulation is scale sensitive, the simulation
//(fluid-fluid interaction) is performed at a physically-correct
//'simulation scale', which is typically much smaller than the 
//'world scale' at which the particles are rendered.
struct FluidParameters
{
	btVector3 m_volumeMin;				//World scale
	btVector3 m_volumeMax;				//World scale
	
	btVector3 m_planeGravity;
	btVector3 m_pointGravityPosition;
	float m_pointGravity;
	
	double m_timeStep;					//Seconds; simulation becomes unstable at > ~0.004s timestep
	
	//
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
	
	//Kernel functions (constant?)
	double m_R2;
	double m_Poly6Kern;
	double m_LapKern;
	double m_SpikyKern;	
};




//Since OpenCL support for doubles is optional, it is necessary to convert all 
//doubles in FluidParameters to float before sending it to OpenCL.
//
//btVector3 will contain doubles if BT_USE_DOUBLE_PRECISION is defined
//"LinearMath/btVector3.h" defines btVector3FloatData, but it is not aligned;
//the OpenCL implementation uses an aligned btVector3(float).
ATTRIBUTE_ALIGNED16(struct) btVector3FloatData_aligned { float m_floats[4]; };

struct FluidParameters_float
{	
	btVector3FloatData_aligned m_volumeMin;		
	btVector3FloatData_aligned m_volumeMax;
	btVector3FloatData_aligned m_planeGravity;
	btVector3FloatData_aligned m_pointGravityPosition;
	float m_pointGravity;
	float m_timeStep;
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
	float m_R2;
	float m_Poly6Kern;
	float m_LapKern;
	float m_SpikyKern;
	
	FluidParameters_float(const FluidParameters &FP)
	{
		for(int i = 0; i < 4; ++i) 
		{
			m_volumeMin.m_floats[i] = static_cast<float>(FP.m_volumeMin.m_floats[i]);
			m_volumeMax.m_floats[i] = static_cast<float>(FP.m_volumeMax.m_floats[i]);
			m_planeGravity.m_floats[i] = static_cast<float>(FP.m_planeGravity.m_floats[i]);
			m_pointGravityPosition.m_floats[i] = static_cast<float>(FP.m_pointGravityPosition.m_floats[i]);
		}
		m_pointGravity = static_cast<float>(FP.m_pointGravity);
		m_timeStep = static_cast<float>(FP.m_timeStep);
		sph_simscale = static_cast<float>(FP.sph_simscale);
		sph_visc = static_cast<float>(FP.sph_visc);
		sph_restdensity = static_cast<float>(FP.sph_restdensity);
		sph_pmass = static_cast<float>(FP.sph_pmass);
		sph_pradius = static_cast<float>(FP.sph_pradius);
		sph_pdist = static_cast<float>(FP.sph_pdist);
		sph_smoothradius = static_cast<float>(FP.sph_smoothradius);
		sph_intstiff = static_cast<float>(FP.sph_intstiff);
		sph_extstiff = static_cast<float>(FP.sph_extstiff);
		sph_extdamp = static_cast<float>(FP.sph_extdamp);
		sph_limit = static_cast<float>(FP.sph_limit);
		m_R2 = static_cast<float>(FP.m_R2);
		m_Poly6Kern = static_cast<float>(FP.m_Poly6Kern);
		m_LapKern = static_cast<float>(FP.m_LapKern);
		m_SpikyKern = static_cast<float>(FP.m_SpikyKern);
	}
};

#endif













