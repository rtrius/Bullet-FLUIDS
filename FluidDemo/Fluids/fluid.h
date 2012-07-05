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
	btScalar m_distances[MAX_NEIGHBORS];
	
public:
	inline unsigned short numNeighbors() const { return m_count; }
	inline unsigned short getNeighborIndex(int index) const { return m_particleIndicies[index]; }
	inline btScalar getDistance(int index) const { return m_distances[index]; }
	inline bool isFilled() const { return (m_count >= MAX_NEIGHBORS); }
	
	inline void clear() { m_count = 0; }
	inline void addNeighbor(unsigned short neighborIndex, btScalar distance)
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
	btAlignedObjectArray<btScalar> m_pressure;					//SPH
	btAlignedObjectArray<btScalar> m_density;					//SPH
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
	btScalar m_pointGravity;
	
	btScalar m_timeStep;				//Seconds; simulation becomes unstable at > ~0.004s timestep
	
	//
	btScalar sph_simscale;				//N*simscale converts N into simulation scale; N/simscale converts N into world scale
	btScalar sph_visc;					//Force calculation
	btScalar sph_restdensity;			//Pressure/density calculation
	btScalar sph_pmass;					//Pressure/density calculation
	btScalar sph_pradius;				//Meters; Collision detection/integration; sim scale
	btScalar sph_pdist;					//Meters; used to determine particle spacing
	btScalar sph_smoothradius;			//Simulation scale; Meters
	btScalar sph_intstiff;				//Pressure/density calculation
	btScalar sph_extstiff;				//Integration; sim scale
	btScalar sph_extdamp;				//Integration; sim scale
	btScalar sph_limit;					//Acceleration/force limit
	
	//Kernel functions (constant?)
	btScalar m_R2;
	btScalar m_Poly6Kern;
	btScalar m_LapKern;
	btScalar m_SpikyKern;	
};

#endif













