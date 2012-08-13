/* FluidParticles.h
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
#ifndef FLUID_PARTICLES_H
#define FLUID_PARTICLES_H

#include "LinearMath/btAlignedObjectArray.h"

class btVector3;

///@brief Stores a table of particle indicies and their distances from a single fluid particle.
///@remarks
///Only particles within the SPH interaction radius are included.
///@par
///The table is generated during the pressure calculation step in order to avoid recalculating
///neighboring particles and their distances during force computation.
class FluidNeighbors
{
	static const int MAX_NEIGHBORS = 80;

	int m_count;
	int m_particleIndicies[MAX_NEIGHBORS];
	btScalar m_distances[MAX_NEIGHBORS];
	
public:
	inline int numNeighbors() const { return m_count; }
	inline int getNeighborIndex(int index) const { return m_particleIndicies[index]; }
	inline btScalar getDistance(int index) const { return m_distances[index]; }
	inline bool isFilled() const { return (m_count >= MAX_NEIGHBORS); }
	
	inline void clear() { m_count = 0; }
	inline void addNeighbor(int neighborIndex, btScalar distance)
	{
		m_particleIndicies[m_count] = neighborIndex;
		m_distances[m_count] = distance;
		++m_count;
	}
};

const int INVALID_PARTICLE_INDEX = -1;

//Value lower than INVALID_PARTICLE_INDEX
//This is required in order to use the optimized loop for iterating through
//a FluidGridIterator when using FluidSortingGrid.
const int INVALID_PARTICLE_INDEX_MINUS_ONE = INVALID_PARTICLE_INDEX - 1; 

///@brief Coordinates the parallel arrays used to store fluid particles.
///@remarks
///Members of this struct should not be accessed directly, except for calling:
/// - btAlignedObjectArray::operator[]() or
/// - btAlignedObjectArray::at()
struct FluidParticles
{
	int m_maxParticles;

	//Parallel arrays
	btAlignedObjectArray<btVector3> m_pos;					///<Current position; world scale; meters.
	btAlignedObjectArray<btVector3> m_vel;					///<'Current + (1/2)*timestep' velocity for leapfrog integration; simulation scale.
	btAlignedObjectArray<btVector3> m_vel_eval;				///<Current velocity; simulation scale; meters.
	btAlignedObjectArray<btVector3> m_sph_force;			///<Sum of pressure and viscosity forces; simulation scale; meters.
	btAlignedObjectArray<btVector3> m_externalAcceleration;	///<Applied during FluidWorld::stepSimulation(), then set to 0; simulation scale; meters.
	btAlignedObjectArray<btScalar> m_pressure;				///<Value of the pressure scalar field at the particle's position.
	btAlignedObjectArray<btScalar> m_density;				///<Value of the density scalar field at the particle's position.
	btAlignedObjectArray<int> m_nextFluidIndex;				///<Index of the next fluid particle in the same grid cell(forward linked list).
	
	btAlignedObjectArray<FluidNeighbors> m_neighborTable;

	
	FluidParticles() : m_maxParticles(0) {}
	
	int	size() const	{ return m_pos.size(); }

	int addParticle(const btVector3 &position);
	void removeParticle(int index);		///<Swaps indicies if index does not correspond to the last index; invalidates grid.
	
	void resize(int newSize);		///<Does not initialize particles if( newSize > size() ).
	
	void setMaxParticles(int maxNumParticles);
	int getMaxParticles() const { return m_maxParticles; }
};


#endif


