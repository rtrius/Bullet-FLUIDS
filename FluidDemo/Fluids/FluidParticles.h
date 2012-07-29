/** FluidParticles.h
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

class Neighbors
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
struct FluidParticles
{
	int m_maxParticles;

	//Parallel arrays
	btAlignedObjectArray<btVector3> m_pos;						//Current position
	btAlignedObjectArray<btVector3> m_vel;						//'Current + (1/2)*timestep' velocity; used for leapfrog integration
	btAlignedObjectArray<btVector3> m_vel_eval;					//Current velocity
	btAlignedObjectArray<btVector3> m_sph_force;				//SPH
	btAlignedObjectArray<btVector3> m_externalAcceleration;		//This is applied during stepSimulation(), then set to 0
	btAlignedObjectArray<btScalar> m_pressure;					//SPH
	btAlignedObjectArray<btScalar> m_density;					//SPH
	btAlignedObjectArray<int> m_nextFluidIndex;					//Index of the next Fluid in the same grid cell(forward linked list)
	
	btAlignedObjectArray<Neighbors>  m_neighborTable;

	
	FluidParticles() : m_maxParticles(0) {}
	
	int	size() const	{ return m_pos.size(); }

	int addParticle(const btVector3 &position);
	void removeParticle(int index);		//Changes indicies; invalidates grid
	
	void resize(int size);		//Does not initialize particles if(size > current_size)
	
	void setMaxParticles(int maxNumParticles);
	//int getMaxParticles() const { return m_maxParticles; }
};


#endif


