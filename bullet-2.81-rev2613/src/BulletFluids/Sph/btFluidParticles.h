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
#ifndef BT_FLUID_PARTICLES_H
#define BT_FLUID_PARTICLES_H

#include "LinearMath/btAlignedObjectArray.h"

class btVector3;

///@brief Coordinates the parallel arrays used to store fluid particles.
///@remarks
///Members of this struct should not be accessed directly, except for calling:
/// - btAlignedObjectArray::operator[]() or
/// - btAlignedObjectArray::at()
struct btFluidParticles
{
	int m_maxParticles;

	//Parallel arrays
	btAlignedObjectArray<btVector3> m_pos;					///<Current position; world scale.
	btAlignedObjectArray<btVector3> m_vel;					///<'Current + (1/2)*timestep' velocity for leapfrog integration; simulation scale.
	btAlignedObjectArray<btVector3> m_vel_eval;				///<Current velocity; simulation scale.
	btAlignedObjectArray<btVector3> m_accumulatedForce;		///<Applied during stepSimulation(), then set to 0; simulation scale.
	
	btAlignedObjectArray<void*> m_userPointer;
	
	btFluidParticles() : m_maxParticles(0) {}
	
	int	size() const	{ return m_pos.size(); }

	int addParticle(const btVector3& position);		///<Returns size() if size() == getMaxParticles().
	void removeParticle(int index);					///<Swaps indicies if index does not correspond to the last index; invalidates grid.
	
	void resize(int newSize);						///<Does not initialize particles if( newSize > size() ).
	
	void setMaxParticles(int maxNumParticles);
	int getMaxParticles() const { return m_maxParticles; }
};


#endif


