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
#include "btFluidParticles.h"

#include "LinearMath/btVector3.h"

int btFluidParticles::addParticle(const btVector3& position)
{
	if( size() < m_maxParticles )
	{
		m_pos.push_back( btVector3() );
		m_vel.push_back( btVector3() );
		m_vel_eval.push_back( btVector3() );
		m_accumulatedForce.push_back( btVector3() );
		m_userPointer.push_back(0);
		
		int index = size() - 1;
		
		m_pos[index] = position;
		m_vel[index].setValue(0,0,0);
		m_vel_eval[index].setValue(0,0,0);
		m_accumulatedForce[index].setValue(0,0,0);
		
		return index;
	}
	
	return size();
}
void btFluidParticles::removeParticle(int index)
{
	btAssert(0 <= index);
	btAssert( index < size() );
	
	int lastIndex = size() - 1;
	
	if(index < lastIndex) 
	{
		m_pos[index] = m_pos[lastIndex];
		m_vel[index] = m_vel[lastIndex];
		m_vel_eval[index] = m_vel_eval[lastIndex];
		m_accumulatedForce[index] = m_accumulatedForce[lastIndex];
		m_userPointer[index] = m_userPointer[lastIndex];
	}
	m_pos.pop_back();
	m_vel.pop_back();
	m_vel_eval.pop_back();
	m_accumulatedForce.pop_back();
	m_userPointer.pop_back();
}

void btFluidParticles::resize(int newSize)
{
	if(newSize > m_maxParticles) m_maxParticles = newSize;

	m_pos.resize(newSize);
	m_vel.resize(newSize);
	m_vel_eval.resize(newSize);
	m_accumulatedForce.resize(newSize);
	m_userPointer.resize(newSize);
}

void btFluidParticles::setMaxParticles(int maxNumParticles)
{
	m_maxParticles = maxNumParticles;
	
	m_pos.reserve(maxNumParticles);
	m_vel.reserve(maxNumParticles);
	m_vel_eval.reserve(maxNumParticles);
	m_accumulatedForce.reserve(maxNumParticles);
	m_userPointer.reserve(maxNumParticles);
}
