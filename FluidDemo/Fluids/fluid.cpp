/** fluid.cpp

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

#include "fluid.h"

int Fluids::addFluid(const btVector3 &position)
{
	int index;
	if( size() < m_maxParticles )
	{
		m_pos.push_back( btVector3() );
		m_vel.push_back( btVector3() );
		m_vel_eval.push_back( btVector3() );
		m_sph_force.push_back( btVector3() );
		m_externalAcceleration.push_back( btVector3() );
		m_prev_pos.push_back( btVector3() );
		m_pressure.push_back(0);
		m_density.push_back(0);
		m_nextFluidIndex.push_back(INVALID_PARTICLE_INDEX);
		
		m_neighborTable.push_back( Neighbors() );
		
		index = size() - 1;
	}
	else
	{
		index = ( size() - 1 ) * rand() / RAND_MAX;		//Random index
	}
	
	m_pos[index] = position;
	m_prev_pos[index] = position;		//CCD_TEST
	m_vel[index].setValue(0,0,0);
	m_vel_eval[index].setValue(0,0,0);
	m_sph_force[index].setValue(0,0,0);
	m_externalAcceleration[index].setValue(0,0,0);
	m_pressure[index] = 0;
	m_density[index] = 0;
	m_nextFluidIndex[index] = INVALID_PARTICLE_INDEX;
	
	m_neighborTable[index].clear();
	
	return index;
}
void Fluids::removeFluid(int index)
{
	if( index >= size() ) return;

	int lastIndex = size() - 1;
	
	if(index < lastIndex) 
	{
		m_pos[index] = m_pos[lastIndex];
		m_vel[index] = m_vel[lastIndex];
		m_vel_eval[index] = m_vel_eval[lastIndex];
		m_sph_force[index] = m_sph_force[lastIndex];
		m_externalAcceleration[index] = m_externalAcceleration[lastIndex];
		m_prev_pos[index] = m_prev_pos[lastIndex];
		m_pressure[index] = m_pressure[lastIndex];
		m_density[index] = m_density[lastIndex];
		m_nextFluidIndex[index] = m_nextFluidIndex[lastIndex];
		
		m_neighborTable[index] = m_neighborTable[lastIndex];
	}
	m_pos.pop_back();
	m_vel.pop_back();
	m_vel_eval.pop_back();
	m_sph_force.pop_back();
	m_externalAcceleration.pop_back();
	m_prev_pos.pop_back();
	m_pressure.pop_back();
	m_density.pop_back();
	m_nextFluidIndex.pop_back();
	
	m_neighborTable.pop_back();
}

void Fluids::resize(int size)
{
	if(size > m_maxParticles) m_maxParticles = size;

	m_pos.resize(size);
	m_vel.resize(size);
	m_vel_eval.resize(size);
	m_sph_force.resize(size);
	m_externalAcceleration.resize(size);
	m_prev_pos.resize(size);
	m_pressure.resize(size);
	m_density.resize(size);
	m_nextFluidIndex.resize(size);
	
	m_neighborTable.resize(size);
}

void Fluids::setMaxParticles(int maxNumParticles)
{
	m_maxParticles = maxNumParticles;
	
	m_pos.reserve(maxNumParticles);
	m_vel.reserve(maxNumParticles);
	m_vel_eval.reserve(maxNumParticles);
	m_sph_force.reserve(maxNumParticles);
	m_externalAcceleration.reserve(maxNumParticles);
	m_prev_pos.reserve(maxNumParticles);
	m_pressure.reserve(maxNumParticles);
	m_density.reserve(maxNumParticles);
	m_nextFluidIndex.reserve(maxNumParticles);
	
	m_neighborTable.reserve(maxNumParticles);
}