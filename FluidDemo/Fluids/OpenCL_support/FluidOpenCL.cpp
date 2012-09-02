/* FluidOpenCL.cpp
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
#include "FluidOpenCL.h"

#include "../FluidParameters.h"
#include "../FluidParticles.h"
#include "../FluidSph.h"

void Fluid_OpenCL::writeToOpenCL(cl_context context, cl_command_queue commandQueue, 
								 const FluidParametersLocal &FL, FluidParticles *particles, bool transferAllData)
{
	if(m_maxParticles != particles->m_maxParticles)
	{
		deallocate();
		allocate(context, particles->m_maxParticles);
	}
	
	m_buffer_localParameters.writeToBuffer( commandQueue, &FL, sizeof(FluidParametersLocal) );
	
	int numParticles = particles->size();
	
	m_buffer_pos.writeToBuffer( commandQueue, &particles->m_pos[0], sizeof(btVector3)*numParticles );
	m_buffer_vel.writeToBuffer( commandQueue, &particles->m_vel[0], sizeof(btVector3)*numParticles );
	m_buffer_vel_eval.writeToBuffer( commandQueue, &particles->m_vel_eval[0], sizeof(btVector3)*numParticles );
	
	m_buffer_externalAcceleration.writeToBuffer( commandQueue, &particles->m_externalAcceleration[0], sizeof(btVector3)*numParticles );
	
	if(transferAllData)
	{
		m_buffer_sph_force.writeToBuffer( commandQueue, &particles->m_sph_force[0], sizeof(btVector3)*numParticles );
		m_buffer_pressure.writeToBuffer( commandQueue, &particles->m_pressure[0], sizeof(btScalar)*numParticles );
		m_buffer_invDensity.writeToBuffer( commandQueue, &particles->m_invDensity[0], sizeof(btScalar)*numParticles );
		m_buffer_nextFluidIndex.writeToBuffer( commandQueue, &particles->m_nextFluidIndex[0], sizeof(int)*numParticles );
		m_buffer_neighborTable.writeToBuffer( commandQueue, &particles->m_neighborTable[0], sizeof(FluidNeighbors)*numParticles );
	}
	
	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}

void Fluid_OpenCL::readFromOpenCL(cl_context context, cl_command_queue commandQueue,
								  FluidParticles *particles, bool transferAllData)
{
	btAssert(m_maxParticles != 0);

	int numParticles = particles->size();

	m_buffer_pos.readFromBuffer( commandQueue, &particles->m_pos[0], sizeof(btVector3)*numParticles );
	m_buffer_vel.readFromBuffer( commandQueue, &particles->m_vel[0], sizeof(btVector3)*numParticles );
	m_buffer_vel_eval.readFromBuffer( commandQueue, &particles->m_vel_eval[0], sizeof(btVector3)*numParticles );
	
	m_buffer_nextFluidIndex.readFromBuffer( commandQueue, &particles->m_nextFluidIndex[0], sizeof(int)*numParticles );
	
	if(transferAllData)
	{
		m_buffer_sph_force.readFromBuffer( commandQueue, &particles->m_sph_force[0], sizeof(btVector3)*numParticles );
		m_buffer_externalAcceleration.readFromBuffer( commandQueue, &particles->m_externalAcceleration[0], sizeof(btVector3)*numParticles );
		m_buffer_pressure.readFromBuffer( commandQueue, &particles->m_pressure[0], sizeof(btScalar)*numParticles );
		m_buffer_invDensity.readFromBuffer( commandQueue, &particles->m_invDensity[0], sizeof(btScalar)*numParticles );
		m_buffer_neighborTable.readFromBuffer( commandQueue, &particles->m_neighborTable[0], sizeof(FluidNeighbors)*numParticles );
	}
	
	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}

Fluid_OpenCLPointers Fluid_OpenCL::getPointers()
{
	Fluid_OpenCLPointers result;
	
	result.m_buffer_localParameters = m_buffer_localParameters.getAddress();
	
	result.m_buffer_pos = m_buffer_pos.getAddress();
	result.m_buffer_vel = m_buffer_vel.getAddress();
	result.m_buffer_vel_eval = m_buffer_vel_eval.getAddress();
	result.m_buffer_sph_force = m_buffer_sph_force.getAddress();
	result.m_buffer_externalAcceleration = m_buffer_externalAcceleration.getAddress();
	result.m_buffer_pressure = m_buffer_pressure.getAddress();
	result.m_buffer_invDensity = m_buffer_invDensity.getAddress();
	result.m_buffer_nextFluidIndex = m_buffer_nextFluidIndex.getAddress();
	
	result.m_buffer_neighborTable = m_buffer_neighborTable.getAddress();

	return result;
}

void Fluid_OpenCL::allocate(cl_context context, int maxParticles)
{
	m_maxParticles = maxParticles;
	
	m_buffer_localParameters.allocate( context, sizeof(FluidParametersLocal) );
	
	m_buffer_pos.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_vel.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_vel_eval.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_sph_force.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_externalAcceleration.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_pressure.allocate( context, sizeof(btScalar) * maxParticles );
	m_buffer_invDensity.allocate( context, sizeof(btScalar) * maxParticles );
	m_buffer_nextFluidIndex.allocate( context, sizeof(int) * maxParticles );
	
	m_buffer_neighborTable.allocate( context, sizeof(FluidNeighbors) * maxParticles );
}
void Fluid_OpenCL::deallocate()
{	
	m_maxParticles = 0;
	
	m_buffer_localParameters.deallocate();
	
	m_buffer_pos.deallocate();
	m_buffer_vel.deallocate();
	m_buffer_vel_eval.deallocate();
	m_buffer_sph_force.deallocate();
	m_buffer_externalAcceleration.deallocate();
	m_buffer_pressure.deallocate();
	m_buffer_invDensity.deallocate();
	m_buffer_nextFluidIndex.deallocate();
	
	m_buffer_neighborTable.deallocate();
}
