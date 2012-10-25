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

void Fluid_OpenCL::writeToOpenCL(cl_context context, cl_command_queue commandQueue, const FluidParametersLocal &FL, FluidParticles *particles)
{
	if(m_maxParticles != particles->m_maxParticles)
	{
		deallocate();
		allocate(context, particles->m_maxParticles);
	}
	
	m_buffer_localParameters.writeToBuffer( commandQueue, &FL, sizeof(FluidParametersLocal) );
	
	int numParticles = particles->size();
	m_buffer_pos.writeToBuffer( commandQueue, &particles->m_pos[0], sizeof(btVector3)*numParticles );
	m_buffer_vel_eval.writeToBuffer( commandQueue, &particles->m_vel_eval[0], sizeof(btVector3)*numParticles );
	
	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}

void Fluid_OpenCL::readFromOpenCL(cl_context context, cl_command_queue commandQueue, FluidParticles *particles)
{
	btAssert(m_maxParticles != 0);

	int numParticles = particles->size();
	m_buffer_sph_force.readFromBuffer( commandQueue, &particles->m_sph_force[0], sizeof(btVector3)*numParticles );
	
	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}

void Fluid_OpenCL::allocate(cl_context context, int maxParticles)
{
	m_maxParticles = maxParticles;
	
	m_buffer_localParameters.allocate( context, sizeof(FluidParametersLocal) );
	
	m_buffer_pos.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_vel_eval.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_sph_force.allocate( context, sizeof(btVector3) * maxParticles );
	m_buffer_pressure.allocate( context, sizeof(btScalar) * maxParticles );
	m_buffer_invDensity.allocate( context, sizeof(btScalar) * maxParticles );
	
	m_buffer_neighborTable.allocate( context, sizeof(FluidNeighbors) * maxParticles );
}
void Fluid_OpenCL::deallocate()
{	
	m_maxParticles = 0;
	
	m_buffer_localParameters.deallocate();
	
	m_buffer_pos.deallocate();
	m_buffer_vel_eval.deallocate();
	m_buffer_sph_force.deallocate();
	m_buffer_pressure.deallocate();
	m_buffer_invDensity.deallocate();
	
	m_buffer_neighborTable.deallocate();
}
