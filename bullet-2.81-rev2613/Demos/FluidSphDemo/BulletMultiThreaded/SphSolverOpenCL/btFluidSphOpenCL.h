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
#ifndef BT_FLUID_SPH_OPENCL
#define BT_FLUID_SPH_OPENCL

#include "btExperimentsOpenCL/btOpenCLArray.h"

class btVector3;
struct btFluidSphParametersLocal;
struct btFluidParticles;

///@brief Manages OpenCL buffers corresponding to btFluidParticles and btFluidSphParametersLocal.
class btFluidSphOpenCL
{
public:
	btOpenCLArray<btFluidSphParametersLocal> m_localParameters;
	
	btOpenCLArray<btVector3> m_pos;
	btOpenCLArray<btVector3> m_vel_eval;
	btOpenCLArray<btVector3> m_sph_force;
	btOpenCLArray<btScalar> m_density;
	btOpenCLArray<int> m_cellIndex;

	btFluidSphOpenCL(cl_context context, cl_command_queue queue) :
		m_localParameters(context, queue),
		m_pos(context, queue),
		m_vel_eval(context, queue),
		m_sph_force(context, queue),
		m_density(context, queue),
		m_cellIndex(context, queue) {}
	
	void writeToOpenCL(cl_command_queue queue, const btFluidSphParametersLocal& FL, btFluidParticles& particles);
	void readFromOpenCL(cl_command_queue queue, btAlignedObjectArray<btVector3>& sphForce);
};

#endif
