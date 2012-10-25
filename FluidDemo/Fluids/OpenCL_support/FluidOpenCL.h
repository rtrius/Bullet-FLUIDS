/* FluidOpenCL.h
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
#ifndef FLUID_OPENCL
#define FLUID_OPENCL

#include "opencl_support.h"

struct FluidParametersLocal;
struct FluidParticles;

///@brief Manages OpenCL buffers corresponding to FluidParticles and a FluidParametersLocal.
class Fluid_OpenCL
{
	int m_maxParticles;
	
public:
	OpenCLBuffer m_buffer_localParameters;			//FluidParametersLocal
	
	OpenCLBuffer m_buffer_pos;						//btVector3[]
	OpenCLBuffer m_buffer_vel_eval;					//btVector3[]
	OpenCLBuffer m_buffer_sph_force;				//btVector3[]
	OpenCLBuffer m_buffer_pressure;					//float[]
	OpenCLBuffer m_buffer_invDensity;				//float[]
	
	OpenCLBuffer m_buffer_neighborTable;			//NeighborTable[]

	Fluid_OpenCL() : m_maxParticles(0) {}
	~Fluid_OpenCL() { deallocate(); }
	
	void writeToOpenCL(cl_context context, cl_command_queue commandQueue, const FluidParametersLocal &FL, FluidParticles *particles);
	void readFromOpenCL(cl_context context, cl_command_queue commandQueue, FluidParticles *particles);
	
private:
	void allocate(cl_context context, int maxParticles);
	void deallocate();
};

#endif
