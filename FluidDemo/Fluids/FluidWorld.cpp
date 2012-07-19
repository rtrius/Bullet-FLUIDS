/** FluidWorld.cpp
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

#include "FluidWorld.h"

FluidWorld::FluidWorld() : m_useOpenCL(false), m_toggledLastFrameOpenCL(false)
{
	m_globalParameters.setDefaultParameters();
		
#ifdef FLUIDS_OPENCL_ENABLED	//#define in "FluidSph.h"
	m_openClSystem.initialize();
#endif
}
FluidWorld::~FluidWorld()
{
#ifdef FLUIDS_OPENCL_ENABLED	//#define in "FluidSph.h"
	m_openClSystem.deactivate();
#endif
}
	
void FluidWorld::stepSimulation()
{		
	if(m_useOpenCL || m_toggledLastFrameOpenCL) 	//GPU Branch
	{
#ifdef BT_USE_DOUBLE_PRECISION
		printf("BT_USE_DOUBLE_PRECISION not supported on OpenCL.\n");
		return;
#endif	

#ifdef FLUIDS_OPENCL_ENABLED 	//#define in "FluidSph.h"
		for(int i = 0; i < m_fluids.size(); ++i) 
		{
			if( !m_fluids[i]->numParticles() ) continue;
			
			m_openClSystem.stepSimulation(m_globalParameters, m_fluids[i], m_toggledLastFrameOpenCL);
		}
#endif
			
		m_toggledLastFrameOpenCL = false;
	}
	else 											//CPU Branch
	{
		for(int i = 0; i < m_fluids.size(); ++i) 
			m_fluids[i]->stepSimulation(m_globalParameters);
		
	}
}
