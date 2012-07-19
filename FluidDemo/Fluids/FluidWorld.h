/** FluidWorld.h
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
#ifndef FLUID_WORLD_H
#define FLUID_WORLD_H

#include "FluidSph.h"

#ifdef FLUIDS_OPENCL_ENABLED	//#define in "FluidSph.h"
	#include "OpenCL_support/fluids_opencl_support.h"
#endif

class FluidWorld
{
	FluidParametersGlobal 			m_globalParameters;

	btAlignedObjectArray<FluidSph*> m_fluids;

	bool							m_useOpenCL;
	bool							m_toggledLastFrameOpenCL;
	
#ifdef FLUIDS_OPENCL_ENABLED	//#define in "FluidSph.h"
	FluidSystem_OpenCL 				m_openClSystem;
#endif
	
public:
	FluidWorld();
	~FluidWorld();
	
	const FluidParametersGlobal& getGlobalParameters() const { return m_globalParameters; }
	void setGlobalParameters(const FluidParametersGlobal &FG) { m_globalParameters = FG; }
	
	void addFluid(FluidSph *fluid) 
	{
		//btAssert( fluid && m_fluids.findLinearSearch(fluid) == m_fluids.size() );
		m_fluids.push_back(fluid);
	}
	void removeFluid(FluidSph *fluid) { m_fluids.remove(fluid);	}	//May swap elements in m_fluids
	int getNumFluids() const { return m_fluids.size(); }
	FluidSph* getFluid(int index) { return m_fluids[index]; }
	
	void stepSimulation();
	
	void toggleOpenCL() 
	{
		m_useOpenCL = !m_useOpenCL; 
		m_toggledLastFrameOpenCL = true;
	}
};


#endif
