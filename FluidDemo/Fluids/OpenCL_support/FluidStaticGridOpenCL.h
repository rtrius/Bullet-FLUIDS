/* FluidStaticGridOpenCL.h
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
#ifndef FLUID_STATIC_GRID_OPENCL_H
#define FLUID_STATIC_GRID_OPENCL_H

#include "opencl_support.h"

class FluidStaticGrid;

struct FluidStaticGrid_OpenCLPointers
{
	void *m_buffer_gridParams;
	
	void *m_buffer_gridCells;
	void *m_buffer_gridCellsNumFluids;
};

///@brief Manages OpenCL buffers corresponding to a FluidStaticGrid.
class FluidStaticGrid_OpenCL
{
	OpenCLBuffer m_buffer_gridParams;			//FluidStaticGridParameters
	
	int m_numGridCells;
	OpenCLBuffer m_buffer_gridCells;			//int[]
	OpenCLBuffer m_buffer_gridCellsNumFluids;	//int[]	
	
public:
	FluidStaticGrid_OpenCL() : m_numGridCells(0) {}
	~FluidStaticGrid_OpenCL() { deallocate(); }
	
	void writeToOpenCL(cl_context context, cl_command_queue commandQueue, FluidStaticGrid *grid);
	void readFromOpenCL(cl_context context, cl_command_queue commandQueue, FluidStaticGrid *grid);
	
	FluidStaticGrid_OpenCLPointers getPointers();
	
private:
	void allocate(cl_context context, int numGridCells);
	void deallocate();
};

#endif