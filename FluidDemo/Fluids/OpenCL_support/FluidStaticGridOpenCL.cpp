/* FluidStaticGridOpenCL.cpp
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
#include "FluidStaticGridOpenCL.h"

#include "../FluidStaticGrid.h"

void FluidStaticGrid_OpenCL::writeToOpenCL(cl_context context, cl_command_queue commandQueue, FluidStaticGrid *grid)
{
	FluidStaticGridParameters GP = grid->getParameters();
	
	int currentNumCells = GP.m_numCells;
	if(m_numGridCells != currentNumCells)
	{
		deallocate();
		allocate(context, currentNumCells);
	}
	
	m_buffer_gridParams.writeToBuffer( commandQueue, &GP, sizeof(FluidStaticGridParameters) );
	
	m_buffer_gridCells.writeToBuffer( commandQueue, grid->internalGetCellsPointer(), sizeof(int)*currentNumCells );
	m_buffer_gridCellsNumFluids.writeToBuffer( commandQueue, grid->internalGetCellsNumFluidsPointer(), sizeof(int)*currentNumCells );

	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}
void FluidStaticGrid_OpenCL::readFromOpenCL(cl_context context, cl_command_queue commandQueue, FluidStaticGrid *grid)
{
	btAssert(m_numGridCells != 0);
	
	int currentNumCells = grid->getParameters().m_numCells;
	
	m_buffer_gridCells.readFromBuffer( commandQueue, grid->internalGetCellsPointer(), sizeof(int)*currentNumCells );
	m_buffer_gridCellsNumFluids.readFromBuffer( commandQueue, grid->internalGetCellsNumFluidsPointer(), sizeof(int)*currentNumCells );

	cl_int error_code = clFinish(commandQueue);
	CHECK_CL_ERROR(error_code);
}

FluidStaticGrid_OpenCLPointers FluidStaticGrid_OpenCL::getPointers()
{
	FluidStaticGrid_OpenCLPointers result;
	
	result.m_buffer_gridParams = m_buffer_gridParams.getAddress();
	
	result.m_buffer_gridCells = m_buffer_gridCells.getAddress();
	result.m_buffer_gridCellsNumFluids = m_buffer_gridCellsNumFluids.getAddress();

	return result;
}

void FluidStaticGrid_OpenCL::allocate(cl_context context, int numGridCells)
{
	m_buffer_gridParams.allocate( context, sizeof(FluidStaticGridParameters) );
	
	m_numGridCells = numGridCells;
	m_buffer_gridCells.allocate( context, sizeof(int) * numGridCells );
	m_buffer_gridCellsNumFluids.allocate( context, sizeof(int) * numGridCells );
}
void FluidStaticGrid_OpenCL::deallocate()
{
	m_buffer_gridParams.deallocate();
	
	m_numGridCells = 0;
	m_buffer_gridCells.deallocate();
	m_buffer_gridCellsNumFluids.deallocate();
}