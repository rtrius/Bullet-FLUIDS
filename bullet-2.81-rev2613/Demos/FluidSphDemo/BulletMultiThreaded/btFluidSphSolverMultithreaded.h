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
#ifndef BT_FLUID_SPH_SOLVER_MULTITHREADED_H
#define BT_FLUID_SPH_SOLVER_MULTITHREADED_H

#include "btParallelFor.h"

#include "BulletFluids/Sph/btFluidSphSolver.h"

struct PF_ComputePressureData
{
	const btFluidSphParametersGlobal& m_globalParameters; 
	const btAlignedObjectArray<int>& m_gridCellGroup;
	const btFluidSortingGrid& m_grid;
	btFluidParticles& m_particles; 
	btFluidSphSolverDefault::SphParticles& m_sphData;
	
	PF_ComputePressureData(const btFluidSphParametersGlobal& FG, const btAlignedObjectArray<int>& gridCellGroup,
							const btFluidSortingGrid& grid, btFluidParticles& particles, btFluidSphSolverDefault::SphParticles& sphData) 
	: m_globalParameters(FG), m_gridCellGroup(gridCellGroup), m_grid(grid), m_particles(particles), m_sphData(sphData) {}
};
inline void PF_ComputePressureFunction(void* parameters, int index)
{
	PF_ComputePressureData* data = static_cast<PF_ComputePressureData*>(parameters);
	
	btFluidSphSolverDefault::calculateSumsInCellSymmetric(data->m_globalParameters, data->m_gridCellGroup[index], 
														data->m_grid, data->m_particles, data->m_sphData);
}

struct PF_ComputeForceData
{
	const btFluidSphParametersGlobal& m_globalParameters; 
	const btScalar m_vterm;
	const btAlignedObjectArray<int>& m_gridCellGroup;
	const btFluidSortingGrid& m_grid;
	btFluidParticles& m_particles;
	btFluidSphSolverDefault::SphParticles& m_sphData;
	
	PF_ComputeForceData(const btFluidSphParametersGlobal& FG, const btScalar vterm, const btAlignedObjectArray<int>& gridCellGroup, 
						const btFluidSortingGrid& grid, btFluidParticles& particles, btFluidSphSolverDefault::SphParticles& sphData) 
	: m_globalParameters(FG), m_vterm(vterm),  m_gridCellGroup(gridCellGroup), 
	m_grid(grid), m_particles(particles), m_sphData(sphData) {}
};
inline void PF_ComputeForceFunction(void* parameters, int index)
{
	PF_ComputeForceData* data = static_cast<PF_ComputeForceData*>(parameters);
	
	btFluidSphSolverDefault::calculateForcesInCellSymmetric(data->m_globalParameters, data->m_vterm, data->m_gridCellGroup[index], 
															data->m_grid, data->m_particles, data->m_sphData);
}

///@brief Multithreaded implementation of btFluidSphSolverDefault.
class btFluidSphSolverMultithreaded : public btFluidSphSolverDefault
{
	btParallelFor m_parallelFor;
	
public:
	///Use a different string for uniqueName if creating multiple instances of btFluidSphSolverMultithreaded
	btFluidSphSolverMultithreaded(int numThreads, const char* uniqueName = "btSphSolver_threads")
	: m_parallelFor(uniqueName, numThreads) {}

	virtual void computeSumsInMultithreadingGroup(const btFluidSphParametersGlobal& FG, const btAlignedObjectArray<int>& multithreadingGroup,
												const btFluidSortingGrid& grid, btFluidParticles& particles, 
												btFluidSphSolverDefault::SphParticles& sphData)
	{
		PF_ComputePressureData PressureData(FG, multithreadingGroup, grid, particles, sphData);
		m_parallelFor.execute( PF_ComputePressureFunction, &PressureData, 0, multithreadingGroup.size() - 1, 1 );
	}
	virtual void computeForcesInMultithreadingGroup(const btFluidSphParametersGlobal& FG, const btScalar vterm, 
													const btAlignedObjectArray<int>& multithreadingGroup, const btFluidSortingGrid& grid, 
													btFluidParticles& particles, btFluidSphSolverDefault::SphParticles& sphData)
	{
		PF_ComputeForceData ForceData(FG, vterm, multithreadingGroup, grid, particles, sphData);
		m_parallelFor.execute( PF_ComputeForceFunction, &ForceData, 0, multithreadingGroup.size() - 1, 1 );
	}
};

#endif


