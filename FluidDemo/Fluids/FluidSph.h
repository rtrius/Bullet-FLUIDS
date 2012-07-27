/*
  FLUIDS v.1 - SPH Fluid Simulator for CPU and GPU
  Copyright (C) 2008. Rama Hoetzlein, http://www.rchoetzlein.com

  ZLib license
  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
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
#ifndef FLUID_SPH_H
#define FLUID_SPH_H

#include "grid.h"
#include "hashgrid.h"

#include "FluidParticles.h"
#include "FluidParameters.h"

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

#define FLUIDS_OPENCL_ENABLED
#ifdef FLUIDS_OPENCL_ENABLED
	#include "OpenCL_support/fluids_opencl_support.h"
#endif


class FluidSph
{
	FluidParametersLocal	m_localParameters;
	
	FluidGrid				*m_grid;
	
	FluidParticles 			m_particles;
	
	btAlignedObjectArray<int> m_removedFluidIndicies;

#ifdef FLUIDS_OPENCL_ENABLED
	Grid_OpenCL				m_grid_OpenCL;
	Fluid_OpenCL 			m_fluid_OpenCL;
#endif

public:
	FluidSph() : m_grid(0) {}

	void initialize(const FluidParametersGlobal &FG, int maxNumParticles, const btVector3 &volumeMin, const btVector3 &volumeMax);
		
	void stepSimulation(const FluidParametersGlobal &FG);
	
	void reset(int maxNumParticles);
	void clear();
	
	int addParticle(const btVector3 &position) { return m_particles.addParticle(position); }
	void setPosition(int index, const btVector3 &position) { m_particles.m_pos[index] = position; }
	void setVelocity(int index, const btVector3 &velocity)
	{
		m_particles.m_vel[index] = velocity;
		m_particles.m_vel_eval[index] = velocity;
	}
	int getNextIndex(int index) const { return  m_particles.m_nextFluidIndex[index]; }
	const btVector3& getPosition(int index) const { return m_particles.m_pos[index]; }
	const btVector3& getPrevPosition(int index) const { return m_particles.m_prev_pos[index]; }
	const btVector3& getVelocity(int index) const { return m_particles.m_vel[index]; }
	const btVector3& getEvalVelocity(int index) const { return m_particles.m_vel_eval[index]; }
	
	void applyAcceleration(int index, const btVector3 &acceleration) { m_particles.m_externalAcceleration[index] += acceleration; }
	
	int	numParticles() const	{ return m_particles.size(); }
	
	
	///Removal occurs on the next call to stepSimulation()
	void markFluidForRemoval(int index) { m_removedFluidIndicies.push_back(index); }
	
	btScalar getEmitterSpacing(const FluidParametersGlobal &FG) const { return m_localParameters.m_particleDist / FG.sph_simscale; }
	
	const FluidGrid* getGrid() const { return m_grid; }
	const btAlignedObjectArray<int>& getNextFluidIndicies() const { return m_particles.m_nextFluidIndex; }
	
	//Parameters
	const FluidParametersLocal& getLocalParameters() const { return m_localParameters; }
	void setLocalParameters(const FluidParametersLocal &FP) { m_localParameters = FP; }

	void setDefaultParameters() { m_localParameters.setDefaultParameters(); }
	
	//Metablobs	
	btScalar getValue(btScalar x, btScalar y, btScalar z) const;
	btVector3 getGradient(btScalar x, btScalar y, btScalar z) const;

#ifdef FLUIDS_OPENCL_ENABLED
	void writeToOpenCl(cl_context context, cl_command_queue commandQueue, bool transferAllData)
	{
		m_grid_OpenCL.writeToOpenCl( context, commandQueue, static_cast<Grid*>(m_grid) );
		m_fluid_OpenCL.writeToOpenCl(context, commandQueue, &m_localParameters, &m_particles, transferAllData);
	}
	void readFromOpenCl(cl_context context, cl_command_queue commandQueue, bool transferAllData)
	{
		m_grid_OpenCL.readFromOpenCl(context, commandQueue, static_cast<Grid*>(m_grid) );
		m_fluid_OpenCL.readFromOpenCl(context, commandQueue, &m_localParameters, &m_particles, transferAllData);
	}
	Grid_OpenCL* getGridOpenCl() { return &m_grid_OpenCL; }
	Fluid_OpenCL* getFluidOpenCl() { return &m_fluid_OpenCL; }
	
	void preOpenClStepSimulation() { m_grid->clear(); }
	void postOpenClStepSimulation() 
	{
		for(int i = 0; i < m_particles.m_externalAcceleration.size(); ++i) 
			m_particles.m_externalAcceleration[i].setValue(0,0,0); 
	}
#endif
	
private:
	void removeMarkedFluids();
	
	void grid_insertParticles();
	
	void sph_computePressureGrid(const FluidParametersGlobal &FG);			//O(kn) - spatial grid
	void sph_computePressureGridReduce(const FluidParametersGlobal &FG);
	
	void sph_computeForceGrid(const FluidParametersGlobal &FG);
	void sph_computeForceGridNC(const FluidParametersGlobal &FG);			//O(cn) - neighbor table
	
	void integrate(const FluidParametersGlobal &FG);
	
#ifdef USE_HASHGRID	
	void grid_insertParticles_HASHGRID()
	{
		BT_PROFILE("FluidSystem::grid_insertParticles_HASHGRID()");
		m_hashgrid.clear();
		m_hashgrid.insertParticles(&m_particles);
	}
	void sph_computePressureGrid_HASHGRID(const FluidParametersGlobal &FG);
#endif
};

struct FluidEmitter
{
	btVector3 m_position;

	btScalar m_velocity;
	
	btScalar m_yaw;
	btScalar m_pitch;
	
	btScalar m_yawSpread;
	btScalar m_pitchSpread;
	
	FluidEmitter() : m_position(0,0,0), m_yaw(0), m_pitch(0), 
					 m_velocity(0), m_yawSpread(0), m_pitchSpread(0) {}
	
	void emit(FluidSph *fluid, int numParticles, btScalar spacing);

	static void addVolume(FluidSph *fluid, const btVector3 &min, const btVector3 &max, btScalar spacing);
};

inline bool isInsideAabb(const btVector3 &min, const btVector3 &max, const btVector3 &point)
{
	return (   min.x() <= point.x() && point.x() <= max.x()
			&& min.y() <= point.y() && point.y() <= max.y()
			&& min.z() <= point.z() && point.z() <= max.z() );
}
struct FluidAbsorber
{
	btVector3 m_min;
	btVector3 m_max;
	
	//int m_maxParticlesRemoved;
	//	add velocity limit / max particles removed, etc.?
	
	FluidAbsorber() : m_min(0,0,0), m_max(0,0,0) {}
	
	void absorb(FluidSph *fluid)
	{
		const FluidGrid *grid = fluid->getGrid();
		const bool isLinkedList = (grid->getGridType() == FT_LinkedList);
		
		btAlignedObjectArray<int> gridCellIndicies;
		grid->getGridCellIndiciesInAabb(m_min, m_max, &gridCellIndicies);
		
		for(int i = 0; i < gridCellIndicies.size(); ++i)
		{
			FluidGridIterator FI = grid->getGridCell( gridCellIndicies[i] );
			
			for( int n = FI.m_firstIndex; FluidGridIterator::isIndexValid(n, FI.m_lastIndex);
					 n = FluidGridIterator::getNextIndex(n, isLinkedList, fluid->getNextFluidIndicies()) )
			{
				if( isInsideAabb( m_min, m_max, fluid->getPosition(n) ) ) fluid->markFluidForRemoval(n);
			}
		}
	}
};

#endif


