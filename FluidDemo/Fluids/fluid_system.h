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

#ifndef DEF_FLUID_SYS
#define DEF_FLUID_SYS

#include "grid.h"

#include "fluid.h"

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro

///USE HASHGRID enables a 255^3 grid that only stores nonempty cells.
///grid_insertParticles_HASHGRID() and HashGrid::findCells()
///are about 10 times slower than the finite grid equivalent.
///CPU only; OpenCL branch continues to use the finite grid.
//#define USE_HASHGRID
#ifdef USE_HASHGRID
	#include "hashgrid.h"
#endif

class FluidSystem
{
	Fluids 					m_fluids;
	
	Grid 					m_grid;	
#ifdef USE_HASHGRID
	HashGrid				m_hashgrid;
#endif
	
	FluidParameters			m_parameters;
	bool					m_useOpenCL;
	bool					m_toggledLastFrameOpenCL;
	
	btAlignedObjectArray<int> m_removedFluidIndicies;

public:
	FluidSystem() : m_useOpenCL(false), m_toggledLastFrameOpenCL(false) {}

	void initialize(int maxNumParticles, const btVector3 &volumeMin, const btVector3 &volumeMax);
		
	void stepSimulation();
	
	void reset(int maxNumParticles);
	void clear();
	
	int addFluid(const btVector3 &position) { return m_fluids.addFluid(position); }
	void setVelocity(int index, const btVector3 &velocity)
	{
		m_fluids.m_vel[index] = velocity;
		m_fluids.m_vel_eval[index] = velocity;
	}
	int getNextIndex(int index) const { return  m_fluids.m_nextFluidIndex[index]; }
	const btVector3& getPosition(int index) const { return m_fluids.m_pos[index]; }
	const btVector3& getPrevPosition(int index) const { return m_fluids.m_prev_pos[index]; }
	const btVector3& getVelocity(int index) const { return m_fluids.m_vel[index]; }
	const btVector3& getEvalVelocity(int index) const { return m_fluids.m_vel_eval[index]; }
	
	void applyAcceleration(int index, const btVector3 &acceleration) { m_fluids.m_externalAcceleration[index] += acceleration; }
	
	int	numParticles() const	{ return m_fluids.size(); }
	
	//Fluid particles are removed on the next call to stepSimulation()
	void markFluidForRemoval(int index) { m_removedFluidIndicies.push_back(index); }
	
	btScalar getEmitterSpacing() const { return m_parameters.sph_pdist / m_parameters.sph_simscale; }
	
	const Grid& getGrid() const { return m_grid; }
	
	//Parameters
	const FluidParameters& getParameters() const { return m_parameters; }
	void setParameters(const FluidParameters &FP) { m_parameters = FP; }
	void setDefaultParameters();
	
	void toggleOpenCL() 
	{ 
		m_useOpenCL = !m_useOpenCL; 
		m_toggledLastFrameOpenCL = true;
	}
	
	//Metablobs	
	btScalar getValue(btScalar x, btScalar y, btScalar z) const;
	btVector3 getGradient(btScalar x, btScalar y, btScalar z) const;
	

private:
	void removeMarkedFluids();
	
	void grid_insertParticles();
	void advance();
	
	//Smoothed Particle Hydrodynamics
	void sph_computePressureSlow();			//O(n^2)
	void sph_computePressureGrid();			//O(kn) - spatial grid
	void sph_computePressureGridReduce();
	
	void sph_computeForceSlow();			//O(n^2)
	void sph_computeForceGrid();			//O(kn) - spatial grid
	void sph_computeForceGridNC();			//O(cn) - neighbor table


#ifdef USE_HASHGRID	
	void grid_insertParticles_HASHGRID()
	{
		BT_PROFILE("FluidSystem::grid_insertParticles_HASHGRID()");
		m_hashgrid.clear();
		m_hashgrid.insertParticles(&m_fluids);
	}
	void sph_computePressureGrid_HASHGRID();
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
	
	void emit(FluidSystem *fluidSystem, int numParticles, btScalar spacing);

	static void addVolume(FluidSystem *fluidSystem, const btVector3 &min, const btVector3 &max, btScalar spacing);
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
	
	void absorb(FluidSystem *fluidSystem)
	{
		const Grid &G = fluidSystem->getGrid();
		const GridParameters &GP = G.getParameters();
		
		int gridMinX, gridMinY, gridMinZ;
		int gridMaxX, gridMaxY, gridMaxZ;
		G.getIndicies(m_min, &gridMinX, &gridMinY, &gridMinZ);
		G.getIndicies(m_max, &gridMaxX, &gridMaxY, &gridMaxZ);
		
		//int numParticlesRemoved = 0;
		for(int z = gridMinZ; z <= gridMaxZ; ++z)
			for(int y = gridMinY; y <= gridMaxY; ++y)
				for(int x = gridMinX; x <= gridMaxX; ++x)
				{
					int lastIndex = G.getLastParticleIndex( (z*GP.m_resolutionY + y)*GP.m_resolutionX + x );
					for(int n = lastIndex; n != INVALID_PARTICLE_INDEX; n = fluidSystem->getNextIndex(n) )
					{
						if( isInsideAabb( m_min, m_max, fluidSystem->getPosition(n) ) )
						{
							fluidSystem->markFluidForRemoval(n);
							//++numParticlesRemoved;
						}
					}
				}
				
		//return numParticlesRemoved;
	}
};


#endif


