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

#include <vector>
#include <algorithm>

#include "grid.h"
#include "fluid.h"

#include "LinearMath/btQuickProf.h"		//BT_PROFILE(name) macro
#include "LinearMath/btAlignedObjectArray.h"


inline bool descendingSort(const int &a, const int &b) { return !(a < b); }		//for std::sort()
class FluidSystem
{
	int m_maxParticles;

	btAlignedObjectArray<Fluid> m_particles;
	
	Grid 					m_grid;	
	
	//Neighbor table of m_particles[i] is at m_neighborTable[i]
	//m_particles.size() == m_neighborTable.size()
	std::vector<Neighbors>  m_neighborTable;	
	
	FluidParameters			m_parameters;
	bool					m_useOpenCL;
	
	std::vector<int> m_removedFluidIndicies;
	
public:
	FluidSystem() : m_useOpenCL(false) {}

	void initialize(int maxNumParticles, const btVector3 &volumeMin, const btVector3 &volumeMax);
		
	void stepSimulation();
	
	void reset(int maxNumParticles);
	/*
	void clear()
	{
		m_particles.resize(0);
		m_neighborTable.clear();
	}
	*/
	
	Fluid* addFluid(const btVector3 &position)	{ return &m_particles[ addPointReuse(position) ]; }
	Fluid* getFluid(int index)	{ return &m_particles[index]; }
	const Fluid& getFluid(int index) const { return m_particles[index]; }
	int	numParticles() const	{ return m_particles.size(); }
	
	//Fluid particles are removed on the next call to stepSimulation()
	void markFluidForRemoval(int index) { m_removedFluidIndicies.push_back(index); }
	
	
	float getEmitterSpacing() const { return m_parameters.sph_pdist / m_parameters.sph_simscale; }
	
	const Grid& getGrid() const { return m_grid; }
	
	//Metablobs	
	float getValue(float x, float y, float z) const;
	btVector3 getGradient(float x, float y, float z) const;
	
	//Parameters
	const FluidParameters& getParameters() const { return m_parameters; }
	void setParameters(const FluidParameters &FP) { m_parameters = FP; }
	
	void toggleOpenCL() { m_useOpenCL = !m_useOpenCL; }
	
private:
	void grid_insertParticles();
	void advance();
	
	//int addPoint();		
	int addPointReuse(const btVector3 &position);

	void removeMarkedFluids();
	void removeFluid(int index);	//Invalidates grid
	
	//Smoothed Particle Hydrodynamics
	void sph_setup();
	void sph_computeKernels();

	void sph_computePressureSlow();			//O(n^2)
	void sph_computePressureGrid();			//O(kn) - spatial grid
	
	void sph_computeForceSlow();			//O(n^2)
	void sph_computeForceGrid();			//O(kn) - spatial grid
	void sph_computeForceGridNC();			//O(cn) - neighbor table
};
	
struct FluidEmitter
{
	btVector3 m_position;

	float m_velocity;
	
	float m_yaw;
	float m_pitch;
	
	float m_yawSpread;
	float m_pitchSpread;
	
	FluidEmitter() : m_position(0,0,0), m_yaw(0), m_pitch(0), 
					 m_velocity(0), m_yawSpread(0), m_pitchSpread(0) {}
	
	void emit(FluidSystem *fluidSystem, int numParticles, float spacing);

	static void addVolume(FluidSystem *fluidSystem, const btVector3 &min, const btVector3 &max, float spacing);
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
					int currentIndex = G.getLastParticleIndex( (z*GP.m_resolutionY + y)*GP.m_resolutionX + x );
					while(currentIndex != INVALID_PARTICLE_INDEX)
					{
						Fluid *f = fluidSystem->getFluid(currentIndex);
					
						if( isInsideAabb(m_min, m_max, f->pos) )
						{
							fluidSystem->markFluidForRemoval(currentIndex);
							//++numParticlesRemoved;
						}
						currentIndex = f->nextFluidIndex;
					}
				}
				
		//return numParticlesRemoved;
	}
};


#endif


