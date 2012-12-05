/*
Bullet-FLUIDS 
Copyright (c) 2012 Jackson Lee

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. 
   If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/
#ifndef BT_FLUID_SOLVER_SPH_SHALLOW_H
#define BT_FLUID_SOLVER_SPH_SHALLOW_H
struct btSweSphParameters
{
	btScalar m_sphSmoothRadius;
	btScalar m_sphRadiusSquared;
	
	//2D kernels
	btScalar m_poly6KernCoeff2d;
	btScalar m_spikyKernGradCoeff2d;
	btScalar m_viscosityKernLapCoeff2d;
	
	void setSphInteractionRadius(btScalar radius)
	{
		m_sphSmoothRadius = radius;
		m_sphRadiusSquared = radius*radius;
		
		m_poly6KernCoeff2d = btScalar(4.0) / ( SIMD_PI * btPow(m_sphSmoothRadius, 8) );
		m_spikyKernGradCoeff2d = btScalar(-30.0) / ( SIMD_PI * btPow(m_sphSmoothRadius, 5) );
		m_viscosityKernLapCoeff2d = btScalar(40.0) / ( SIMD_PI * btPow(m_sphSmoothRadius, 5) );
	}
};

inline void resolveAabbCollision2(btScalar stiff, btScalar damp, const btVector3& vel_eval,
								 btVector3 *acceleration, const btVector3& normal, btScalar penetrationDepth)
{
	const btScalar COLLISION_EPSILON = btScalar(0.00001);

	if(penetrationDepth > COLLISION_EPSILON)
	{
		btScalar adj = stiff * penetrationDepth - damp * normal.dot(vel_eval);
		
		btVector3 collisionAcceleration = normal;
		collisionAcceleration *= adj;	

		*acceleration += collisionAcceleration;
	}
}

///Highly experimental; do not use; solves the Shallow Water Equations(SWE) using SPH
///@remarks
///A 2D SPH simulation is performed, with density interpreted as height.
class btFluidSolverSphShallow : public btFluidSolver
{
public:
	virtual void updateGridAndCalculateSphForces(const btFluidParametersGlobal& FG, btFluidSph** fluids, int numFluids)
	{
		BT_PROFILE("btFluidSolverSphShallow::updateGridAndCalculateSphForces()");
		
		btSweSphParameters SP;
		SP.setSphInteractionRadius(FG.m_sphSmoothRadius);
		
		for(int i = 0; i < numFluids; ++i)
		{
			//Hack - set height to 0.f so that the 2d particles will be aligned in the 3d grid
			//Particles would incorrectly be excluded from interaction otherwise
			btFluidParticles& particles = fluids[i]->internalGetParticles();
			for(int n = 0; n < particles.size(); ++n) particles.m_pos[n].setY(0.0);
			
			fluids[i]->insertParticlesIntoGrid();
			
			sphComputePressure( FG, SP, fluids[i] );
			
			sphComputeForce( FG, SP, fluids[i] );
			
			integrate( FG, fluids[i]->getLocalParameters(), fluids[i]->internalGetParticles() );
			
			//Update grid to correct y position, which is set in integrate()
			fluids[i]->insertParticlesIntoGrid();
		}
	}
	
protected:

	virtual void sphComputePressure(const btFluidParametersGlobal& FG, const btSweSphParameters& SP, btFluidSph* fluid)
	{
		BT_PROFILE("sphComputePressure()");
		
		const btFluidParametersLocal& FL = fluid->getLocalParameters();
		const btFluidSortingGrid& grid = fluid->getGrid();
		btFluidParticles& particles = fluid->internalGetParticles();
		
		for(int i = 0; i < fluid->numParticles(); ++i)
		{
			//btScalar sum = 0.0;	
			btScalar sum = SP.m_sphRadiusSquared*SP.m_sphRadiusSquared*SP.m_sphRadiusSquared;	
			particles.m_neighborTable[i].clear();

			btFluidSortingGrid::FoundCells foundCells;
			grid.findCells( grid.getDiscretePosition(particles.m_pos[i]), foundCells );
			
			for(int cell = 0; cell < btFluidSortingGrid::NUM_FOUND_CELLS; cell++) 
			{
				btFluidGridIterator& FI = foundCells.m_iterators[cell];
				
				for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
				{
					if(i == n) continue;
					
					//Simulation-scale distance
					btScalar difference_x = ( particles.m_pos[i].x() - particles.m_pos[n].x() ) * FG.m_simulationScale;	
					btScalar difference_z = ( particles.m_pos[i].z() - particles.m_pos[n].z() ) * FG.m_simulationScale;	
					btScalar distanceSquared = difference_x*difference_x + difference_z*difference_z;
					
					if(SP.m_sphRadiusSquared > distanceSquared) 
					{
						btScalar c = SP.m_sphRadiusSquared - distanceSquared;
						sum += c * c * c;
						
						if( !particles.m_neighborTable[i].isFilled() ) particles.m_neighborTable[i].addNeighbor( n, btSqrt(distanceSquared) );
						else break;
					}
				}
			}
			
			//btScalar density = sum * FL.m_particleMass * SP.m_poly6KernCoeff2d;
			//particles.m_invDensity[i] = density;	//Store density in 'inv density'(inaccurate naming)
			
			const btScalar STIFFNESS = 0.5;	//default 0.5
			const btScalar MASS_CORRECTION = ( btScalar(4.0)/btScalar(3.0) ) * FG.m_sphSmoothRadius;	//3D to 2D
			btScalar density = sum * FL.m_particleMass * SP.m_poly6KernCoeff2d;
			particles.m_pressure[i] = (density - FL.m_restDensity*MASS_CORRECTION) * STIFFNESS;
			particles.m_invDensity[i] = 1.0f / density;
		}
	}
	
	virtual void sphComputeForce(const btFluidParametersGlobal& FG, const btSweSphParameters& SP, btFluidSph* fluid)
	{
		BT_PROFILE("sphComputeForce()");
		
		const btFluidParametersLocal& FL = fluid->getLocalParameters();
		btFluidParticles& particles = fluid->internalGetParticles();
		
		const btScalar VISCOSITY = 0.2;		//default 0.2
		btScalar vterm = SP.m_viscosityKernLapCoeff2d * VISCOSITY;
		
		for(int i = 0; i < fluid->numParticles(); ++i)
		{
			btScalar force_x = 0.0f;
			btScalar force_z = 0.0f;
		
			for(int j = 0; j < particles.m_neighborTable[i].numNeighbors(); ++j)
			{
				int n = particles.m_neighborTable[i].getNeighborIndex(j);
				
				btScalar distance_x = ( particles.m_pos[i].x() - particles.m_pos[n].x() ) * FG.m_simulationScale;
				btScalar distance_z = ( particles.m_pos[i].z() - particles.m_pos[n].z() ) * FG.m_simulationScale;
				
				btScalar c = SP.m_sphSmoothRadius - particles.m_neighborTable[i].getDistance(j);
				btScalar pterm = -0.5f * c * SP.m_spikyKernGradCoeff2d 
							 * ( particles.m_pressure[i] + particles.m_pressure[n]) / particles.m_neighborTable[i].getDistance(j);
				btScalar dterm = c * particles.m_invDensity[i] * particles.m_invDensity[n];

				force_x += (pterm * distance_x + vterm * (particles.m_vel_eval[n].x() - particles.m_vel_eval[i].x())) * dterm;
				force_z += (pterm * distance_z + vterm * (particles.m_vel_eval[n].z() - particles.m_vel_eval[i].z())) * dterm;
			}
			
			particles.m_sph_force[i] = btVector3(force_x, 0, force_z) * FL.m_particleMass;
		}
	}
	
	virtual void integrate(const btFluidParametersGlobal& FG, const btFluidParametersLocal& FL, btFluidParticles& particles)
	{
		BT_PROFILE("integrate()");
		
		//const btScalar speedLimit2d = FG.m_speedLimit * ( btScalar(2.0) / btScalar(3.0) );
	
		const btScalar SPEED_LIMIT = 1.11;
		const btScalar SPEED_LIMIT_SQUARED = SPEED_LIMIT*SPEED_LIMIT;
		btScalar R2 = 2.0f * FL.m_particleRadius * FG.m_simulationScale;
	
		const btScalar ss = FG.m_simulationScale;
	
		const btScalar stiff = FL.m_boundaryStiff;
		const btScalar damp = FL.m_boundaryDamp;
		
		const btVector3& min = FL.m_volumeMin;
		const btVector3& max = FL.m_volumeMax;
		for(int i = 0; i < particles.size(); ++i)
		{
			//Compute Acceleration		
			btVector3 accel = particles.m_sph_force[i];
			accel *= FL.m_particleMass;

			//Limit speed
			btScalar speedSquared = accel.x()*accel.x() + accel.z()*accel.z();
			if(speedSquared > SPEED_LIMIT_SQUARED) accel *= SPEED_LIMIT / btSqrt(speedSquared);

			//Apply acceleration to keep particles in the FluidSystem's AABB
			resolveAabbCollision2( stiff, damp, particles.m_vel_eval[i], &accel,
									btVector3( 1.0, 0.0, 0.0), R2 - ( particles.m_pos[i].x() - min.x() )*ss );
			resolveAabbCollision2( stiff, damp, particles.m_vel_eval[i], &accel, 
									btVector3(-1.0, 0.0, 0.0), R2 - ( max.x() - particles.m_pos[i].x() )*ss );
			resolveAabbCollision2( stiff, damp, particles.m_vel_eval[i], &accel, 
									btVector3(0.0, 0.0,  1.0), R2 - ( particles.m_pos[i].z() - min.z() )*ss );
			resolveAabbCollision2( stiff, damp, particles.m_vel_eval[i], &accel, 
									btVector3(0.0, 0.0, -1.0), R2 - ( max.z() - particles.m_pos[i].z() )*ss );
			
			//Plane gravity
			//accel += FL.m_gravity;

			//Apply external forces
			//accel += particles.m_accumulatedForce[i] / FL.m_particleMass;
			//particles.m_accumulatedForce[i].setValue(0, 0, 0);

			//Integrate velocity
			particles.m_vel[i].setY(0.0);
			particles.m_vel_eval[i].setY(0.0);
			
			btVector3 vnext = particles.m_vel[i] + accel * FG.m_timeStep;			//v(t+1/2) = v(t-1/2) + a(t) dt	
			particles.m_vel_eval[i] = (particles.m_vel[i] + vnext) * btScalar(0.5);	//v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
			particles.m_vel[i] = vnext;
			
			//Integrate position
			particles.m_pos[i] += particles.m_vel[i]*(FG.m_timeStep / ss);
		
			//Set height
			const btScalar HEIGHT_SCALE = btScalar(60.0);
			btScalar height = (1.0f / particles.m_invDensity[i]) / FL.m_restDensity;
			particles.m_pos[i].setY(height * HEIGHT_SCALE);
		}
	}
};

#endif
