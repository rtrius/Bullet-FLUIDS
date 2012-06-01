/** fluids.cl
	Fluids v.2 OpenCL Port

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

#ifdef cl_amd_printf
	#pragma OPENCL EXTENSION cl_amd_printf : enable
#endif

//	replace with float4?
typedef struct { float x, y, z, w; } __attribute__((aligned(16))) btVector3;

inline float btVector3_dot(btVector3 a, btVector3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline btVector3 btVector3_normalize(btVector3 v)
{ 
	float length = sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
	
	if(length != 0.f)
	{
		v.x /= length;
		v.y /= length;
		v.z /= length;
	}
	
	return v;
}

//Syncronize with 'struct Fluid' in "fluid.h"
#define INVALID_PARTICLE_INDEX -1
typedef struct
{
	btVector3		pos;
	btVector3		vel;			
	btVector3		vel_eval;
	btVector3		sph_force;
	btVector3		externalAcceleration;	
	btVector3		prev_pos;			//CCD_TEST
	float			pressure;
	float			density;	
	int				nextFluidIndex;
	
} Fluid;

//Syncronize with 'class Neighbors' in "fluid.h"
#define MAX_NEIGHBORS 80
typedef struct
{
	unsigned short m_count;
	unsigned short m_particleIndicies[MAX_NEIGHBORS];
	float m_distances[MAX_NEIGHBORS];
	
} Neighbors;

//Syncronize with 'struct FluidParameters_float' in "fluid.h"
typedef struct
{
	btVector3 m_volumeMin;
	btVector3 m_volumeMax;
	btVector3 m_planeGravity;
	btVector3 m_pointGravityPosition;
	float m_pointGravity;
	float sph_simscale;
	float sph_visc;
	float sph_restdensity;
	float sph_pmass;
	float sph_pradius;
	float sph_pdist;
	float sph_smoothradius;
	float sph_intstiff;
	float sph_extstiff;
	float sph_extdamp;
	float sph_limit;
	float m_timeStep;
	float m_R2;
	float m_Poly6Kern;
	float m_LapKern;
	float m_SpikyKern;
	
} FluidParameters_float;

//Syncronize with 'struct GridParameters' in "grid.h"
typedef struct
{
	btVector3 m_min;
	btVector3 m_max;
	float m_gridCellSize;
	int m_resolutionX;
	int m_resolutionY;
	int m_resolutionZ;
	int m_numCells;

} GridParameters;

////////////////////////////////////////////////////////////////////////////////
/// class FluidSystem
////////////////////////////////////////////////////////////////////////////////
__kernel void grid_insertParticles(__global Fluid *fluids, __global GridParameters *gridParams, 
								   __global volatile int *gridCells, __global volatile int *gridCellsNumFluids)	
{
	//Current implementation assumes that all values in
	//gridCells[] are set to INVALID_PARTICLE_INDEX,
	//and all values in gridCellsNumFluids[] are set to 0
	//before this function is called.

	__global GridParameters *GP = gridParams;
	
	//
	int particleIndex = get_global_id(0);
	__global Fluid *f = &fluids[particleIndex];
	
	//Reset particles
	f->nextFluidIndex = INVALID_PARTICLE_INDEX;
	
	//Load into grid
	int index_x = convert_int( (f->pos.x - GP->m_min.x) / GP->m_gridCellSize );
	int index_y = convert_int( (f->pos.y - GP->m_min.y) / GP->m_gridCellSize );
	int index_z = convert_int( (f->pos.z - GP->m_min.z) / GP->m_gridCellSize );
	
	int cellIndex = (index_z*GP->m_resolutionY + index_y)*GP->m_resolutionX + index_x;
	
	if(0 <= cellIndex && cellIndex < GP->m_numCells) 
	{
		//Add particle to linked list
		//Equivalent to:
		//		f->nextFluidIndex = grid->cells[cellIndex];
		//		grid->cells[cellIndex] = particleIndex;
		f->nextFluidIndex = atomic_xchg(&gridCells[cellIndex], particleIndex);
		
		//Equivalent to:
		//		grid->cells_num_fluids[cellIndex]++;
		atomic_inc(&gridCellsNumFluids[cellIndex]);
	}
}


//Grid::findCells()
#define RESULTS_PER_GRID_SEARCH 8
inline void findCells(__global GridParameters *gridParams, btVector3 position, float radius, int8 *out_cells)
{
	__global GridParameters *GP = gridParams;

	//Store a 2x2x2 grid cell query result in m_findCellsResult,
	//where m_findCellsResult.m_indicies[0], the cell with the lowest index,
	//corresponds to the minimum point of the sphere's AABB
	
	//Determine the grid cell index at the minimum point of the particle's AABB
	int index_x = convert_int( (-radius + position.x - GP->m_min.x) / GP->m_gridCellSize );
	int index_y = convert_int( (-radius + position.y - GP->m_min.y) / GP->m_gridCellSize );
	int index_z = convert_int( (-radius + position.z - GP->m_min.z) / GP->m_gridCellSize );
	
	//Clamp index to grid bounds
	if(index_x < 0) index_x = 0;
	if(index_y < 0) index_y = 0;
	if(index_z < 0) index_z = 0;
	
		//Since a 2x2x2 volume is accessed, subtract 2 from the upper index bounds
			//Subtract 1 as a 2x2x2 volume is accessed, and the index we want is the 'min' index
			//Subtract 1 again as indicies start from 0 (GP->m_resolutionX/Y/Z is out of bounds)
	if(index_x >= GP->m_resolutionX - 2) index_x = GP->m_resolutionX - 2;
	if(index_y >= GP->m_resolutionY - 2) index_y = GP->m_resolutionY - 2;
	if(index_z >= GP->m_resolutionZ - 2) index_z = GP->m_resolutionZ - 2;
	
	//Load indicies
	const int stride_x = 1;
	const int stride_y = GP->m_resolutionX;
	const int stride_z = GP->m_resolutionX*GP->m_resolutionY;
	
	(*out_cells).s0 = (index_z * GP->m_resolutionY + index_y) * GP->m_resolutionX + index_x ;
	(*out_cells).s1 = (*out_cells).s0 + stride_x;
	(*out_cells).s2 = (*out_cells).s0 + stride_y;
	(*out_cells).s3 = (*out_cells).s0 + stride_y + stride_x;

	(*out_cells).s4 = (*out_cells).s0 + stride_z;
	(*out_cells).s5 = (*out_cells).s1 + stride_z;
	(*out_cells).s6 = (*out_cells).s2 + stride_z;
	(*out_cells).s7 = (*out_cells).s3 + stride_z;
}

__kernel void sph_computePressure(__global FluidParameters_float *fluidParams, __global Fluid *fluids, __global Neighbors *neighbors,
								  __global GridParameters *gridParams,  __global int *gridCells)
{	
	__global FluidParameters_float *FP = fluidParams;
	
	float searchRadius = FP->sph_smoothradius / FP->sph_simscale;


	int i = get_global_id(0);
	__global Fluid* p = &fluids[i];
	
	float sum = 0.0f;	
	neighbors[i].m_count = 0;	//m_neighborTable[i].clear();

	int8 grid_query_result;
	findCells(gridParams, p->pos, searchRadius, &grid_query_result);
	
	int* query_result = (int*) &grid_query_result;
	for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; ++cell) 
	{
		int pndx = gridCells[ query_result[cell] ];				
		while(pndx != INVALID_PARTICLE_INDEX) 
		{
			__global Fluid* pcurr = &fluids[pndx];					
			if(pcurr == p) 
			{
				pndx = pcurr->nextFluidIndex; 
				continue; 
			}
			
			float dx = (p->pos.x - pcurr->pos.x) * FP->sph_simscale;		// dist in cm
			float dy = (p->pos.y - pcurr->pos.y) * FP->sph_simscale;
			float dz = (p->pos.z - pcurr->pos.z) * FP->sph_simscale;
			float dsq = dx*dx + dy*dy + dz*dz;
			if(FP->m_R2 > dsq) 
			{
				float c = FP->m_R2 - dsq;
				sum += c * c * c;
				
				if(neighbors[i].m_count < MAX_NEIGHBORS)	//if( !m_neighborTable[i].isFilled() ) 
				{	
					//m_neighborTable[i].addNeighbor( pndx, sqrt(dsq) );
					neighbors[i].m_particleIndicies[ neighbors[i].m_count ] = pndx;
					neighbors[i].m_distances[ neighbors[i].m_count ] = sqrt(dsq);
					++neighbors[i].m_count;
				}
			}
			
			pndx = pcurr->nextFluidIndex;
		}
	}
	
	p->density = sum * FP->sph_pmass * FP->m_Poly6Kern;	
	p->pressure = (p->density - FP->sph_restdensity) * FP->sph_intstiff;		
	p->density = 1.0f / p->density;		
}

__kernel void sph_computeForce(__global FluidParameters_float *fluidParams, __global Fluid *fluids, __global Neighbors *neighbors)
{
	__global FluidParameters_float *FP = fluidParams;

	float vterm = FP->m_LapKern * FP->sph_visc;
	
	int i = get_global_id(0);
	
	__global Fluid *p = &fluids[i];
	
	btVector3 force = {0, 0, 0};
	for(int j = 0; j < neighbors[i].m_count; ++j) 
	{
		__global Fluid *pcurr = &fluids[ neighbors[i].m_particleIndicies[j] ];
		float dx = (p->pos.x - pcurr->pos.x) * FP->sph_simscale;		// dist in cm
		float dy = (p->pos.y - pcurr->pos.y) * FP->sph_simscale;
		float dz = (p->pos.z - pcurr->pos.z) * FP->sph_simscale;
		float c = FP->sph_smoothradius - neighbors[i].m_distances[j];
		float pterm = -0.5f * c * FP->m_SpikyKern * ( p->pressure + pcurr->pressure) / neighbors[i].m_distances[j];
		float dterm = c * p->density * pcurr->density;
		force.x += ( pterm * dx + vterm * (pcurr->vel_eval.x - p->vel_eval.x) ) * dterm;
		force.y += ( pterm * dy + vterm * (pcurr->vel_eval.y - p->vel_eval.y) ) * dterm;
		force.z += ( pterm * dz + vterm * (pcurr->vel_eval.z - p->vel_eval.z) ) * dterm;
	}
	
	p->sph_force = force;
}

inline void resolveAabbCollision(float stiff, float damp, btVector3 vel_eval,
							 	 btVector3 *acceleration, btVector3 normal, float depthOfPenetration)
{
	const float COLLISION_EPSILON = 0.00001f;
	
	if(depthOfPenetration > COLLISION_EPSILON)
	{
		float adj = stiff * depthOfPenetration - damp * btVector3_dot(normal, vel_eval);
		acceleration->x += adj * normal.x; 
		acceleration->y += adj * normal.y; 
		acceleration->z += adj * normal.z;					
	}
}
__kernel void advance(__global FluidParameters_float *fluidParams, __global Fluid *fluids)
{
	__global FluidParameters_float *FP = fluidParams;
	
	float SL = FP->sph_limit;
	float SL2 = SL*SL;
	
	float stiff = FP->sph_extstiff;
	float damp = FP->sph_extdamp;
	float radius = FP->sph_pradius;
	float R2 = 2.0f * radius;
	float ss = FP->sph_simscale;
	
	btVector3 min = FP->m_volumeMin;
	btVector3 max = FP->m_volumeMax;
	
	bool planeGravityEnabled = ( FP->m_planeGravity.x != 0.0f 
							  || FP->m_planeGravity.y != 0.0f 
							  || FP->m_planeGravity.z != 0.0f );
	
	int i = get_global_id(0);
	__global Fluid *p = &fluids[i];

	//CCD_TEST
	p->prev_pos = p->pos;
	
	//Compute Acceleration		
	btVector3 accel = p->sph_force;
	accel.x *= FP->sph_pmass;
	accel.y *= FP->sph_pmass;
	accel.z *= FP->sph_pmass;

	//Velocity limiting 
	float speed = accel.x*accel.x + accel.y*accel.y + accel.z*accel.z;
	if(speed > SL2) 
	{
		float scaled_speed = SL / sqrt(speed);
		accel.x *= scaled_speed;
		accel.y *= scaled_speed;
		accel.z *= scaled_speed;
	}

	//Apply acceleration to keep particles in the FluidSystem's AABB
	resolveAabbCollision( stiff, damp, p->vel_eval, &accel, (btVector3){1.0f, 0.0f, 0.0f}, R2 - (p->pos.x - min.x)*ss );
	resolveAabbCollision( stiff, damp, p->vel_eval, &accel, (btVector3){-1.0f, 0.0f, 0.0f}, R2 - (max.x - p->pos.x)*ss );
	resolveAabbCollision( stiff, damp, p->vel_eval, &accel, (btVector3){0.0f, 1.0f, 0.0f}, R2 - (p->pos.y - min.y)*ss );
	resolveAabbCollision( stiff, damp, p->vel_eval, &accel, (btVector3){0.0f, -1.0f, 0.0f}, R2 - (max.y - p->pos.y)*ss );
	resolveAabbCollision( stiff, damp, p->vel_eval, &accel, (btVector3){0.0f, 0.0f, 1.0f}, R2 - (p->pos.z - min.z)*ss );
	resolveAabbCollision( stiff, damp, p->vel_eval, &accel, (btVector3){0.0f, 0.0f, -1.0f}, R2 - (max.z - p->pos.z)*ss );
	
	//Plane gravity
	if(planeGravityEnabled)
	{
		accel.x += FP->m_planeGravity.x;
		accel.y += FP->m_planeGravity.y;
		accel.z += FP->m_planeGravity.z;
	}

	//Point gravity
	if(FP->m_pointGravity > 0.0f) 
	{
		btVector3 norm;
		norm.x = p->pos.x - FP->m_pointGravityPosition.x;
		norm.y = p->pos.y - FP->m_pointGravityPosition.y;
		norm.z = p->pos.z - FP->m_pointGravityPosition.z;
		norm = btVector3_normalize(norm);
		norm.x *= FP->m_pointGravity;
		norm.y *= FP->m_pointGravity;
		norm.z *= FP->m_pointGravity;
		accel.x -= norm.x;
		accel.y -= norm.y;
		accel.z -= norm.z;
	}
	
	//Apply external forces
	accel.x += p->externalAcceleration.x;
	accel.y += p->externalAcceleration.y;
	accel.z += p->externalAcceleration.z;
	p->externalAcceleration = (btVector3){0.0f, 0.0f, 0.0f};

	// Leapfrog Integration ----------------------------
	btVector3 vnext = accel;							
	vnext.x *= FP->m_timeStep;
	vnext.y *= FP->m_timeStep;
	vnext.z *= FP->m_timeStep;
	vnext.x += p->vel.x;						// v(t+1/2) = v(t-1/2) + a(t) dt
	vnext.y += p->vel.y;						// v(t+1/2) = v(t-1/2) + a(t) dt
	vnext.z += p->vel.z;						// v(t+1/2) = v(t-1/2) + a(t) dt
	p->vel_eval = p->vel;
	p->vel_eval.x += vnext.x;
	p->vel_eval.y += vnext.y;
	p->vel_eval.z += vnext.z;
	p->vel_eval.x *= 0.5f;						// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
	p->vel_eval.y *= 0.5f;						// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
	p->vel_eval.z *= 0.5f;						// v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5		used to compute forces later
	p->vel = vnext;
	vnext.x *= FP->m_timeStep/ss;
	vnext.y *= FP->m_timeStep/ss;
	vnext.z *= FP->m_timeStep/ss;
	p->pos.x += vnext.x;						// p(t+1) = p(t) + v(t+1/2) dt
	p->pos.y += vnext.y;						// p(t+1) = p(t) + v(t+1/2) dt
	p->pos.z += vnext.z;						// p(t+1) = p(t) + v(t+1/2) dt
}




////////////////////////////////////////////////////////////////////////////////
/// Marching Cubes
////////////////////////////////////////////////////////////////////////////////

//	draft; incomplete and untested
//fluid_rendering.h - getVertex()
inline btVector3 getVertex(const btVector3 min, const btVector3 cell_size, int index_x, int index_y, int index_z)
{
	btVector3 v;
	v.x = min.x + cell_size.x * (float)index_x; 
	v.y = min.y + cell_size.y * (float)index_y;
	v.z = min.z + cell_size.z * (float)index_z;
						
	return v;
}
//FluidSystem::getValue()
inline float getValue(__global Fluid *fluids, __global GridParameters *gridParams, __global int *gridCells, 
					  float gridCellSize, float x, float y, float z)
{
	float sum = 0.0f;
	
	const float searchRadius = gridCellSize * 0.5f;
	const float R2 = 1.8f*1.8f;
	//const float R2 = 0.8f*0.8f;		//	marching cubes rendering test

	btVector3 position; 
	position.x = x;
	position.y = y;
	position.z = z;
	
	int8 grid_query_result;
	findCells(gridParams, position, searchRadius, &grid_query_result);
	
	int* query_result = (int*) &grid_query_result;
	for(int cell = 0; cell < RESULTS_PER_GRID_SEARCH; ++cell) 
	{
		int pndx = gridCells[ query_result[cell] ];
		while(pndx != INVALID_PARTICLE_INDEX) 
		{		
			__global Fluid *pFluid = &fluids[pndx];
			float dx = x - pFluid->pos.x;
			float dy = y - pFluid->pos.y;
			float dz = z - pFluid->pos.z;
			float dsq = dx*dx+dy*dy+dz*dz;
				
			if(dsq < R2) sum += R2 / dsq;
			
			pndx = pFluid->nextFluidIndex;
		}	
	}
	
	return sum;	
}

__kernel void loadMarchingCubeVertex(const btVector3 min, const btVector3 cell_size,
									 __global Fluid *fluids,  __global GridParameters *gridParams, __global int *gridCells, float gridCellSize,
									 __global float *march_cells, int march_cells_per_edge)
{
	int index_x = get_global_id(0);
	int index_y = get_global_id(1);
	int index_z = get_global_id(2);
	
	btVector3 position = getVertex(min, cell_size, index_x, index_y, index_z);
	float value = getValue(fluids, gridParams, gridCells, gridCellSize, position.x, position.y, position.z);
	
	march_cells[index_x + index_y * march_cells_per_edge + index_z * march_cells_per_edge* march_cells_per_edge] = value;
}

