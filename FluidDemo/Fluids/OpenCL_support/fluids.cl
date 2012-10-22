/* fluids.cl
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

typedef float btScalar;
typedef float4 btVector3;

//Note that these are vector3 functions -- OpenCL functions are vector4 functions
inline btScalar btVector3_length2(btVector3 v) { return v.x*v.x + v.y*v.y + v.z*v.z; }
inline btScalar btVector3_dot(btVector3 a, btVector3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline btVector3 btVector3_normalize(btVector3 v)
{
	btScalar length2 = btVector3_length2(v);
	if( length2 != (btScalar)0.0 ) v /= sqrt(length2);
	
	return v;
}

//Defined in "FluidSortingGrid.cpp" -- FluidSortingGrid::findAdjacentGridCells()
#define INVALID_FIRST_INDEX -1
#define INVALID_LAST_INDEX -2

//Syncronize with 'class FluidNeighbors' in "FluidParticles.h"
#define MAX_NEIGHBORS 80
typedef struct
{
	int m_count;
	int m_particleIndicies[MAX_NEIGHBORS];
	btScalar m_distances[MAX_NEIGHBORS];
	
} FluidNeighbors;

//Syncronize with 'struct FluidParametersGlobal' in "FluidParameters.h"
typedef struct
{
	btVector3 m_planeGravity;
	btVector3 m_pointGravityPosition;
	btScalar m_pointGravity;
	btScalar m_timeStep;
	btScalar m_simulationScale;
	btScalar m_particleRadius;
	btScalar m_speedLimit;
	btScalar m_sphSmoothRadius;
	btScalar m_sphRadiusSquared;
	btScalar m_poly6KernCoeff;
	btScalar m_spikyKernGradCoeff;
	btScalar m_viscosityKernLapCoeff;
} FluidParametersGlobal;

//Syncronize with 'struct FluidParametersLocal' in "FluidParameters.h"
typedef struct
{
	btVector3 m_volumeMin;
	btVector3 m_volumeMax;
	btScalar m_viscosity;
	btScalar m_restDensity;
	btScalar m_particleMass;
	btScalar m_stiffness;
	btScalar m_boundaryStiff;
	btScalar m_boundaryDamp;
	btScalar m_particleDist;
} FluidParametersLocal;


//#define SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED	//Ensure that this is also #defined in "FluidSortingGrid.h"
#ifdef SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED	
	typedef unsigned long SortGridUint64;
	typedef SortGridUint64 SortGridValue;		//Range must contain SORT_GRID_INDEX_RANGE^3
	//typedef short int SortGridIndex;
	//#define SORT_GRID_INDEX_RANGE 65536		//2^( 8*sizeof(SortGridIndex) )
	typedef int SortGridIndex;
	#define SORT_GRID_INDEX_RANGE 2097152		//2^21
#else
	typedef unsigned int SortGridValue;			//Range must contain SORT_GRID_INDEX_RANGE^3
	typedef char SortGridIndex;
	#define SORT_GRID_INDEX_RANGE 256			//2^( 8*sizeof(SortGridIndex) )
#endif

//SortGridIndex_large is a signed type with range including all values in:
//	[-HALVED_SORT_GRID_INDEX_RANGE, (HALVED_SORT_GRID_INDEX_RANGE - 1) + HALVED_SORT_GRID_INDEX_RANGE]
//
//	e.g if SORT_GRID_INDEX_RANGE == 256, HALVED_SORT_GRID_INDEX_RANGE == 128
//	then SortGridIndex_large must contain [-128, 127 + 128] == [-128, 255].
//	It is used to convert from SortGridIndex(signed, small range), to SortGridValue(unsigned, large range).
typedef int SortGridIndex_large;	

#define HALVED_SORT_GRID_INDEX_RANGE SORT_GRID_INDEX_RANGE/2


typedef struct
{
	int m_firstIndex;
	int m_lastIndex;
	
} FluidGridIterator;

typedef struct 
{
	SortGridValue m_value;
	int m_index;
	
} ValueIndexPair;

typedef struct
{
	SortGridIndex x;		
	SortGridIndex y;
	SortGridIndex z;
	SortGridIndex padding;
	
} SortGridIndicies;

SortGridIndicies getSortGridIndicies(btScalar cellSize, btVector3 position)	//FluidSortingGrid::generateIndicies()
{
	SortGridIndicies result;
	
	btVector3 discretePosition = position / cellSize;
	result.x = (SortGridIndex)( (position.x >= 0.0f) ? discretePosition.x : floor(discretePosition.x) );
	result.y = (SortGridIndex)( (position.y >= 0.0f) ? discretePosition.y : floor(discretePosition.y) );
	result.z = (SortGridIndex)( (position.z >= 0.0f) ? discretePosition.z : floor(discretePosition.z) );
	
	return result;
}
SortGridValue getSortGridValue(SortGridIndicies quantizedPosition)	//SortGridIndicies::getValue()
{
	SortGridIndex_large signedX = (SortGridIndex_large)quantizedPosition.x + HALVED_SORT_GRID_INDEX_RANGE;
	SortGridIndex_large signedY = (SortGridIndex_large)quantizedPosition.y + HALVED_SORT_GRID_INDEX_RANGE;
	SortGridIndex_large signedZ = (SortGridIndex_large)quantizedPosition.z + HALVED_SORT_GRID_INDEX_RANGE;
	
	SortGridValue unsignedX = (SortGridValue)signedX;
	SortGridValue unsignedY = (SortGridValue)signedY * SORT_GRID_INDEX_RANGE;
	SortGridValue unsignedZ = (SortGridValue)signedZ * SORT_GRID_INDEX_RANGE * SORT_GRID_INDEX_RANGE;
	
	return unsignedX + unsignedY + unsignedZ;
}

__kernel void generateValueIndexPairs(btScalar cellSize, __global btVector3 *fluidPositions, __global ValueIndexPair *out_pairs)
{
	int index = get_global_id(0);
	
	ValueIndexPair result;
	result.m_index = index;
	result.m_value = getSortGridValue( getSortGridIndicies(cellSize, fluidPositions[index]) );
	
	out_pairs[index] = result;
}

__kernel void rearrangeParticleArrays(__global ValueIndexPair *sortedPairs, __global btVector3 *rearrange, __global btVector3 *temporary)
{
	int index = get_global_id(0);
	
	//
	int oldIndex = sortedPairs[index].m_index;
	int newIndex = index;
	
	temporary[newIndex] = rearrange[oldIndex];
}

__kernel void generateUniques(__global ValueIndexPair *sortedPairs, int numSortedPairs,
							  __global SortGridValue *out_activeCells, __global FluidGridIterator *out_cellContents,
							  __global int *out_numActiveCells )
{
	//Assuming that out_activeCells[] is large enough to contain
	//all active cells( out_activeCells.size() >= numSortedPairs ).

	//Iterate from sortedPairs[0] to sortedPairs[numSortedPairs-1],
	//adding unique SortGridValue(s) and FluidGridIterator(s) to 
	//out_activeCells and out_cellContents, respectively.
	
	if( get_global_id(0) == 0 )
	{
		int numActiveCells = 0;
		
		if( numSortedPairs ) 
		{
			out_activeCells[numActiveCells] = sortedPairs[0].m_value;
			out_cellContents[numActiveCells] = (FluidGridIterator){-1, -2};
			++numActiveCells;
			
			out_cellContents[0].m_firstIndex = 0;
			out_cellContents[0].m_lastIndex = 0;
			
			for(int i = 1; i < numSortedPairs; ++i)
			{
				if( sortedPairs[i].m_value != sortedPairs[i - 1].m_value )
				{
					out_activeCells[numActiveCells] = sortedPairs[i].m_value;
					out_cellContents[numActiveCells] = (FluidGridIterator){-1, -2};
					++numActiveCells;
			
					int lastIndex = numActiveCells - 1;
					out_cellContents[lastIndex].m_firstIndex = i;
					out_cellContents[lastIndex].m_lastIndex = i;
					
					//
					out_cellContents[lastIndex - 1].m_lastIndex = i - 1;
				}
			}
			
			int valuesLastIndex = numSortedPairs - 1;
			if( sortedPairs[valuesLastIndex].m_value == sortedPairs[valuesLastIndex - 1].m_value )
			{
				int uniqueLastIndex = numActiveCells - 1;
				out_cellContents[uniqueLastIndex].m_lastIndex = valuesLastIndex;
			}
		}
		
		*out_numActiveCells = numActiveCells;
	}
}


//#define GRID_CELL_SIZE_2R		//Ensure that this is also #defined in "FluidSortingGrid.h"
#ifdef GRID_CELL_SIZE_2R	
	#define FluidSortingGrid_NUM_FOUND_CELLS 8
#else
	#define FluidSortingGrid_NUM_FOUND_CELLS 27
#endif

inline __global FluidGridIterator* findCell(int numActiveCells, __global SortGridValue *cellValues, __global FluidGridIterator *cellContents,
										   SortGridValue value)
{
	//#define USE_LINEAR_SEARCH
	#ifdef USE_LINEAR_SEARCH

		for(int i = 0; i < numActiveCells; ++i)
			if( value == cellValues[i] ) return &cellContents[i];
		
	#else

		//From btAlignedObjectArray::findBinarySearch()
		//Assumes cellValues[] is sorted
		
		int first = 0;
		int last = numActiveCells - 1;
		
		while(first <= last) 
		{
			int mid = (first + last) / 2;
			if(value > cellValues[mid]) first = mid + 1;
			else if(value < cellValues[mid]) last = mid - 1;
			else return &cellContents[mid];
		}
		
	#endif

	return 0;
}

inline void findCells(int numActiveCells, __global SortGridValue *cellValues, __global FluidGridIterator *cellContents, btScalar cellSize,
					  btVector3 position, btScalar radius, FluidGridIterator *out_cells)
{
	SortGridIndicies cellIndicies[FluidSortingGrid_NUM_FOUND_CELLS];	//	may be allocated in global memory(slow)
	
#ifdef GRID_CELL_SIZE_2R	
	btVector3 sphereMin = position - radius;
	SortGridIndicies indicies = getSortGridIndicies(cellSize, sphereMin);
	
	for(int i = 0; i < 4; ++i) cellIndicies[i] = indicies;
	cellIndicies[1].x++;
	cellIndicies[2].y++;
	cellIndicies[3].x++;
	cellIndicies[3].y++;
	
	for(int i = 0; i < 4; ++i)
	{
		cellIndicies[i+4] = cellIndicies[i];
		cellIndicies[i+4].z++;
	}
#else
	SortGridIndicies indicies = getSortGridIndicies(cellSize, position);

	for(int i = 0; i < 9; ++i) cellIndicies[i] = indicies;
	cellIndicies[1].x++;
	cellIndicies[2].y++;
	cellIndicies[3].x++;
	cellIndicies[3].y++;
	
	cellIndicies[4].x--;
	cellIndicies[5].y--;
	cellIndicies[6].x--;
	cellIndicies[6].y--;
	
	cellIndicies[7].x++;
	cellIndicies[7].y--;
	
	cellIndicies[8].x--;
	cellIndicies[8].y++;
	
	for(int i = 0; i < 9; ++i)
	{
		cellIndicies[i+9] = cellIndicies[i];
		cellIndicies[i+9].z++;
	}
	for(int i = 0; i < 9; ++i)
	{
		cellIndicies[i+18] = cellIndicies[i];
		cellIndicies[i+18].z--;
	}
#endif

	for(int i = 0; i < FluidSortingGrid_NUM_FOUND_CELLS; ++i)
	{
		__global FluidGridIterator *cell = findCell( numActiveCells, cellValues, cellContents, getSortGridValue(cellIndicies[i]) );
		
		out_cells[i].m_firstIndex = (cell) ? cell->m_firstIndex : INVALID_FIRST_INDEX;
		out_cells[i].m_lastIndex = (cell) ? cell->m_lastIndex : INVALID_LAST_INDEX;
	}
}

//
__kernel void sphComputePressure(__global FluidParametersGlobal *FG,  __global FluidParametersLocal *FL,
								  __global btVector3 *fluidPosition, __global btScalar *fluidPressure, 
								  __global btScalar *fluidInvDensity,  __global FluidNeighbors *fluidNeighbors,
								  __global int *numActiveCells, __global SortGridValue *cellValues, 
								  __global FluidGridIterator *cellContents, btScalar cellSize)
{
	btScalar searchRadius = FG->m_sphSmoothRadius / FG->m_simulationScale;

	int i = get_global_id(0);
	
	btScalar sum = 0.0f;	
	int neighborCount = 0;	//m_neighborTable[i].clear();
	
	FluidGridIterator foundCells[FluidSortingGrid_NUM_FOUND_CELLS];	//	may be allocated in global memory(slow)
	findCells(*numActiveCells, cellValues, cellContents, cellSize, fluidPosition[i], searchRadius, foundCells);
	
	for(int cell = 0; cell < FluidSortingGrid_NUM_FOUND_CELLS; ++cell) 
	{
		FluidGridIterator foundCell = foundCells[cell];
		
		for(int n = foundCell.m_firstIndex; n <= foundCell.m_lastIndex; ++n)
		{		
			if(i == n) continue; 
			
			btVector3 distance = (fluidPosition[i] - fluidPosition[n]) * FG->m_simulationScale;	//Simulation scale distance
			btScalar distanceSquared = btVector3_length2(distance);
			
			if(FG->m_sphRadiusSquared > distanceSquared) 
			{
				btScalar c = FG->m_sphRadiusSquared - distanceSquared;
				sum += c * c * c;
				
				if(neighborCount < MAX_NEIGHBORS)	//if( !m_neighborTable[i].isFilled() ) 
				{	
					//m_neighborTable[i].addNeighbor( n, sqrt(distanceSquared) );
					fluidNeighbors[i].m_particleIndicies[neighborCount] = n;
					fluidNeighbors[i].m_distances[neighborCount] = sqrt(distanceSquared);
					++neighborCount;
				}
				else break;
			}
		}
	}
	
	btScalar density = sum * FL->m_particleMass * FG->m_poly6KernCoeff;
	fluidPressure[i] = (density - FL->m_restDensity) * FL->m_stiffness;
	fluidInvDensity[i] = 1.0f / density;
	
	fluidNeighbors[i].m_count = neighborCount;
}


__kernel void sphComputeForce(__global FluidParametersGlobal *FG, __global FluidParametersLocal *FL,
							   __global btVector3 *fluidPosition, __global btVector3 *fluidVelEval, 
							   __global btVector3 *fluidSphForce, __global btScalar *fluidPressure, 
							   __global btScalar *fluidInvDensity, __global FluidNeighbors *fluidNeighbors)
{
	btScalar vterm = FG->m_viscosityKernLapCoeff * FL->m_viscosity;
	
	int i = get_global_id(0);
	
	btVector3 force = {0.0f, 0.0f, 0.0f, 0.0f};
	for(int j = 0; j < fluidNeighbors[i].m_count; ++j) 
	{
		int n = fluidNeighbors[i].m_particleIndicies[j];
	
		btVector3 distance = (fluidPosition[i] - fluidPosition[n]) * FG->m_simulationScale;	//Simulation scale distance
		
		btScalar c = FG->m_sphSmoothRadius - fluidNeighbors[i].m_distances[j];
		btScalar pterm = -0.5f * c * FG->m_spikyKernGradCoeff * ( fluidPressure[i] + fluidPressure[n] ) / fluidNeighbors[i].m_distances[j];
		btScalar dterm = c * fluidInvDensity[i] * fluidInvDensity[n];
		
		//force += (distance * pterm + (fluidVelEval[n] - fluidVelEval[i]) * vterm) * dterm;
		force.x += ( pterm * distance.x + vterm * (fluidVelEval[n].x - fluidVelEval[i].x) ) * dterm;
		force.y += ( pterm * distance.y + vterm * (fluidVelEval[n].y - fluidVelEval[i].y) ) * dterm;
		force.z += ( pterm * distance.z + vterm * (fluidVelEval[n].z - fluidVelEval[i].z) ) * dterm;
	}
	
	fluidSphForce[i] = force;
}


