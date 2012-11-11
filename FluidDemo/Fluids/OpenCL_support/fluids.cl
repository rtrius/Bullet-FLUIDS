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

//Defined in "btFluidSortingGrid.h"
#define INVALID_FIRST_INDEX -1
#define INVALID_LAST_INDEX -2


//Syncronize with 'struct btFluidParametersGlobal' in "btFluidParameters.h"
typedef struct
{
	btScalar m_timeStep;
	btScalar m_simulationScale;
	btScalar m_speedLimit;
	btScalar m_sphSmoothRadius;
	btScalar m_sphRadiusSquared;
	btScalar m_poly6KernCoeff;
	btScalar m_spikyKernGradCoeff;
	btScalar m_viscosityKernLapCoeff;
} btFluidParametersGlobal;

//Syncronize with 'struct btFluidParametersLocal' in "btFluidParameters.h"
typedef struct
{
	btVector3 m_volumeMin;
	btVector3 m_volumeMax;
	btVector3 m_gravity;
	btScalar m_viscosity;
	btScalar m_restDensity;
	btScalar m_particleMass;
	btScalar m_stiffness;
	btScalar m_particleRadius;
	btScalar m_boundaryStiff;
	btScalar m_boundaryDamp;
	btScalar m_particleDist;
} btFluidParametersLocal;


//#define SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED	//Ensure that this is also #defined in "btFluidSortingGrid.h"
#ifdef SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED	
	typedef unsigned long btSortGridUint64;
	typedef btSortGridUint64 btSortGridValue;		//Range must contain SORT_GRID_INDEX_RANGE^3
	typedef int btSortGridIndex;
	#define SORT_GRID_INDEX_RANGE 2097152		//2^21
#else
	typedef unsigned int btSortGridValue;			//Range must contain SORT_GRID_INDEX_RANGE^3
	typedef char btSortGridIndex;
	#define SORT_GRID_INDEX_RANGE 256			//2^( 8*sizeof(btSortGridIndex) )
#endif

//btSortGridIndex_large is a signed type with range including all values in:
//	[-HALVED_SORT_GRID_INDEX_RANGE, (HALVED_SORT_GRID_INDEX_RANGE - 1) + HALVED_SORT_GRID_INDEX_RANGE]
//
//	e.g if SORT_GRID_INDEX_RANGE == 256, HALVED_SORT_GRID_INDEX_RANGE == 128
//	then btSortGridIndex_large must contain [-128, 127 + 128] == [-128, 255].
//	It is used to convert from btSortGridIndex(signed, small range), to btSortGridValue(unsigned, large range).
typedef int btSortGridIndex_large;	

#define HALVED_SORT_GRID_INDEX_RANGE SORT_GRID_INDEX_RANGE/2


typedef struct
{
	int m_firstIndex;
	int m_lastIndex;
	
} btFluidGridIterator;

typedef struct 
{
	btSortGridValue m_value;
	int m_index;
	
} btValueIndexPair;

typedef struct
{
	btSortGridIndex x;		
	btSortGridIndex y;
	btSortGridIndex z;
	btSortGridIndex padding;
	
} btSortGridIndicies;

btSortGridIndicies getDiscretePosition(btScalar cellSize, btVector3 position)	//btFluidSortingGrid::getDiscretePosition()
{
	btSortGridIndicies result;
	
	btVector3 discretePosition = position / cellSize;
	result.x = (btSortGridIndex)( (position.x >= 0.0f) ? discretePosition.x : floor(discretePosition.x) );
	result.y = (btSortGridIndex)( (position.y >= 0.0f) ? discretePosition.y : floor(discretePosition.y) );
	result.z = (btSortGridIndex)( (position.z >= 0.0f) ? discretePosition.z : floor(discretePosition.z) );
	
	return result;
}
btSortGridValue getSortGridValue(btSortGridIndicies quantizedPosition)	//btSortGridIndicies::getValue()
{
	btSortGridIndex_large signedX = (btSortGridIndex_large)quantizedPosition.x + HALVED_SORT_GRID_INDEX_RANGE;
	btSortGridIndex_large signedY = (btSortGridIndex_large)quantizedPosition.y + HALVED_SORT_GRID_INDEX_RANGE;
	btSortGridIndex_large signedZ = (btSortGridIndex_large)quantizedPosition.z + HALVED_SORT_GRID_INDEX_RANGE;
	
	btSortGridValue unsignedX = (btSortGridValue)signedX;
	btSortGridValue unsignedY = (btSortGridValue)signedY * SORT_GRID_INDEX_RANGE;
	btSortGridValue unsignedZ = (btSortGridValue)signedZ * SORT_GRID_INDEX_RANGE * SORT_GRID_INDEX_RANGE;
	
	return unsignedX + unsignedY + unsignedZ;
}

__kernel void generateValueIndexPairs(__global btVector3* fluidPositions, __global btValueIndexPair* out_pairs, btScalar cellSize)
{
	int index = get_global_id(0);
	
	btValueIndexPair result;
	result.m_index = index;
	result.m_value = getSortGridValue( getDiscretePosition(cellSize, fluidPositions[index]) );
	
	out_pairs[index] = result;
}

__kernel void rearrangeParticleArrays(__global btValueIndexPair* sortedPairs, __global btVector3* rearrange, __global btVector3* temporary)
{
	int index = get_global_id(0);
	
	//
	int oldIndex = sortedPairs[index].m_index;
	int newIndex = index;
	
	temporary[newIndex] = rearrange[oldIndex];
}

__kernel void generateUniques(__global btValueIndexPair* sortedPairs, 
							  __global btSortGridValue* out_activeCells, __global btFluidGridIterator* out_cellContents,
							  __global int* out_numActiveCells, int numSortedPairs )
{
	//Assuming that out_activeCells[] is large enough to contain
	//all active cells( out_activeCells.size() >= numSortedPairs ).

	//Iterate from sortedPairs[0] to sortedPairs[numSortedPairs-1],
	//adding unique btSortGridValue(s) and btFluidGridIterator(s) to 
	//out_activeCells and out_cellContents, respectively.
	
	if( get_global_id(0) == 0 )
	{
		int numActiveCells = 0;
		
		if( numSortedPairs ) 
		{
			out_activeCells[numActiveCells] = sortedPairs[0].m_value;
			out_cellContents[numActiveCells] = (btFluidGridIterator){INVALID_FIRST_INDEX, INVALID_LAST_INDEX};
			++numActiveCells;
			
			out_cellContents[0].m_firstIndex = 0;
			out_cellContents[0].m_lastIndex = 0;
			
			for(int i = 1; i < numSortedPairs; ++i)
			{
				if( sortedPairs[i].m_value != sortedPairs[i - 1].m_value )
				{
					out_activeCells[numActiveCells] = sortedPairs[i].m_value;
					out_cellContents[numActiveCells] = (btFluidGridIterator){INVALID_FIRST_INDEX, INVALID_LAST_INDEX};
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


//Note that this value differs from btFluidSortingGrid::NUM_FOUND_CELLS in "btFluidSortingGrid.h"
//
//Since the hash function used to determine the 'value' of particles is simply 
//(x + y*CELLS_PER_ROW + z*CELLS_PER_PLANE), adjacent cells have a value 
//that is 1 greater and lesser than the current cell. 
//This makes it possible to query 3 cells simultaneously(as a 3 cell bar extended along the x-axis) 
//by using a 'binary range search' in the range [current_cell_value-1, current_cell_value+1]. 
//Furthermore, as the 3 particle index ranges returned are also adjacent, it is also possible to 
//stitch them together to form a single index range.
#define btFluidSortingGrid_NUM_FOUND_CELLS 9

inline void binaryRangeSearch(int numActiveCells, __global btSortGridValue* cellValues,
							  btSortGridValue lowerValue, btSortGridValue upperValue, int* out_lowerIndex, int* out_upperIndex)
{
	int first = 0;
	int last = numActiveCells - 1;
	
	while(first <= last)
	{
		int mid = (first + last) / 2;
		if( lowerValue > cellValues[mid] )
		{
			first = mid + 1;
		}
		else if( upperValue < cellValues[mid] )
		{
			last = mid - 1;
		}
		else 
		{
			//At this point, (lowerValue <= cellValues[mid] <= upperValue)
			//Perform a linear search to find the lower and upper index range
		
			int lowerIndex = mid;
			int upperIndex = mid;
			while(lowerIndex-1 >= 0 && cellValues[lowerIndex-1] >= lowerValue) lowerIndex--;
			while(upperIndex+1 < numActiveCells && cellValues[upperIndex+1] <= upperValue) upperIndex++;
		
			*out_lowerIndex = lowerIndex;
			*out_upperIndex = upperIndex;
			return;
		}
	}

	*out_lowerIndex = numActiveCells;
	*out_upperIndex = numActiveCells;
}

inline void findCells(int numActiveCells, __global btSortGridValue* cellValues, __global btFluidGridIterator* cellContents, 
						btScalar cellSize, btVector3 position, btFluidGridIterator* out_cells)
{
	btSortGridIndicies cellIndicies[btFluidSortingGrid_NUM_FOUND_CELLS];	//	may be allocated in global memory(slow)
	
	btSortGridIndicies indicies = getDiscretePosition(cellSize, position);

	for(int i = 0; i < 9; ++i) cellIndicies[i] = indicies;
	cellIndicies[1].y++;
	cellIndicies[2].z++;
	cellIndicies[3].y++;
	cellIndicies[3].z++;
	
	cellIndicies[4].y--;
	cellIndicies[5].z--;
	cellIndicies[6].y--;
	cellIndicies[6].z--;
	
	cellIndicies[7].y++;
	cellIndicies[7].z--;
	
	cellIndicies[8].y--;
	cellIndicies[8].z++;
	
	for(int i = 0; i < btFluidSortingGrid_NUM_FOUND_CELLS; ++i) out_cells[i] = (btFluidGridIterator){INVALID_FIRST_INDEX, INVALID_LAST_INDEX};
	for(int i = 0; i < 9; ++i)
	{
		btSortGridIndicies lower = cellIndicies[i];
		lower.x--;
	
		btSortGridIndicies upper = cellIndicies[i];
		upper.x++;
		
		int lowerIndex, upperIndex;
		binaryRangeSearch(numActiveCells, cellValues, getSortGridValue(lower), getSortGridValue(upper), &lowerIndex, &upperIndex);
		
		if(lowerIndex != numActiveCells)
		{
			out_cells[i] = (btFluidGridIterator){cellContents[lowerIndex].m_firstIndex, cellContents[upperIndex].m_lastIndex};
		}
	}
}

//
#define MAX_NEIGHBORS 80
__kernel void sphComputePressure(__global btFluidParametersGlobal* FG,  __global btFluidParametersLocal* FL,
								  __global btVector3* fluidPosition, __global btScalar* fluidDensity,
								  __global int* numActiveCells, __global btSortGridValue* cellValues, 
								  __global btFluidGridIterator* cellContents, btScalar cellSize)
{
	int i = get_global_id(0);
	
	btScalar sum = 0.0f;
	int neighborCount = 0;
	
	btFluidGridIterator foundCells[btFluidSortingGrid_NUM_FOUND_CELLS];
	findCells(*numActiveCells, cellValues, cellContents, cellSize, fluidPosition[i], foundCells);
	
	for(int cell = 0; cell < btFluidSortingGrid_NUM_FOUND_CELLS; ++cell) 
	{
		btFluidGridIterator foundCell = foundCells[cell];
		
		for(int n = foundCell.m_firstIndex; n <= foundCell.m_lastIndex; ++n)
		{
			if(i == n) continue; 
			
			btVector3 delta = (fluidPosition[i] - fluidPosition[n]) * FG->m_simulationScale;	//Simulation scale distance
			btScalar distanceSquared = btVector3_length2(delta);
			
			if(FG->m_sphRadiusSquared > distanceSquared) 
			{
				btScalar c = FG->m_sphRadiusSquared - distanceSquared;
				sum += c * c * c;
				
				if(++neighborCount >= MAX_NEIGHBORS) break;
			}
		}
	}
	
	fluidDensity[i] = sum * FL->m_particleMass * FG->m_poly6KernCoeff;
}


__kernel void sphComputeForce(__global btFluidParametersGlobal* FG, __global btFluidParametersLocal* FL,
							   __global btVector3* fluidPosition, __global btVector3* fluidVelEval, 
							   __global btVector3* fluidSphForce, __global btScalar* fluidDensity,
							   __global int* numActiveCells, __global btSortGridValue* cellValues, 
							   __global btFluidGridIterator* cellContents, btScalar cellSize)
{
	btScalar vterm = FG->m_viscosityKernLapCoeff * FL->m_viscosity;
	
	int i = get_global_id(0);
	btScalar density_i = fluidDensity[i];
	btScalar invDensity_i = 1.0f / density_i;
	btScalar pressure_i = (density_i - FL->m_restDensity) * FL->m_stiffness;
	
	btVector3 force = {0.0f, 0.0f, 0.0f, 0.0f};
	int neighborCount = 0;
	
	btFluidGridIterator foundCells[btFluidSortingGrid_NUM_FOUND_CELLS];
	findCells(*numActiveCells, cellValues, cellContents, cellSize, fluidPosition[i], foundCells);
	
	for(int cell = 0; cell < btFluidSortingGrid_NUM_FOUND_CELLS; ++cell) 
	{
		btFluidGridIterator foundCell = foundCells[cell];
		
		for(int n = foundCell.m_firstIndex; n <= foundCell.m_lastIndex; ++n)
		{
			if(i == n) continue; 
			
			btVector3 delta = (fluidPosition[i] - fluidPosition[n]) * FG->m_simulationScale;	//Simulation scale distance
			btScalar distanceSquared = btVector3_length2(delta);
			
			if(FG->m_sphRadiusSquared > distanceSquared) 
			{
				btScalar density_n = fluidDensity[n];
				btScalar invDensity_n = 1.0f / density_n;
				btScalar pressure_n = (density_n - FL->m_restDensity) * FL->m_stiffness;
			
				btScalar distance = sqrt(distanceSquared);
				btScalar c = FG->m_sphSmoothRadius - distance;
				btScalar pterm = -0.5f * c * FG->m_spikyKernGradCoeff * (pressure_i + pressure_n) / distance;
				btScalar dterm = c * invDensity_i * invDensity_n;
				
				//force += (delta * pterm + (fluidVelEval[n] - fluidVelEval[i]) * vterm) * dterm;
				force.x += ( pterm * delta.x + vterm * (fluidVelEval[n].x - fluidVelEval[i].x) ) * dterm;
				force.y += ( pterm * delta.y + vterm * (fluidVelEval[n].y - fluidVelEval[i].y) ) * dterm;
				force.z += ( pterm * delta.z + vterm * (fluidVelEval[n].z - fluidVelEval[i].z) ) * dterm;
		
				if(++neighborCount >= MAX_NEIGHBORS) break;
			}
		}
	}
	
	fluidSphForce[i] = force * FL->m_particleMass;
}


