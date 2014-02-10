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

#ifdef cl_amd_printf
	#pragma OPENCL EXTENSION cl_amd_printf : enable
#endif

typedef float btScalar;
typedef float4 btVector3;
#define btMax max
#define btMin min


//Note that these are vector3 functions -- OpenCL functions are vector4 functions
inline btScalar btVector3_length2(btVector3 v) { return v.x*v.x + v.y*v.y + v.z*v.z; }
inline btScalar btVector3_dot(btVector3 a, btVector3 b) { return a.x*b.x + a.y*b.y + a.z*b.z; }
inline btVector3 btVector3_normalize(btVector3 v)
{
	btScalar length2 = btVector3_length2(v);
	if( length2 != (btScalar)0.0f ) v /= sqrt(length2);
	
	return v;
}
//Defined in btFluidSortingGrid.h
#define INVALID_FIRST_INDEX -1
#define INVALID_LAST_INDEX -2


//Syncronize with 'struct btFluidSphParametersGlobal' in btFluidSphParameters.h
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
} btFluidSphParametersGlobal;

//Syncronize with 'struct btFluidSphParametersLocal' in btFluidSphParameters.h
typedef struct
{
	btVector3 m_aabbBoundaryMin;
	btVector3 m_aabbBoundaryMax;
	int m_enableAabbBoundary;
	btVector3 m_gravity;
	btScalar m_viscosity;
	btScalar m_restDensity;
	btScalar m_sphParticleMass;
	btScalar m_stiffness;
	btScalar m_initialSum;
	btScalar m_surfaceTension;
	btScalar m_particleDist;
	btScalar m_particleRadius;
	btScalar m_particleMargin;
	btScalar m_particleMass;
	btScalar m_particleRadiusExpansion;
	btScalar m_boundaryStiff;
	btScalar m_boundaryDamp;
	btScalar m_boundaryFriction;
	btScalar m_boundaryRestitution;
	btScalar m_boundaryErp;
} btFluidSphParametersLocal;


#define BT_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT	//Ensure that this is also #defined in btFluidSortingGrid.h
#ifdef BT_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT	
	typedef unsigned long btFluidGridUint64;
	typedef btFluidGridUint64 btFluidGridCombinedPos;	//Range must contain BT_FLUID_GRID_COORD_RANGE^3
	#define BT_FLUID_GRID_COORD_RANGE 2097152		//2^21
	
	inline void splitCombinedPosition(btFluidGridUint64 resolutionX, btFluidGridUint64 resolutionY, 
										btFluidGridUint64 value, int* out_x, int* out_y, int* out_z)
	{
		btFluidGridUint64 cellsPerLine = resolutionX;
		btFluidGridUint64 cellsPerPlane = resolutionX * resolutionY;
		
		btFluidGridUint64 x = value % cellsPerLine;
		btFluidGridUint64 z = value / cellsPerPlane;
		btFluidGridUint64 y = (value - z*cellsPerPlane) / cellsPerLine;
		
		*out_x = (int)x;
		*out_z = (int)z;
		*out_y = (int)y;
	}
#else
	typedef unsigned int btFluidGridCombinedPos;		//Range must contain BT_FLUID_GRID_COORD_RANGE^3
	#define BT_FLUID_GRID_COORD_RANGE 1024			//2^10	
	
	inline void splitCombinedPosition(int resolutionX, int resolutionY, int value, int* out_x, int* out_y, int* out_z)
	{
		int x = value % resolutionX;
		int z = value / (resolutionX*resolutionY);
		int y = (value - z*resolutionX*resolutionY) / resolutionX;
		
		*out_x = (int)x;
		*out_z = (int)z;
		*out_y = (int)y;
	}
#endif

typedef int btFluidGridCoordinate;
#define BT_FLUID_GRID_COORD_RANGE_HALVED BT_FLUID_GRID_COORD_RANGE/2



typedef struct
{
	int m_firstIndex;
	int m_lastIndex;
	
} btFluidGridIterator;


//Since the hash function used to determine the 'value' of particles is simply 
//(x + y*CELLS_PER_ROW + z*CELLS_PER_PLANE), adjacent cells have a value 
//that is 1 greater and lesser than the current cell. 
//This makes it possible to query 3 cells simultaneously(as a 3 cell bar extended along the x-axis) 
//by using a 'binary range search' in the range [current_cell_value-1, current_cell_value+1]. 
//Furthermore, as the 3 particle index ranges returned are also adjacent, it is also possible to 
//stitch them together to form a single index range.
#define btFluidSortingGrid_NUM_FOUND_CELLS_GPU 9

typedef struct
{
	btFluidGridIterator m_iterators[btFluidSortingGrid_NUM_FOUND_CELLS_GPU];
	
} btFluidSortingGridFoundCellsGpu;		//btFluidSortingGrid::FoundCellsGpu in btFluidSortingGrid.h

typedef struct 
{
	btFluidGridCombinedPos m_value;
	int m_index;
	
} btFluidGridValueIndexPair;

typedef struct
{
	btFluidGridCoordinate x;		
	btFluidGridCoordinate y;
	btFluidGridCoordinate z;
	btFluidGridCoordinate padding;
	
} btFluidGridPosition;

btFluidGridPosition getDiscretePosition(btScalar cellSize, btVector3 position)	//btFluidSortingGrid::getDiscretePosition()
{
	btVector3 discretePosition = position / cellSize;
	
	btFluidGridPosition result;
	result.x = (btFluidGridCoordinate)( (position.x >= 0.0f) ? discretePosition.x : floor(discretePosition.x) );
	result.y = (btFluidGridCoordinate)( (position.y >= 0.0f) ? discretePosition.y : floor(discretePosition.y) );
	result.z = (btFluidGridCoordinate)( (position.z >= 0.0f) ? discretePosition.z : floor(discretePosition.z) );
	
	return result;
}
btFluidGridCombinedPos getCombinedPosition(btFluidGridPosition quantizedPosition)	//btFluidGridPosition::getCombinedPosition()
{
	btFluidGridCoordinate signedX = quantizedPosition.x + BT_FLUID_GRID_COORD_RANGE_HALVED;
	btFluidGridCoordinate signedY = quantizedPosition.y + BT_FLUID_GRID_COORD_RANGE_HALVED;
	btFluidGridCoordinate signedZ = quantizedPosition.z + BT_FLUID_GRID_COORD_RANGE_HALVED;
	
	btFluidGridCombinedPos unsignedX = (btFluidGridCombinedPos)signedX;
	btFluidGridCombinedPos unsignedY = (btFluidGridCombinedPos)signedY * BT_FLUID_GRID_COORD_RANGE;
	btFluidGridCombinedPos unsignedZ = (btFluidGridCombinedPos)signedZ * BT_FLUID_GRID_COORD_RANGE * BT_FLUID_GRID_COORD_RANGE;
	
	return unsignedX + unsignedY + unsignedZ;
}

__kernel void generateValueIndexPairs(__global btVector3* fluidPositions, __global btFluidGridValueIndexPair* out_pairs, 
										btScalar cellSize, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	btFluidGridValueIndexPair result;
	result.m_index = index;
	result.m_value = getCombinedPosition( getDiscretePosition(cellSize, fluidPositions[index]) );
	
	out_pairs[index] = result;
}

__kernel void rearrangeParticleArrays(__global btFluidGridValueIndexPair* sortedPairs, __global btVector3* rearrange, 
										__global btVector3* temporary, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	//
	int oldIndex = sortedPairs[index].m_index;
	int newIndex = index;
	
	temporary[newIndex] = rearrange[oldIndex];
}


__kernel void markUniques(__global btFluidGridValueIndexPair* valueIndexPairs, __global int* out_retainValueAtThisIndex, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	int lastValidIndex = numFluidParticles - 1;
	
	//Retain if the next particle has a different btFluidGridCombinedPos(is in another cell)
	int isRetained = (index < lastValidIndex) ? (valueIndexPairs[index].m_value != valueIndexPairs[index+1].m_value) : 1;
	
	out_retainValueAtThisIndex[index] = isRetained;
}
__kernel void storeUniques(__global btFluidGridValueIndexPair* valueIndexPairs, __global int* retainValue, __global int* scanResults, 
							__global btFluidGridCombinedPos* out_sortGridValues, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	if(retainValue[index])
	{
		int scannedIndex = scanResults[index];
		
		out_sortGridValues[scannedIndex] = valueIndexPairs[index].m_value;
	}
}
__kernel void setZero(__global int* array, int numUniques)
{
	int index = get_global_id(0);
	if(index >= numUniques) return;
	
	array[index] = 0;
}
inline int binarySearch(__global btFluidGridCombinedPos *sortGridValues, int sortGridValuesSize, btFluidGridCombinedPos value)
{
	//From btAlignedObjectArray::findBinarySearch()
	//Assumes sortGridValues[] is sorted
	
	int first = 0;
	int last = sortGridValuesSize - 1;
	
	while(first <= last) 
	{
		int mid = (first + last) / 2;
		if(value > sortGridValues[mid]) first = mid + 1;
		else if(value < sortGridValues[mid]) last = mid - 1;
		else return mid;
	}

	return sortGridValuesSize;
}
__kernel void countUniques(__global btFluidGridValueIndexPair* valueIndexPairs, __global btFluidGridCombinedPos* sortGridValues, 
							__global int* out_valuesCount, int numUniqueValues, int numFluidParticles)
{
	int index = get_global_id(0);
	if(index >= numFluidParticles) return;
	
	btFluidGridCombinedPos particleValue = valueIndexPairs[index].m_value;
	
	int countArrayIndex = binarySearch(sortGridValues, numUniqueValues, particleValue);
	
	//particleValue should exist in sortGridValues; this check is not necessary
	if(countArrayIndex != numUniqueValues) 
		atomic_inc( &out_valuesCount[countArrayIndex] );
}
__kernel void generateIndexRanges(__global int* scanResults, __global btFluidGridIterator* out_iterators, int numActiveCells, int numParticles)
{
	int index = get_global_id(0);
	if(index >= numActiveCells) return;
	
	int lowerIndex, upperIndex;
	if(index < numActiveCells-1)
	{
		lowerIndex = scanResults[index];
		upperIndex = scanResults[index+1] - 1;
	}
	else
	{
		lowerIndex = scanResults[index];
		upperIndex = numParticles - 1;
	}
	
	out_iterators[index] = (btFluidGridIterator){ lowerIndex, upperIndex };
}

__kernel void generateUniques(__global btFluidGridValueIndexPair* sortedPairs, 
							  __global btFluidGridCombinedPos* out_activeCells, __global btFluidGridIterator* out_cellContents,
							  __global int* out_numActiveCells, int numSortedPairs )
{
	//Assuming that out_activeCells[] is large enough to contain
	//all active cells( out_activeCells.size() >= numSortedPairs ).

	//Iterate from sortedPairs[0] to sortedPairs[numSortedPairs-1],
	//adding unique btFluidGridCombinedPos(s) and btFluidGridIterator(s) to 
	//out_activeCells and out_cellContents, respectively.
	
	if( get_global_id(0) == 0 )
	{
		int numActiveCells = 0;
		
		if( numSortedPairs ) 
		{
			//Crashes on compiling with Catalyst 13.1 if
			//(btFluidGridIterator){INVALID_FIRST_INDEX, INVALID_FIRST_INDEX} is used directly
			int invalidLowerIndex = INVALID_FIRST_INDEX;
			int invalidUpperIndex = INVALID_LAST_INDEX;
		
			out_activeCells[numActiveCells] = sortedPairs[0].m_value;
			//out_cellContents[numActiveCells] = (btFluidGridIterator){INVALID_FIRST_INDEX, INVALID_FIRST_INDEX};
			out_cellContents[numActiveCells] = (btFluidGridIterator){invalidLowerIndex, invalidUpperIndex};
			++numActiveCells;
			
			out_cellContents[0].m_firstIndex = 0;
			out_cellContents[0].m_lastIndex = 0;
			
			for(int i = 1; i < numSortedPairs; ++i)
			{
				if( sortedPairs[i].m_value != sortedPairs[i - 1].m_value )
				{
					out_activeCells[numActiveCells] = sortedPairs[i].m_value;
					//out_cellContents[numActiveCells] = (btFluidGridIterator){INVALID_FIRST_INDEX, INVALID_FIRST_INDEX};
					out_cellContents[numActiveCells] = (btFluidGridIterator){invalidLowerIndex, invalidUpperIndex};
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

inline void binaryRangeSearch(int numActiveCells, __global btFluidGridCombinedPos* cellValues,
							  btFluidGridCombinedPos lowerValue, btFluidGridCombinedPos upperValue, int* out_lowerIndex, int* out_upperIndex)
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

inline void findCellsFromGridPosition(int numActiveCells, __global btFluidGridCombinedPos* cellValues, __global btFluidGridIterator* cellContents, 
										btFluidGridPosition combinedPosition, btFluidGridIterator* out_cells)
{
	btFluidGridPosition cellIndicies[btFluidSortingGrid_NUM_FOUND_CELLS_GPU];	//	may be allocated in global memory(slow)
	
	btFluidGridPosition indicies = combinedPosition;

	for(int i = 0; i < btFluidSortingGrid_NUM_FOUND_CELLS_GPU; ++i) cellIndicies[i] = indicies;
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
	
	for(int i = 0; i < btFluidSortingGrid_NUM_FOUND_CELLS_GPU; ++i) 
	{
		//Crashes on compiling with Catalyst 13.1 if
		//(btFluidGridIterator){INVALID_FIRST_INDEX, INVALID_FIRST_INDEX} is used directly
		int invalidLowerIndex = INVALID_FIRST_INDEX;
		int invalidUpperIndex = INVALID_LAST_INDEX;
		out_cells[i] = (btFluidGridIterator){invalidLowerIndex, invalidUpperIndex};
		//out_cells[i] = (btFluidGridIterator){INVALID_FIRST_INDEX, INVALID_LAST_INDEX};
	}
	for(int i = 0; i < btFluidSortingGrid_NUM_FOUND_CELLS_GPU; ++i)
	{
	
		btFluidGridPosition lower = cellIndicies[i];
		lower.x--;
	
		btFluidGridPosition upper = cellIndicies[i];
		upper.x++;
		
		int lowerIndex, upperIndex;
		binaryRangeSearch(numActiveCells, cellValues, getCombinedPosition(lower), getCombinedPosition(upper), &lowerIndex, &upperIndex);
		
		if(lowerIndex != numActiveCells)
		{
			out_cells[i] = (btFluidGridIterator){cellContents[lowerIndex].m_firstIndex, cellContents[upperIndex].m_lastIndex};
		}
	
	}
}


__kernel void findNeighborCellsPerCell( __constant int* numActiveCells, __global btFluidGridCombinedPos* cellValues, 
										__global btFluidGridIterator* cellContents, __global btFluidSortingGridFoundCellsGpu* out_foundCells)
{
	int gridCellIndex = get_global_id(0);
	if(gridCellIndex >= *numActiveCells) return;
	
	btFluidGridCombinedPos combinedPosition = cellValues[gridCellIndex];
	
	btFluidGridPosition splitPosition;
	splitCombinedPosition(BT_FLUID_GRID_COORD_RANGE, BT_FLUID_GRID_COORD_RANGE, combinedPosition, 
							&splitPosition.x, &splitPosition.y, &splitPosition.z);
	splitPosition.x -= BT_FLUID_GRID_COORD_RANGE_HALVED;
	splitPosition.y -= BT_FLUID_GRID_COORD_RANGE_HALVED;
	splitPosition.z -= BT_FLUID_GRID_COORD_RANGE_HALVED;
	btFluidGridIterator foundCells[btFluidSortingGrid_NUM_FOUND_CELLS_GPU];
	findCellsFromGridPosition(*numActiveCells, cellValues, cellContents, splitPosition, foundCells);
	for(int cell = 0; cell < btFluidSortingGrid_NUM_FOUND_CELLS_GPU; ++cell) out_foundCells[gridCellIndex].m_iterators[cell] = foundCells[cell];
}

__kernel void findGridCellIndexPerParticle(__constant int* numActiveCells, __global btFluidGridIterator* cellContents, 
											__global int* out_gridCellIndicies)
{
	int gridCellIndex = get_global_id(0);
	if(gridCellIndex >= *numActiveCells) return;
	
	btFluidGridIterator foundCell = cellContents[gridCellIndex];
	for(int n = foundCell.m_firstIndex; n <= foundCell.m_lastIndex; ++n) out_gridCellIndicies[n] = gridCellIndex;
}

//
#define SIMD_EPSILON FLT_EPSILON
__kernel void sphComputePressure(__constant btFluidSphParametersGlobal* FG,  __constant btFluidSphParametersLocal* FL,
								  __global btVector3* fluidPosition, __global btScalar* fluidDensity,
								  __global btFluidSortingGridFoundCellsGpu* foundCells, __global int* foundCellIndex, int numFluidParticles)
{
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	btScalar poly6ZeroDistance = FG->m_sphRadiusSquared * FG->m_sphRadiusSquared * FG->m_sphRadiusSquared;
	btScalar sum = poly6ZeroDistance * FL->m_initialSum;
	
	for(int cell = 0; cell < btFluidSortingGrid_NUM_FOUND_CELLS_GPU; ++cell) 
	{
		btFluidGridIterator foundCell = foundCells[ foundCellIndex[i] ].m_iterators[cell];
		
		for(int n = foundCell.m_firstIndex; n <= foundCell.m_lastIndex; ++n)
		{
			btVector3 delta = (fluidPosition[i] - fluidPosition[n]) * FG->m_simulationScale;	//Simulation scale distance
			btScalar distanceSquared = btVector3_length2(delta);
			
			btScalar c = FG->m_sphRadiusSquared - distanceSquared;
			sum += (c > 0.0f && i != n) ? c*c*c : 0.0f;		//If c is positive, the particle is within interaction radius(poly6 kernel radius)
		}
	}
	
	fluidDensity[i] = sum * FL->m_sphParticleMass * FG->m_poly6KernCoeff;
}


__kernel void sphComputeForce(__constant btFluidSphParametersGlobal* FG, __constant btFluidSphParametersLocal* FL,
							   __global btVector3* fluidPosition, __global btVector3* fluidVelEval, 
							   __global btVector3* fluidSphForce, __global btScalar* fluidDensity,
							   __global btFluidSortingGridFoundCellsGpu* foundCells, __global int* foundCellIndex, int numFluidParticles)
{
	btScalar vterm = FG->m_viscosityKernLapCoeff * FL->m_viscosity;
	
	int i = get_global_id(0);
	if(i >= numFluidParticles) return;
	
	btScalar density_i = fluidDensity[i];
	btScalar invDensity_i = 1.0f / density_i;
	btScalar pressure_i = (density_i - FL->m_restDensity) * FL->m_stiffness;
	
	btVector3 force = {0.0f, 0.0f, 0.0f, 0.0f};
	
	for(int cell = 0; cell < btFluidSortingGrid_NUM_FOUND_CELLS_GPU; ++cell) 
	{
		btFluidGridIterator foundCell = foundCells[ foundCellIndex[i] ].m_iterators[cell];
		
		for(int n = foundCell.m_firstIndex; n <= foundCell.m_lastIndex; ++n)
		{	
			btVector3 delta = (fluidPosition[i] - fluidPosition[n]) * FG->m_simulationScale;	//Simulation scale distance
			btScalar distanceSquared = btVector3_length2(delta);
			
			if(FG->m_sphRadiusSquared > distanceSquared && i != n)
			{
				btScalar density_n = fluidDensity[n];
				btScalar invDensity_n = 1.0f / density_n;
				btScalar pressure_n = (density_n - FL->m_restDensity) * FL->m_stiffness;
			
				btScalar distance = sqrt(distanceSquared);
				btScalar c = FG->m_sphSmoothRadius - distance;
				btScalar pterm = -0.5f * c * FG->m_spikyKernGradCoeff * (pressure_i + pressure_n);
				pterm /= (distance < SIMD_EPSILON) ? SIMD_EPSILON : distance;
				
				btScalar dterm = c * invDensity_i * invDensity_n;
				
				force += (delta * pterm + (fluidVelEval[n] - fluidVelEval[i]) * vterm) * dterm;
			}
		}
	}
	
	fluidSphForce[i] = force * FL->m_sphParticleMass;
}

