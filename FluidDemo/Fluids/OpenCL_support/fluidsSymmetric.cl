/* fluidsSymmetric.cl

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


#define FluidSortingGrid_NUM_FOUND_CELLS 14

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

inline void binaryRangeSearch(int numActiveCells, __global SortGridValue *cellValues,
							  SortGridValue lowerValue, SortGridValue upperValue, int *out_lowerIndex, int *out_upperIndex)
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

inline void findCells(int numActiveCells, __global SortGridValue *cellValues, __global FluidGridIterator *cellContents, 
					  btScalar cellSize, btVector3 position, FluidGridIterator *out_cells)
{
	//Due to symmetry, only 14 of 27 adjacent cells are found
	for(int i = 0; i < FluidSortingGrid_NUM_FOUND_CELLS; ++i) out_cells[i] = (FluidGridIterator){-1, -2};

	SortGridIndicies indicies = getSortGridIndicies(cellSize, position);
	SortGridIndicies centers[6];	//Center of the 6 rows
	for(int i = 0; i < 6; ++i) centers[i] = indicies;
	
	//centers[0] and centers[1] contain 2-cell rows
	centers[1].y++;
	
	//centers[2] to centers[4] contain 3-cell rows
	centers[2].z--;
	
	centers[3].z--;
	centers[3].y++;
	
	centers[4].z--;
	centers[4].y--;
	
	//centers[5] contains 1 cell
	centers[5].y--;
	centers[5].x--;
	
	//centers[0]
	{
		SortGridIndicies lower = centers[0];
		lower.x--;
	
		SortGridIndicies upper = centers[0];
		SortGridValue centerValue = getSortGridValue(centers[0]);
		
		int lowerIndex, upperIndex;
		binaryRangeSearch(numActiveCells, cellValues, getSortGridValue(lower), getSortGridValue(upper), &lowerIndex, &upperIndex);
		if(lowerIndex != numActiveCells)
		{
			//out_cells[0] must be the center grid cell if it exists, and INVALID_ITERATOR((FluidGridIterator){-1, -2}) otherwise
			if(lowerIndex != upperIndex) 
			{
				out_cells[0] = cellContents[upperIndex];
				out_cells[1] = cellContents[lowerIndex];
			}
			else
			{
				if(centerValue == cellValues[lowerIndex])out_cells[0] = cellContents[lowerIndex];
				else out_cells[1] = cellContents[lowerIndex];
			}
		}
	}
	
	//centers[1]
	{
		SortGridIndicies lower = centers[1];
		lower.x--;
	
		SortGridIndicies upper = centers[1];
		
		int lowerIndex, upperIndex;
		binaryRangeSearch(numActiveCells, cellValues, getSortGridValue(lower), getSortGridValue(upper), &lowerIndex, &upperIndex);
		if(lowerIndex != numActiveCells)
		{
			out_cells[2] = cellContents[lowerIndex];
			if(lowerIndex != upperIndex) out_cells[3] = cellContents[upperIndex];
		}
	}
	
	//centers[2-4]
	for(int i = 0; i < 3; ++i)
	{
		SortGridIndicies lower = centers[i+2];
		lower.x--;
	
		SortGridIndicies upper = centers[i+2];
		upper.x++;
		
		int lowerIndex, upperIndex;
		binaryRangeSearch(numActiveCells, cellValues, getSortGridValue(lower), getSortGridValue(upper), &lowerIndex, &upperIndex);
		
		if(lowerIndex != numActiveCells)
		{
			int range = upperIndex - lowerIndex;
			
			for(int n = 0; n <= range; ++n) out_cells[4 + i*3 + n] = cellContents[lowerIndex + n];
		}
	}
	
	//centers[5]
	{
		SortGridIndicies lowerAndUpper = centers[5];
		
		int lowerIndex, upperIndex;
		binaryRangeSearch(numActiveCells, cellValues, getSortGridValue(lowerAndUpper), getSortGridValue(lowerAndUpper), &lowerIndex, &upperIndex);
		if(lowerIndex != numActiveCells) out_cells[13] = cellContents[lowerIndex];
	}
	
	for(int i = 14; i < FluidSortingGrid_NUM_FOUND_CELLS; ++i) out_cells[i] = (FluidGridIterator){-1, -2};
}

//
__kernel void clearInvDensityAndNeighbors(__global btScalar *fluidInvDensity, __global FluidNeighbors *fluidNeighbors)
{
	int i = get_global_id(0);
	
	fluidInvDensity[i] = 0.0f;
	fluidNeighbors[i].m_count = 0;
}
__kernel void sphComputePressure(__global FluidParametersGlobal *FG,  __global FluidParametersLocal *FL,
								  __global btVector3 *fluidPosition, __global btScalar *fluidPressure, 
								  __global btScalar *fluidInvDensity,  __global FluidNeighbors *fluidNeighbors,
								  __global int *numActiveCells, __global SortGridValue *cellValues, 
								  __global FluidGridIterator *cellContents, btScalar cellSize, __global int *cellProcessingGroup)
{
	int gridCellIndex = cellProcessingGroup[ get_global_id(0) ];
	
	FluidGridIterator currentCell = cellContents[gridCellIndex];
	if(currentCell.m_firstIndex <= currentCell.m_lastIndex)	//if cell is not empty
	{
		FluidGridIterator foundCells[FluidSortingGrid_NUM_FOUND_CELLS];	//	may be allocated in global memory(slow)
		findCells(*numActiveCells, cellValues, cellContents, cellSize, fluidPosition[currentCell.m_firstIndex], foundCells);
		
		for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
		{
			//Remove particle, with index i, from grid cell
			++foundCells[0].m_firstIndex; //local cell
			
			for(int cell = 0; cell < FluidSortingGrid_NUM_FOUND_CELLS; ++cell) 
			{
				FluidGridIterator foundCell = foundCells[cell];
				
				for(int n = foundCell.m_firstIndex; n <= foundCell.m_lastIndex; ++n)
				{
					btVector3 distance = (fluidPosition[i] - fluidPosition[n]) * FG->m_simulationScale;	//Simulation scale distance
					btScalar distanceSquared = btVector3_length2(distance);
					
					if(FG->m_sphRadiusSquared > distanceSquared) 
					{
						btScalar c = FG->m_sphRadiusSquared - distanceSquared;
						btScalar c_cubed = c * c * c;
						fluidInvDensity[i] += c_cubed;
						fluidInvDensity[n] += c_cubed;
						
						if(fluidNeighbors[i].m_count < MAX_NEIGHBORS)	//if( !m_neighborTable[i].isFilled() ) 
						{	
							//m_neighborTable[i].addNeighbor( n, sqrt(distanceSquared) );
							fluidNeighbors[i].m_particleIndicies[ fluidNeighbors[i].m_count ] = n;
							fluidNeighbors[i].m_distances[ fluidNeighbors[i].m_count ] = sqrt(distanceSquared);
							++fluidNeighbors[i].m_count;
						}
						else if(fluidNeighbors[n].m_count < MAX_NEIGHBORS)
						{
							fluidNeighbors[n].m_particleIndicies[ fluidNeighbors[n].m_count ] = i;
							fluidNeighbors[n].m_distances[ fluidNeighbors[n].m_count ] = sqrt(distanceSquared);
							++fluidNeighbors[n].m_count;
						}
						else break;
					}
				}
			}
		}
	}
}
__kernel void computePressureAndInvDensity(__global FluidParametersGlobal *FG, __global FluidParametersLocal *FL,
										   __global btScalar *fluidPressure, __global btScalar *fluidInvDensity)
{
	int i = get_global_id(0);
	btScalar density = fluidInvDensity[i] * FL->m_particleMass * FG->m_poly6KernCoeff;
	
	fluidPressure[i] = (density - FL->m_restDensity) * FL->m_stiffness;
	fluidInvDensity[i] = 1.0f / density;
}

__kernel void clearSphForce( __global btVector3 *fluidSphForce)
{
	int i = get_global_id(0);
	
	fluidSphForce[i] = (btVector3){0.0f, 0.0f, 0.0f, 0.0f};
}
__kernel void sphComputeForce(__global FluidParametersGlobal *FG, __global FluidParametersLocal *FL,
							   __global btVector3 *fluidPosition, __global btVector3 *fluidVelEval, 
							   __global btVector3 *fluidSphForce, __global btScalar *fluidPressure, 
							   __global btScalar *fluidInvDensity, __global FluidNeighbors *fluidNeighbors,
							   __global FluidGridIterator *cellContents, __global int *cellProcessingGroup)
{
	btScalar vterm = FG->m_viscosityKernLapCoeff * FL->m_viscosity;
	
	int gridCellIndex = cellProcessingGroup[ get_global_id(0) ];
	
	FluidGridIterator currentCell = cellContents[gridCellIndex];
	for(int i = currentCell.m_firstIndex; i <= currentCell.m_lastIndex; ++i)
	{
		for(int j = 0; j < fluidNeighbors[i].m_count; j++ ) 
		{
			int n = fluidNeighbors[i].m_particleIndicies[j];
			
			btVector3 distance = (fluidPosition[i] - fluidPosition[n]) * FG->m_simulationScale;	//Simulation scale distance
			
			btScalar c = FG->m_sphSmoothRadius - fluidNeighbors[i].m_distances[j];
			btScalar pterm = -0.5f * c * FG->m_spikyKernGradCoeff 
							* ( fluidPressure[i] + fluidPressure[n] ) / fluidNeighbors[i].m_distances[j];
			btScalar dterm = c * fluidInvDensity[i] * fluidInvDensity[n];
			
			btVector3 sphForce;
			//sphForce = (distance * pterm + (fluidVelEval[n] - fluidVelEval[i]) * vterm) * dterm;
			sphForce.x = ( pterm * distance.x + vterm * (fluidVelEval[n].x - fluidVelEval[i].x) ) * dterm;
			sphForce.y = ( pterm * distance.y + vterm * (fluidVelEval[n].y - fluidVelEval[i].y) ) * dterm;
			sphForce.z = ( pterm * distance.z + vterm * (fluidVelEval[n].z - fluidVelEval[i].z) ) * dterm;
			
			fluidSphForce[i] += sphForce;
			fluidSphForce[n] += -sphForce;
		}
	}
}


