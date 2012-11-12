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
#ifndef BT_FLUID_SORTING_GRID_H
#define BT_FLUID_SORTING_GRID_H

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btQuickProf.h"

#include "btFluidParticles.h"

class btVector3;
struct btFluidParticles;


///@brief Used to iterate through all particles in a single btFluidSortingGrid cell.
///@remarks
///The standard method for iterating through a grid cell is:
///@code
///		btFluidGridIterator FI;
///		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
///@endcode
struct btFluidGridIterator
{
	int m_firstIndex;
	int m_lastIndex;
	
	btFluidGridIterator() {}
	btFluidGridIterator(int firstIndex, int lastIndex) : m_firstIndex(firstIndex), m_lastIndex(lastIndex) {}
};


//#define SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED	//Ensure that this is also #defined in "fluids.cl"
#ifdef SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED
	typedef unsigned long long int btSortGridUint64;
	typedef btSortGridUint64 btSortGridValue;			//Range must contain SORT_GRID_INDEX_RANGE^3
	typedef int btSortGridIndex;
	const btSortGridValue SORT_GRID_INDEX_RANGE = 2097152;	//2^21
#else
	typedef unsigned int btSortGridValue;				//Range must contain SORT_GRID_INDEX_RANGE^3
	typedef char btSortGridIndex;
	const btSortGridValue SORT_GRID_INDEX_RANGE = 256;			//2^( 8*sizeof(btSortGridIndex) )
#endif

//btSortGridIndex_large is a signed type with range including all values in:
//	[-HALVED_SORT_GRID_INDEX_RANGE, (HALVED_SORT_GRID_INDEX_RANGE - 1) + HALVED_SORT_GRID_INDEX_RANGE]
//
//	e.g if SORT_GRID_INDEX_RANGE == 256, HALVED_SORT_GRID_INDEX_RANGE == 128
//	then btSortGridIndex_large must contain [-128, 127 + 128] == [-128, 255].
//	It is used to convert from btSortGridIndex(signed, small range), to btSortGridValue(unsigned, large range).
typedef int btSortGridIndex_large;	

const btSortGridIndex_large HALVED_SORT_GRID_INDEX_RANGE = SORT_GRID_INDEX_RANGE/2;
//const int BITS_PER_SORT_GRID_INDEX = 8 * sizeof(btSortGridIndex);

struct btValueIndexPair
{
	btSortGridValue m_value;
	int m_index;
	
	btValueIndexPair() {}
	btValueIndexPair(btSortGridValue value, int index) : m_value(value), m_index(index) {}
};

///@brief Contains a world scale position quantized to units of btFluidSortingGrid.m_gridCellSize.
struct btSortGridIndicies
{
	btSortGridIndex x;		
	btSortGridIndex y;
	btSortGridIndex z;
private:
	btSortGridIndex padding;
	
public:
	bool operator==(const btSortGridIndicies& GI) const { return (x == GI.x && y == GI.y && z == GI.z); }
	bool operator!=(const btSortGridIndicies& GI) const { return (x != GI.x || y != GI.y || z != GI.z); }
	
	bool operator>(const btSortGridIndicies& GI) const
	{
		if(z != GI.z) return (z > GI.z);
		if(y != GI.y) return (y > GI.y);
		return (x > GI.x);
	}	
	bool operator<(const btSortGridIndicies& GI) const
	{
		if(z != GI.z) return (z < GI.z);
		if(y != GI.y) return (y < GI.y);
		return (x < GI.x);
	}
	
	btSortGridValue getValue() const
	{
		//Convert range 
		//from [-HALVED_SORT_GRID_INDEX_RANGE, HALVED_SORT_GRID_INDEX_RANGE - 1] 
		//  to [0, SORT_GRID_INDEX_RANGE - 1] before combining
		//e.g. from [-128, 127] to [0, 255] before combining	
		btSortGridIndex_large signedX = static_cast<btSortGridIndex_large>(x) + HALVED_SORT_GRID_INDEX_RANGE;
		btSortGridIndex_large signedY = static_cast<btSortGridIndex_large>(y) + HALVED_SORT_GRID_INDEX_RANGE;
		btSortGridIndex_large signedZ = static_cast<btSortGridIndex_large>(z) + HALVED_SORT_GRID_INDEX_RANGE;
		
		btSortGridValue unsignedX = static_cast<btSortGridValue>(signedX);
		btSortGridValue unsignedY = static_cast<btSortGridValue>(signedY) * SORT_GRID_INDEX_RANGE;
		btSortGridValue unsignedZ = static_cast<btSortGridValue>(signedZ) * SORT_GRID_INDEX_RANGE * SORT_GRID_INDEX_RANGE;
		
		return unsignedX + unsignedY + unsignedZ;
	}
};



///@brief Uniform grid broadphase for fluid particles.
///@remarks
///A fundamental operation in SPH fluid simulations is the detection
///of collisions between fluid particles.
///@par
///Since testing each fluid pair would require O(n^2) operations, a grid based 
///broadphase is implemented to accelerate the search to O(kn), where k
///is the average number of particles in each cell. Each grid cell has a size 
///of r, where r is the SPH interaction radius at world scale. When a particle 
///is queried, it searches a 3x3x3 grid cell volume surrounding its position.
///Note that, for each particle, a 'collision' is detected if the center of
///other particles is within r; that is, the effective radius of collision,
///were all particles to be treated as spheres(and not as points), is r/2.
///@par
///Particles are stored as a set of index ranges. First, the world scale positions
///of particles are quantized into grid cell coordinates(btSortGridIndicies). Those
///integer coordinates are then converted into single values that are used for sorting
///(btSortGridValue). After sorting the particles by the grid cell they are contained in, 
///the lower and upper indicies of each nonempty cell is detected and stored(btFluidGridIterator).
///@par
///Effective size: SORT_GRID_INDEX_RANGE^3, which is currently 255^3 
///or 2^21^3(with #define SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED) grid cells.
class btFluidSortingGrid
{
	//INVALID_LAST_INDEX must be lower than INVALID_FIRST_INDEX,
	//such that the below loop will not execute.
	//		for(int i = FI.m_firstIndex; i <= FI.m_lastIndex; ++i)
	static const int INVALID_FIRST_INDEX = -1;
	static const int INVALID_LAST_INDEX = INVALID_FIRST_INDEX - 1;
public:
	static const int NUM_MULTITHREADING_GROUPS = 27; 	///<Number of grid cells that may be accessed when iterating through a single grid cell
	static const int NUM_FOUND_CELLS = 27;				///<Number of grid cells returned from btFluidSortingGrid::findCells()
	static const int NUM_FOUND_CELLS_SYMMETRIC = 14;	///<Number of grid cells returned from btFluidSortingGrid::findCellsSymmetric()

	struct FoundCells { btFluidGridIterator m_iterators[btFluidSortingGrid::NUM_FOUND_CELLS]; }; ///<Contains results of btFluidSortingGrid::findCells()
	
private:
	btVector3 m_pointMin;	//AABB calculated from the center of fluid particles, without considering particle radius
	btVector3 m_pointMax;
	
	///Each array contains a set of grid cell indicies that may be simultaneously processed
	btAlignedObjectArray<int> m_multithreadingGroups[btFluidSortingGrid::NUM_MULTITHREADING_GROUPS];
	
	btScalar m_gridCellSize;

	btAlignedObjectArray<btSortGridValue> m_activeCells;		//Stores the value of each nonempty grid cell
	btAlignedObjectArray<btFluidGridIterator> m_cellContents;	//Stores the range of indicies that correspond to the values in m_activeCells
	
	btAlignedObjectArray<btValueIndexPair> m_tempPairs;
	btAlignedObjectArray<btVector3> m_tempBufferVector;
	btAlignedObjectArray<void*> m_tempBufferVoid;
public:
	btFluidSortingGrid() {}
	btFluidSortingGrid(btScalar simulationScale, btScalar sphSmoothRadius) { setCellSize(simulationScale, sphSmoothRadius); }

	void insertParticles(btFluidParticles& fluids);
	
	void clear() 
	{ 
		m_activeCells.resize(0);
		m_cellContents.resize(0);
		for(int i = 0; i < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++i) m_multithreadingGroups[i].resize(0);
	}
	
	///Returns a 3x3x3 group of btFluidGridIterator, which is the maximum extent of cells
	///that may interact with an AABB defined by (position - radius, position + radius). 
	///Where radius is the SPH smoothing radius, in btFluidParametersGlobal, converted to world scale.
	///@param position Center of the AABB defined by (position - radius, position + radius).
	void findCells(const btVector3& position, btFluidSortingGrid::FoundCells& out_gridCells) const
	{
		findAdjacentGridCells( getDiscretePosition(position), out_gridCells );
	}
	
	///Returns 14 grid cells, with out_gridCells->m_iterator[0] as the center cell corresponding to position.
	void findCellsSymmetric(const btVector3& position, btFluidSortingGrid::FoundCells& out_gridCells) const
	{
		findAdjacentGridCellsSymmetric( getDiscretePosition(position), out_gridCells );
	}
	
	int getNumGridCells() const { return m_activeCells.size(); }	///<Returns the number of nonempty grid cells.
	btFluidGridIterator getGridCell(int gridCellIndex) const { return m_cellContents[gridCellIndex]; }
	void getGridCellIndiciesInAabb(const btVector3& min, const btVector3& max, btAlignedObjectArray<int>& out_indicies) const;
	
	btScalar getCellSize() const { return m_gridCellSize; }
	void setCellSize(btScalar simulationScale, btScalar sphSmoothRadius) 
	{
		m_gridCellSize = sphSmoothRadius / simulationScale; 	//Divide by simulationScale to convert to world scale
	}
	
	///Returns the AABB calculated from the center of each fluid particle in the grid, without considering particle radius
	void getPointAabb(btVector3& out_pointMin, btVector3& out_pointMax) const
	{
		out_pointMin = m_pointMin;
		out_pointMax = m_pointMax;
	}
	
	btAlignedObjectArray<btSortGridValue>& internalGetActiveCells() { return m_activeCells; }
	btAlignedObjectArray<btFluidGridIterator>& internalGetCellContents() { return m_cellContents; }
	
	const btAlignedObjectArray<int>& internalGetMultithreadingGroup(int index) const { return m_multithreadingGroups[index]; }
	
private:
	btSortGridIndicies getDiscretePosition(const btVector3& position) const;

	void findAdjacentGridCells(btSortGridIndicies indicies, btFluidSortingGrid::FoundCells& out_gridCells) const;
	void findAdjacentGridCellsSymmetric(btSortGridIndicies indicies, btFluidSortingGrid::FoundCells& out_gridCells) const;
	
	void generateMultithreadingGroups();
};

#endif



