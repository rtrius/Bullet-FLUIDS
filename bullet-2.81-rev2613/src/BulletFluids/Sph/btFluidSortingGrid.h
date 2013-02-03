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


#define BT_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT	//Ensure that this is also #defined in "fluidSph.cl"
#ifdef BT_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT
	typedef unsigned long long int btFluidGridUint64;
	typedef btFluidGridUint64 btFluidGridCombinedPos;					//Range must contain BT_FLUID_GRID_COORD_RANGE^3
	const btFluidGridCombinedPos BT_FLUID_GRID_COORD_RANGE = 2097152;	//2^21
#else
	typedef unsigned int btFluidGridCombinedPos;						//Range must contain BT_FLUID_GRID_COORD_RANGE^3
	const btFluidGridCombinedPos BT_FLUID_GRID_COORD_RANGE = 1024;		//2^10
#endif

typedef int btFluidGridCoordinate;
const btFluidGridCoordinate BT_FLUID_GRID_COORD_RANGE_HALVED = BT_FLUID_GRID_COORD_RANGE/2;

///For sorting; contains a btFluidSortingGrid grid cell id and fluid particle index.
struct btFluidGridValueIndexPair
{
	btFluidGridCombinedPos m_value;		///<Grid cell id
	int m_index;						///<Fluid particle index
	
	btFluidGridValueIndexPair() {}
	btFluidGridValueIndexPair(btFluidGridCombinedPos value, int index) : m_value(value), m_index(index) {}
};

///@brief Contains a world scale position quantized to units of btFluidSortingGrid.m_gridCellSize.
struct btFluidGridPosition
{
	btFluidGridCoordinate x;		
	btFluidGridCoordinate y;
	btFluidGridCoordinate z;
private:
	btFluidGridCoordinate padding;
	
public:
	bool operator==(const btFluidGridPosition& GI) const { return (x == GI.x && y == GI.y && z == GI.z); }
	bool operator!=(const btFluidGridPosition& GI) const { return (x != GI.x || y != GI.y || z != GI.z); }
	
	bool operator>(const btFluidGridPosition& GI) const
	{
		if(z != GI.z) return (z > GI.z);
		if(y != GI.y) return (y > GI.y);
		return (x > GI.x);
	}	
	bool operator<(const btFluidGridPosition& GI) const
	{
		if(z != GI.z) return (z < GI.z);
		if(y != GI.y) return (y < GI.y);
		return (x < GI.x);
	}
	
	btFluidGridCombinedPos getCombinedPosition() const
	{
		//Convert range 
		//from [-BT_FLUID_GRID_COORD_RANGE_HALVED, BT_FLUID_GRID_COORD_RANGE_HALVED - 1] 
		//  to [0, BT_FLUID_GRID_COORD_RANGE - 1] before combining
		//e.g. from [-512, 511] to [0, 1023] before combining	
		btFluidGridCoordinate signedX = x + BT_FLUID_GRID_COORD_RANGE_HALVED;
		btFluidGridCoordinate signedY = y + BT_FLUID_GRID_COORD_RANGE_HALVED;
		btFluidGridCoordinate signedZ = z + BT_FLUID_GRID_COORD_RANGE_HALVED;
		
		btFluidGridCombinedPos unsignedX = static_cast<btFluidGridCombinedPos>(signedX);
		btFluidGridCombinedPos unsignedY = static_cast<btFluidGridCombinedPos>(signedY) * BT_FLUID_GRID_COORD_RANGE;
		btFluidGridCombinedPos unsignedZ = static_cast<btFluidGridCombinedPos>(signedZ) * BT_FLUID_GRID_COORD_RANGE * BT_FLUID_GRID_COORD_RANGE;
		
		return unsignedX + unsignedY + unsignedZ;
	}
};



///@brief Uniform grid broadphase for btFluidSph particles.
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
///of particles are quantized into grid cell coordinates(btFluidGridPosition). Those
///integer coordinates are then converted into single values that are used for sorting
///(btFluidGridCombinedPos). After sorting the particles by the grid cell they are contained in, 
///the lower and upper indicies of each nonempty cell is detected and stored(btFluidGridIterator).
///@par
///Effective size: BT_FLUID_GRID_COORD_RANGE^3, which is currently 1024^3 
///or 2^21^3(with #define BT_ENABLE_FLUID_SORTING_GRID_LARGE_WORLD_SUPPORT) grid cells.
///Worlds larger than 2^21^3 are unsupported.
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

	btAlignedObjectArray<btFluidGridCombinedPos> m_activeCells;		//Stores the value of each nonempty grid cell
	btAlignedObjectArray<btFluidGridIterator> m_cellContents;	//Stores the range of indicies that correspond to the values in m_activeCells
	
	btAlignedObjectArray<btFluidGridValueIndexPair> m_tempPairs;
	btAlignedObjectArray<btVector3> m_tempBufferVector;
	btAlignedObjectArray<void*> m_tempBufferVoid;
	
public:
	btFluidSortingGrid() : m_pointMin(0,0,0), m_pointMax(0,0,0), m_gridCellSize(1) {}

	void insertParticles(btFluidParticles& fluids);
	
	void clear() 
	{ 
		m_activeCells.resize(0);
		m_cellContents.resize(0);
		for(int i = 0; i < btFluidSortingGrid::NUM_MULTITHREADING_GROUPS; ++i) m_multithreadingGroups[i].resize(0);
	}
	
	///Returns a 3x3x3 group of btFluidGridIterator, which is the maximum extent of cells
	///that may interact with an AABB defined by (position - radius, position + radius). 
	///Where radius is the SPH smoothing radius, in btFluidSphParametersGlobal, converted to world scale.
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
	
	struct AabbCallback
	{
		AabbCallback() {}
		virtual ~AabbCallback() {}
		
		///btFluidSortingGrid::forEachGridCell() will continue calling processParticles() if this returns true
		virtual bool processParticles(const btFluidGridIterator FI, const btVector3& aabbMin, const btVector3& aabbMax) = 0;
	};
	void forEachGridCell(const btVector3& aabbMin, const btVector3& aabbMax, btFluidSortingGrid::AabbCallback& callback) const;

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
	
	btAlignedObjectArray<btFluidGridCombinedPos>& internalGetActiveCells() { return m_activeCells; }
	btAlignedObjectArray<btFluidGridIterator>& internalGetCellContents() { return m_cellContents; }
	
	const btAlignedObjectArray<int>& internalGetMultithreadingGroup(int index) const { return m_multithreadingGroups[index]; }
	
	btVector3& internalGetPointAabbMin() { return m_pointMin; }
	btVector3& internalGetPointAabbMax() { return m_pointMax; }
	
private:
	btFluidGridPosition getDiscretePosition(const btVector3& position) const;

	void findAdjacentGridCells(btFluidGridPosition indicies, btFluidSortingGrid::FoundCells& out_gridCells) const;
	void findAdjacentGridCellsSymmetric(btFluidGridPosition indicies, btFluidSortingGrid::FoundCells& out_gridCells) const;
	
	void generateMultithreadingGroups();
};

#endif



