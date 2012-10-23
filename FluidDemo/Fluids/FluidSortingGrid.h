/* FluidSortingGrid.h

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
#ifndef FLUID_SORTING_GRID_H
#define FLUID_SORTING_GRID_H

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"
#include "LinearMath/btQuickProf.h"

#include "FluidParticles.h"

class btVector3;
struct FluidParticles;


///@brief Used to iterate through all particles in a single FluidSortingGrid cell.
///@remarks
///The standard method for iterating through a grid cell is:
///@code
///		FluidGridIterator FI;
///		for(int n = FI.m_firstIndex; n <= FI.m_lastIndex; ++n)
///@endcode
struct FluidGridIterator
{
	int m_firstIndex;
	int m_lastIndex;
	
	FluidGridIterator() {}
	FluidGridIterator(int firstIndex, int lastIndex) : m_firstIndex(firstIndex), m_lastIndex(lastIndex) {}
};


//#define SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED	//Ensure that this is also #defined in "fluids.cl"
#ifdef SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED
	typedef unsigned long long int SortGridUint64;
	typedef SortGridUint64 SortGridValue;			//Range must contain SORT_GRID_INDEX_RANGE^3
	typedef int SortGridIndex;
	const SortGridValue SORT_GRID_INDEX_RANGE = 2097152;	//2^21
#else
	typedef unsigned int SortGridValue;				//Range must contain SORT_GRID_INDEX_RANGE^3
	typedef char SortGridIndex;
	const SortGridValue SORT_GRID_INDEX_RANGE = 256;			//2^( 8*sizeof(SortGridIndex) )
#endif

//SortGridIndex_large is a signed type with range including all values in:
//	[-HALVED_SORT_GRID_INDEX_RANGE, (HALVED_SORT_GRID_INDEX_RANGE - 1) + HALVED_SORT_GRID_INDEX_RANGE]
//
//	e.g if SORT_GRID_INDEX_RANGE == 256, HALVED_SORT_GRID_INDEX_RANGE == 128
//	then SortGridIndex_large must contain [-128, 127 + 128] == [-128, 255].
//	It is used to convert from SortGridIndex(signed, small range), to SortGridValue(unsigned, large range).
typedef int SortGridIndex_large;	

const SortGridIndex_large HALVED_SORT_GRID_INDEX_RANGE = SORT_GRID_INDEX_RANGE/2;
//const int BITS_PER_SORT_GRID_INDEX = 8 * sizeof(SortGridIndex);

struct ValueIndexPair
{
	SortGridValue m_value;
	int m_index;
	
	ValueIndexPair() {}
	ValueIndexPair(SortGridValue value, int index) : m_value(value), m_index(index) {}
};

///@brief Contains a world scale position quantized to units of FluidSortingGrid.m_gridCellSize.
struct SortGridIndicies
{
	SortGridIndex x;		
	SortGridIndex y;
	SortGridIndex z;
private:
	SortGridIndex padding;
	
public:
	bool operator==(const SortGridIndicies& GI) const { return (x == GI.x && y == GI.y && z == GI.z); }
	bool operator!=(const SortGridIndicies& GI) const { return (x != GI.x || y != GI.y || z != GI.z); }
	
	bool operator>(const SortGridIndicies& GI) const
	{
		if(z != GI.z) return (z > GI.z);
		if(y != GI.y) return (y > GI.y);
		return (x > GI.x);
	}	
	bool operator<(const SortGridIndicies& GI) const
	{
		if(z != GI.z) return (z < GI.z);
		if(y != GI.y) return (y < GI.y);
		return (x < GI.x);
	}
	
	SortGridValue getValue() const
	{
		//Convert range 
		//from [-HALVED_SORT_GRID_INDEX_RANGE, HALVED_SORT_GRID_INDEX_RANGE - 1] 
		//  to [0, SORT_GRID_INDEX_RANGE - 1] before combining
		//e.g. from [-128, 127] to [0, 255] before combining	
		SortGridIndex_large signedX = static_cast<SortGridIndex_large>(x) + HALVED_SORT_GRID_INDEX_RANGE;
		SortGridIndex_large signedY = static_cast<SortGridIndex_large>(y) + HALVED_SORT_GRID_INDEX_RANGE;
		SortGridIndex_large signedZ = static_cast<SortGridIndex_large>(z) + HALVED_SORT_GRID_INDEX_RANGE;
		
		SortGridValue unsignedX = static_cast<SortGridValue>(signedX);
		SortGridValue unsignedY = static_cast<SortGridValue>(signedY) * SORT_GRID_INDEX_RANGE;
		SortGridValue unsignedZ = static_cast<SortGridValue>(signedZ) * SORT_GRID_INDEX_RANGE * SORT_GRID_INDEX_RANGE;
		
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
///of particles are quantized into grid cell coordinates(SortGridIndicies). Those
///integer coordinates are then converted into single values that are used for sorting
///(SortGridValue). After sorting the particles by the grid cell they are contained in, 
///the lower and upper indicies of each nonempty cell is detected and stored(FluidGridIterator).
///@par
///Effective size: SORT_GRID_INDEX_RANGE^3, which is currently 255^3 
///or 2^21^3(with #define SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED) grid cells.
class FluidSortingGrid
{
	//INVALID_LAST_INDEX must be lower than INVALID_FIRST_INDEX,
	//such that the below loop will not execute.
	//		for(int i = FI.m_firstIndex; i <= FI.m_lastIndex; ++i)
	static const int INVALID_FIRST_INDEX = -1;
	static const int INVALID_LAST_INDEX = INVALID_FIRST_INDEX - 1;
public:
	static const int NUM_CELL_PROCESSING_GROUPS = 27; 	///<Number of grid cells that may be accessed when iterating through a single grid cell
	static const int NUM_FOUND_CELLS = 27;				///<Number of grid cells returned from FluidSortingGrid::findCells()
	static const int NUM_FOUND_CELLS_SYMMETRIC = 14;	///<Number of grid cells returned from FluidSortingGrid::findCellsSymmetric()

	struct FoundCells { FluidGridIterator m_iterators[FluidSortingGrid::NUM_FOUND_CELLS]; }; ///<Contains results of FluidSortingGrid::findCells()
	
private:
	btVector3 m_pointMin;	//AABB calculated from the center of fluid particles, without considering particle radius
	btVector3 m_pointMax;
	
	btAlignedObjectArray<int> m_cellProcessingGroups[FluidSortingGrid::NUM_CELL_PROCESSING_GROUPS];
	
	btScalar m_gridCellSize;

	btAlignedObjectArray<SortGridValue> m_activeCells;		//Stores the value of each nonempty grid cell
	btAlignedObjectArray<FluidGridIterator> m_cellContents;	//Stores the range of indicies that correspond to the values in m_activeCells
	
public:
	FluidSortingGrid() {}
	FluidSortingGrid(btScalar simulationScale, btScalar sphSmoothRadius) { setup(simulationScale, sphSmoothRadius); }

	void setup(btScalar simulationScale, btScalar sphSmoothRadius) 
	{	
		btScalar worldCellSize = sphSmoothRadius / simulationScale;
		
		m_gridCellSize = worldCellSize;
	}

	void clear() 
	{ 
		m_activeCells.resize(0);
		m_cellContents.resize(0);
	}
	void insertParticles(FluidParticles *fluids);
	
	///Returns a 3x3x3 group of FluidGridIterator, which is the maximum extent of cells
	///that may interact with an AABB defined by (position - radius, position + radius). 
	///Where radius is the SPH smoothing radius, in FluidParametersGlobal, converted to world scale.
	///@param position Center of the AABB defined by (position - radius, position + radius).
	void findCells(const btVector3 &position, FluidSortingGrid::FoundCells *out_gridCells) const
	{
		findAdjacentGridCells( generateIndicies(position), out_gridCells );
	}
	
	///Returns 14 grid cells, with out_gridCells->m_iterator[0] as the center cell corresponding to position.
	void findCellsSymmetric(const btVector3 &position, FluidSortingGrid::FoundCells *out_gridCells) const
	{
		findAdjacentGridCellsSymmetric( generateIndicies(position), out_gridCells );
	}
	
	int getNumGridCells() const { return m_activeCells.size(); }	///<Returns the number of nonempty grid cells.
	FluidGridIterator getGridCell(int gridCellIndex) const
	{
		return FluidGridIterator(m_cellContents[gridCellIndex].m_firstIndex, m_cellContents[gridCellIndex].m_lastIndex);
	}
	void getGridCellIndiciesInAabb(const btVector3 &min, const btVector3 &max, btAlignedObjectArray<int> *out_indicies) const;
	
	btScalar getCellSize() const { return m_gridCellSize; }
	
	///Returns the AABB calculated from the center of each fluid particle in the grid, without considering particle radius
	void getPointAabb(btVector3 *out_pointMin, btVector3 *out_pointMax) const
	{
		*out_pointMin = m_pointMin;
		*out_pointMax = m_pointMax;
	}
	
	void internalGetIndiciesReduce(int gridCellIndex, int *out_x, int *out_y, int *out_z) const
	{
		splitIndex(SORT_GRID_INDEX_RANGE, SORT_GRID_INDEX_RANGE, m_activeCells[gridCellIndex], out_x, out_y, out_z);
	}
	
	btAlignedObjectArray<SortGridValue>& internalGetActiveCells() { return m_activeCells; }
	btAlignedObjectArray<FluidGridIterator>& internalGetCellContents() { return m_cellContents; }
	
	btAlignedObjectArray<int>& internalGetCellProcessingGroup(int index) { return m_cellProcessingGroups[index]; }
	const btAlignedObjectArray<int>& internalGetCellProcessingGroup(int index) const { return m_cellProcessingGroups[index]; }
	
private:
	const FluidGridIterator* getCell(SortGridValue value) const
	{
		int index = m_activeCells.findBinarySearch(value);	//findBinarySearch() returns m_activeCells.size() on failure
		return ( index != m_activeCells.size() ) ? &m_cellContents[index] : 0;
	}

	SortGridIndicies generateIndicies(const btVector3 &position) const;

	void findAdjacentGridCells(SortGridIndicies indicies, FluidSortingGrid::FoundCells *out_gridCells) const;
	void findAdjacentGridCellsSymmetric(SortGridIndicies indicies, FluidSortingGrid::FoundCells *out_gridCells) const;
	
	void generateCellProcessingGroups();
	
	
	void resetPointAabb()
	{
		m_pointMin.setValue(BT_LARGE_FLOAT, BT_LARGE_FLOAT, BT_LARGE_FLOAT);
		m_pointMax.setValue(-BT_LARGE_FLOAT, -BT_LARGE_FLOAT, -BT_LARGE_FLOAT);
	}
	
	void updatePointAabb(const btVector3 &position)
	{
		for(int i = 0; i < 3; ++i)
		{
			if(position.m_floats[i] < m_pointMin.m_floats[i]) m_pointMin.m_floats[i] = position.m_floats[i];
			if(position.m_floats[i] > m_pointMax.m_floats[i]) m_pointMax.m_floats[i] = position.m_floats[i];
		}
	}
	
	
#ifndef SORTING_GRID_LARGE_WORLD_SUPPORT_ENABLED
	//Extracts the (x,y,z) indicies from combinedIndex, where
	//combinedIndex == x + y*resolutionX + z*resolutionX*resolutionY
	static void splitIndex(int resolutionX, int resolutionY, int combinedIndex, int *out_x, int *out_y, int *out_z)
	{
		int x = combinedIndex % resolutionX;
		int z = combinedIndex / (resolutionX*resolutionY);
		int y = (combinedIndex - z*resolutionX*resolutionY) / resolutionX;
				
		*out_x = x;
		*out_z = z;
		*out_y = y;
	}
#else
	static void splitIndex(unsigned long long int resolutionX, unsigned long long int resolutionY, 
						   unsigned long long int combinedIndex, int *out_x, int *out_y, int *out_z)
	{
		unsigned long long int cellsPerLine = resolutionX;
		unsigned long long int cellsPerPlane = resolutionX * resolutionY;
											 
		unsigned long long int x = combinedIndex % cellsPerLine;
		unsigned long long int z = combinedIndex / cellsPerPlane;
		unsigned long long int y = (combinedIndex - z*cellsPerPlane) / cellsPerLine;
				
		*out_x = static_cast<int>(x);
		*out_z = static_cast<int>(z);
		*out_y = static_cast<int>(y);
	}
#endif
};

#endif



