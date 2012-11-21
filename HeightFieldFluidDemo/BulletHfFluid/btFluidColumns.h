#ifndef BT_FLUID_COLUMNS_H
#define BT_FLUID_COLUMNS_H

#include "LinearMath/btScalar.h"
#include "LinearMath/btAlignedObjectArray.h"

///Stores configuration of a single btHfFluid.
struct btFluidHfParameters
{
	btScalar m_globalVelocityU;
	btScalar m_globalVelocityV;
	btScalar m_gravity;
	btScalar m_volumeDisplacementScale;
	btScalar m_horizontalVelocityScale;

	btScalar m_epsHeight;
	btScalar m_epsEta;		//If m_eta is lower than this, the cell is considered empty/inactive
	
	btFluidHfParameters()
	{
		m_globalVelocityU = btScalar(0.0);
		m_globalVelocityV = btScalar(0.0);
		m_gravity = btScalar(-9.8);

		m_volumeDisplacementScale = btScalar(0.5);
		m_horizontalVelocityScale = btScalar(0.5);
		
		m_epsHeight = btScalar(0.001);
		m_epsEta = btScalar(0.01);
	}
};

///Coordinates the parallel arrays used to represent a heightfield fluid.
struct btFluidColumns
{
	int m_numNodesWidth;	//X
	int m_numNodesLength;	//Z

	btScalar m_gridCellWidth;
	btScalar m_gridCellWidthInv;
	
	int m_rIndex;								///Index of m_r[], toggled every frame
	
	btAlignedObjectArray<btScalar> m_temp;
	btAlignedObjectArray<btScalar> m_height;	///Combined height; fluid + ground
	btAlignedObjectArray<btScalar> m_ground;	///Heightfield
	btAlignedObjectArray<btScalar> m_eta; 		///Depth of fluid; height - ground
	btAlignedObjectArray<btScalar> m_u;			///x velocity
	btAlignedObjectArray<btScalar> m_v;			///z velocity
	btAlignedObjectArray<btScalar> m_r[2];		///Displaced fluid
	btAlignedObjectArray<btScalar> m_fillRatio;
	btAlignedObjectArray<bool> m_flags;			///If false, the cell is a 'dry' cell with no fluid
	
	inline int getIndex(int i, int j) const
	{
		btAssert (i >= 0);
		btAssert (i < m_numNodesWidth);
		btAssert (j >= 0);
		btAssert (j < m_numNodesLength);
		int index = i + (j * m_numNodesWidth);
		return index;
	}
	
	void setGridDimensions(btScalar gridCellWidth, int numNodesWidth, int numNodesLength)
	{
		m_gridCellWidth = gridCellWidth;
		m_gridCellWidthInv = btScalar(1.0) / gridCellWidth;
		
		m_numNodesWidth = numNodesWidth;
		m_numNodesLength = numNodesLength;

		//
		m_rIndex = 0;
		
		//
		int numNodes = numNodesWidth * numNodesLength;
		m_temp.resize(numNodes);
		m_height.resize(numNodes);
		m_ground.resize(numNodes);
		m_eta.resize(numNodes);
		m_u.resize(numNodes);
		m_v.resize(numNodes);
		m_r[0].resize(numNodes);
		m_r[1].resize(numNodes);
		m_fillRatio.resize(numNodes);
		m_flags.resize(numNodes);

		for (int i = 0; i < numNodes; i++)
		{
			m_temp[i] = btScalar(0.0);
			m_height[i] = btScalar(0.0);
			m_eta[i] = btScalar(0.0);
			m_u[i] = btScalar(0.0);
			m_v[i] = btScalar(0.0);
			m_r[0][i] = btScalar(0.0);
			m_r[1][i] = btScalar(0.0);
			m_ground[i] = btScalar(0.0);
			m_fillRatio[i] = btScalar(0.0);
			m_flags[i] = false;
		}
	}
};

#endif