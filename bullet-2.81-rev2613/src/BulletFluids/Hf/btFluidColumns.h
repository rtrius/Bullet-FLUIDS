/*
Bullet Continuous Collision Detection and Physics Library
Copyright (c) 2003-2009 Erwin Coumans  http://bulletphysics.com

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising from the use of this software.
Permission is granted to anyone to use this software for any purpose, 
including commercial applications, and to alter it and redistribute it freely, 
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must not claim that you wrote the original software. If you use this software in a product, an acknowledgment in the product documentation would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.

Experimental Buoyancy fluid demo written by John McCutchan
*/
//This is an altered source version based on the HeightFieldFluidDemo included with Bullet Physics 2.80(bullet-2.80-rev2531).

#ifndef BT_FLUID_COLUMNS_H
#define BT_FLUID_COLUMNS_H

#include "LinearMath/btScalar.h"
#include "LinearMath/btAlignedObjectArray.h"

///Stores configuration of a single btFluidHf.
struct btFluidHfParameters
{
	///You can enforce a global velocity at the surface of the fluid; default: 0.0 and 0.0
	btScalar m_globalVelocityX;
	btScalar m_globalVelocityZ;
	
	///Force of gravity, should match physics world; must be below 0; default: -9.8
	btScalar m_gravity;
	
	///Controls the amount of fluid displaced into adjacent cells when a rigid body is submerged.
	///For example, a value of 2.0 means that the rigid body will displace twice its volume; default: 0.5.
	btScalar m_volumeDisplacementScale;
	
	///Influence of the fluid's horizontal velocity on submerged bodies; range [0.0, 1.0]; default: 0.5
	btScalar m_horizontalVelocityScale;

	btScalar m_heightEpsilon;
	btScalar m_dryCellEpsilon;		///<If btFluidColumns.m_fluidDepth is lower than this, the cell is considered empty/inactive
	
	btScalar m_density;				///<Density of the fluid when interacting with btRigidBody; kg / m^3.
	
	///@name Below parameters apply only to btFluidHfSolverExperimental
	///@{
	btScalar m_heightDamping; 		///<Scales the rate at which height changes(lower values increase viscosity); range (0.0, 1.0]; def: 1.0
	btScalar m_fluidFlowThreshold; 	///<A column must have at least this much fluid(height) to flow to another column; def: 0.0
	///@}
	
	btFluidHfParameters()
	{
		m_globalVelocityX = btScalar(0.0);
		m_globalVelocityZ = btScalar(0.0);
		m_gravity = btScalar(-9.8);

		m_volumeDisplacementScale = btScalar(0.25);
		m_horizontalVelocityScale = btScalar(0.5);
		
		m_heightEpsilon = btScalar(0.001);
		m_dryCellEpsilon = btScalar(0.01);
		
		//Avoid using realistic values, as fluid-rigid interaction becomes unstable with large mass ratios(water == 1000kg / m^3)
		//m_density = btScalar(0.35);
		m_density = btScalar(2.35);
		
		m_heightDamping = btScalar(1.0);
		m_fluidFlowThreshold = btScalar(0.0);
	}
};

///Coordinates the parallel arrays used to represent a heightfield fluid.
struct btFluidColumns
{
	int m_numNodesX;
	int m_numNodesZ;

	btScalar m_gridCellWidth;
	btScalar m_gridCellWidthInv;
	
	int m_displacedIndex;								///<Index of m_displaced[], toggled every frame
	
	btAlignedObjectArray<btScalar> m_temp;
	
	btAlignedObjectArray<btScalar> m_combinedHeight;	///<Combined height; fluid + ground
	btAlignedObjectArray<btScalar> m_ground;			///<Heightfield
	btAlignedObjectArray<btScalar> m_fluidDepth; 		///<Depth of fluid; combinedHeight - ground
	btAlignedObjectArray<btScalar> m_vel_x;				///<x velocity
	btAlignedObjectArray<btScalar> m_vel_z;				///<z velocity
	btAlignedObjectArray<btScalar> m_displaced[2];		///<Fluid displaced from rigid body interaction
	btAlignedObjectArray<btScalar> m_fillRatio;
	btAlignedObjectArray<bool> m_active;				///<If false, the cell is a 'dry' cell with no fluid and is not updated
	
	inline int getIndex(int i, int j) const
	{
		btAssert (i >= 0);
		btAssert (i < m_numNodesX);
		btAssert (j >= 0);
		btAssert (j < m_numNodesZ);
		int index = i + (j * m_numNodesX);
		return index;
	}
	
	int numColumns() const { return m_numNodesX * m_numNodesZ; }
	
	//Adjust dimensions of btFluidHf AABB after calling this
	void setGridDimensions(btScalar gridCellWidth, int numNodesWidth, int numNodesLength)
	{
		m_gridCellWidth = gridCellWidth;
		m_gridCellWidthInv = btScalar(1.0) / gridCellWidth;
		
		m_numNodesX = numNodesWidth;
		m_numNodesZ = numNodesLength;

		//
		m_displacedIndex = 0;
		
		//
		int numNodes = numNodesWidth * numNodesLength;
		m_temp.resize(numNodes);
		m_combinedHeight.resize(numNodes);
		m_ground.resize(numNodes);
		m_fluidDepth.resize(numNodes);
		m_vel_x.resize(numNodes);
		m_vel_z.resize(numNodes);
		m_displaced[0].resize(numNodes);
		m_displaced[1].resize(numNodes);
		m_fillRatio.resize(numNodes);
		m_active.resize(numNodes);

		for (int i = 0; i < numNodes; i++)
		{
			m_temp[i] = btScalar(0.0);
			m_combinedHeight[i] = btScalar(0.0);
			m_fluidDepth[i] = btScalar(0.0);
			m_vel_x[i] = btScalar(0.0);
			m_vel_z[i] = btScalar(0.0);
			m_displaced[0][i] = btScalar(0.0);
			m_displaced[1][i] = btScalar(0.0);
			m_ground[i] = btScalar(0.0);
			m_fillRatio[i] = btScalar(0.0);
			m_active[i] = false;
		}
	}
};

#endif