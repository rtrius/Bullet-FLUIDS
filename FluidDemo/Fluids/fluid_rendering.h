/** fluid_rendering.h
	Copyright (C) 2012 Jackson Lee

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
#ifndef FLUID_RENDERING_H_INCLUDED
#define FLUID_RENDERING_H_INCLUDED

#include "fluid_system.h"

#include "GlutStuff.h"
#include "LinearMath/btVector3.h"


///The indicies of a MarchingCube must correspond to this arrangement;
///that is, the position of vertices[i] relative to other vertices in
///the cube should match the index in the diagram.
///		
///		(certain edges removed in the 2nd cube for clarity)
///		
///		   4-------5           4-------5
///		  /|      /|           |       |
/// 	 / |     / |           |       |
///		7-------6  |        7-------6  |
///		|  0----|--1        |  0----|--1
///		| /     | /         |       | /
///		|/      |/          |       |/
///		3-------2           3-------2
struct MarchingCube
{
	btVector3 m_vertices[8];
	float m_scalars[8];
};

extern const int edgeTable[256];
extern const int triangleTable[256][16];
class MarchingCubes
{
	int m_cellsPerEdge;
	btAlignedObjectArray<float> m_scalarField;	//Scalar field of dimension m_cellsPerEdge^3
	
	btAlignedObjectArray<float> m_vertices;
	
public:
	MarchingCubes() : m_cellsPerEdge(0) {}
	
	void initialize(int cellsPerEdge)
	{
		m_cellsPerEdge = cellsPerEdge;
		m_scalarField.resize(cellsPerEdge * cellsPerEdge * cellsPerEdge);
		for(int i = 0; i < m_scalarField.size(); ++i)m_scalarField[i] = 0.f;
	}
	
	void generateMesh(const FluidSystem &FS)
	{
		BT_PROFILE("MarchingCubes::generateMesh()");
	
		btVector3 cellSize;
		loadScalarField(FS, m_cellsPerEdge, &m_scalarField, &cellSize);
		
		m_vertices.resize(0);
		marchingCubes(FS.getGrid().getParameters().m_min, cellSize, m_scalarField, m_cellsPerEdge, &m_vertices);
	}
	
	const btAlignedObjectArray<float>& getTriangleVertices() const { return m_vertices; }
	
private:
	static void loadScalarField(const FluidSystem &FS, int cellsPerEdge, btAlignedObjectArray<float> *out_scalarField, btVector3 *out_cellSize);
	static void marchingCubes(const btVector3 &gridMin, const btVector3 &cellSize, const btAlignedObjectArray<float> &scalarField, 
							  int cellsPerEdge, btAlignedObjectArray<float> *out_vertices);
	static void generateVertices(const MarchingCube &C, btAlignedObjectArray<float> *out_vertices);
	
	static inline btVector3 getVertex(const btVector3 &gridMin, const btVector3 &cellSize,
									  int index_x, int index_y, int index_z)
	{
		return btVector3 ( 	gridMin.x() + cellSize.x() * static_cast<float>(index_x), 
							gridMin.y() + cellSize.y() * static_cast<float>(index_y),
							gridMin.z() + cellSize.z() * static_cast<float>(index_z) );
	}

	static inline btVector3 vertexLerp(float isolevel, float scalar1, float scalar2, 
									   const btVector3 &p1, const btVector3 &p2)
	{
		//P = P1 + (isolevel - V1) (P2 - P1) / (V2 - V1) 	
		//P == vertex, V == scalar at that vertex

		float scalar = (isolevel - scalar1) / (scalar2 - scalar1);

		return btVector3( p1.x() + scalar * (p2.x() - p1.x()),
						  p1.y() + scalar * (p2.y() - p1.y()),
						  p1.z() + scalar * (p2.z() - p1.z())  );
	}
};

GLuint generateSphereList(float radius);
void drawSphere(GLuint glSphereList, const btVector3 &position, float velocity);

#endif



