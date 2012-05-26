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

#include <vector>

#include "vector3df.h"
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
	Vector3DF m_vertices[8];
	float m_scalars[8];
};

extern const int edgeTable[256];
extern const int triangleTable[256][16];
class MarchingCubes
{
	int m_cellsPerEdge;
	float *m_scalarField;	//Scalar field of dimension cellsPerEdge^3
	
	std::vector<float> m_vertices;
	
public:
	MarchingCubes() : m_scalarField(0) {}
	
	void initialize(int cellsPerEdge)
	{
		if(!m_scalarField)
		{
			m_cellsPerEdge = cellsPerEdge;
		
			int numCells = cellsPerEdge * cellsPerEdge * cellsPerEdge;
			
			m_scalarField = new float[numCells];
			//m_scalarField = new (std::nothrow) float[numCells]; 	//	requires #include <new>
			for(int i = 0; i < numCells; ++i)m_scalarField[i] = 0.f;
		}
		else
		{
			deactivate();
			initialize(cellsPerEdge);
		}
	}
	void deactivate()
	{
		if(m_scalarField) 
		{
			delete[] m_scalarField;
			m_scalarField = 0;
		}
	}
	
	void generateMesh(const FluidSystem &FS)
	{
		BT_PROFILE("MarchingCubes::generateMesh()");
	
		Vector3DF cellSize;
		loadScalarField(FS, m_scalarField, m_cellsPerEdge, &cellSize);
		
		m_vertices.clear();
		marchingCubes(FS.getGrid().getParameters().m_min, cellSize, m_scalarField, m_cellsPerEdge, &m_vertices);
	}
	
	const std::vector<float>& getTriangleVertices() const { return m_vertices; }
	
private:
	static void loadScalarField(const FluidSystem &FS, float *scalarField, int cellsPerEdge, Vector3DF *out_cellSize);
	static void marchingCubes(const Vector3DF &gridMin, const Vector3DF &cellSize,
							  float *scalarField, int cellsPerEdge, std::vector<float> *out_vertices);
	static void generateVertices(const MarchingCube &C, std::vector<float> *out_vertices);
	
	static inline Vector3DF getVertex(const Vector3DF &gridMin, const Vector3DF &cellSize,
									  int index_x, int index_y, int index_z)
	{
		return Vector3DF ( 	gridMin.x() + cellSize.x() * static_cast<float>(index_x), 
							gridMin.y() + cellSize.y() * static_cast<float>(index_y),
							gridMin.z() + cellSize.z() * static_cast<float>(index_z) );
	}

	static inline Vector3DF vertexLerp(float isolevel, float scalar1, float scalar2, 
									   const Vector3DF &p1, const Vector3DF &p2)
	{
		//P = P1 + (isolevel - V1) (P2 - P1) / (V2 - V1) 	
		//P == vertex, V == scalar at that vertex

		float scalar = (isolevel - scalar1) / (scalar2 - scalar1);

		return Vector3DF( p1.x() + scalar * (p2.x() - p1.x()),
						  p1.y() + scalar * (p2.y() - p1.y()),
						  p1.z() + scalar * (p2.z() - p1.z())  );
	}
};

GLuint generateSphereList(float radius);
void drawSphere(GLuint glSphereList, const Vector3DF &position, float velocity);

#endif



