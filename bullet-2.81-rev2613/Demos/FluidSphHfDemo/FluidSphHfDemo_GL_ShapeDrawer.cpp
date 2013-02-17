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

#include "FluidSphHfDemo_GL_ShapeDrawer.h"

#include "LinearMath/btIDebugDraw.h"

#include "GlutStuff.h"

#include "BulletCollision/CollisionShapes/btSphereShape.h"

#include "BulletFluids/Hf/btFluidHfBuoyantConvexShape.h"
#include "BulletFluids/Hf/btFluidHf.h"
#include "BulletFluids/Hf/btFluidHfCollisionShape.h"
#include "BulletFluids/Sph/btFluidSph.h"
#include "BulletFluids/Sph/btFluidSphCollisionShape.h"

class GlDrawcallback : public btTriangleCallback
{
public:
	bool m_wireframe;

	GlDrawcallback() : m_wireframe(false) {}

	virtual void processTriangle(btVector3* triangle,int partId, int triangleIndex)
	{
		(void)triangleIndex;
		(void)partId;

		if (m_wireframe)
		{
			glBegin(GL_LINES);
			glColor3f(1, 0, 0);
			glVertex3d(triangle[0].getX(), triangle[0].getY(), triangle[0].getZ());
			glVertex3d(triangle[1].getX(), triangle[1].getY(), triangle[1].getZ());
			glColor3f(0, 1, 0);
			glVertex3d(triangle[2].getX(), triangle[2].getY(), triangle[2].getZ());
			glVertex3d(triangle[1].getX(), triangle[1].getY(), triangle[1].getZ());
			glColor3f(0, 0, 1);
			glVertex3d(triangle[2].getX(), triangle[2].getY(), triangle[2].getZ());
			glVertex3d(triangle[0].getX(), triangle[0].getY(), triangle[0].getZ());
			glEnd();
		} else
		{
			glBegin(GL_TRIANGLES);
			//glColor3f(1, 1, 1);
			
			glVertex3d(triangle[0].getX(), triangle[0].getY(), triangle[0].getZ());
			glVertex3d(triangle[1].getX(), triangle[1].getY(), triangle[1].getZ());
			glVertex3d(triangle[2].getX(), triangle[2].getY(), triangle[2].getZ());

			glVertex3d(triangle[2].getX(), triangle[2].getY(), triangle[2].getZ());
			glVertex3d(triangle[1].getX(), triangle[1].getY(), triangle[1].getZ());
			glVertex3d(triangle[0].getX(), triangle[0].getY(), triangle[0].getZ());
			glEnd();
		}
	}
};

inline void drawTriangle(const btVector3 &a, const btVector3 &b, const btVector3 &c)
{
	glVertex3f( a.x(), a.y(), a.z() );
	glVertex3f( b.x(), b.y(), b.z() );
	glVertex3f( c.x(), c.y(), c.z() );
}
inline void drawQuad(const btVector3 &rightReverse, const btVector3 &rightForward, const btVector3 &leftForward, const btVector3 &leftReverse)
{
	//drawTriangle(leftForward, rightForward, leftReverse);
	//drawTriangle(rightForward, rightReverse, leftReverse);
	drawTriangle(leftReverse, rightForward, leftForward);
	drawTriangle(leftReverse, rightReverse, rightForward);
}
void drawSolidBox(const btVector3 &min, const btVector3 &max)
{
	//LR - Left/Right 		(-x / +x)
	//TB - Top/Bottom 		(+y / -y)
	//fr - Forward/Reverse	(-z / +z)
	btVector3 LTf( min.x(), max.y(), min.z() );
	btVector3 LTr( min.x(), max.y(), max.z() );
	btVector3 LBf( min.x(), min.y(), min.z() );
	btVector3 LBr( min.x(), min.y(), max.z() );
	btVector3 RTf( max.x(), max.y(), min.z() );
	btVector3 RTr( max.x(), max.y(), max.z() );
	btVector3 RBf( max.x(), min.y(), min.z() );
	btVector3 RBr( max.x(), min.y(), max.z() );
	
	glBegin(GL_TRIANGLES);
		drawQuad(RTr, RTf, LTf, LTr);		//Top
		drawQuad(RBr, RBf, LBf, LBr);		//Bottom
		drawQuad(RTf, RBf, LBf, LTf);		//Forward
		drawQuad(RTr, RBr, LBr, LTr);		//Rear
		drawQuad(LTr, LBr, LBf, LTf);		//Left
		drawQuad(RTr, RBr, RBf, RTf);		//Right
	glEnd();

}
void drawHfFluidAsColumns(const btFluidHf* fluid)
{
	glDisable(GL_CULL_FACE);
	
	//const btVector3 &origin = fluid->getWorldTransform().getOrigin();
	
	for(int i = 0; i < fluid->getNumNodesX(); i++)
		for(int j = 0; j < fluid->getNumNodesZ(); j++)
		{
			int index = fluid->arrayIndex(i, j);

			if( (i % 2 && (j+1) % 2) || ((i+1) % 2 && j % 2) ) glColor4f(0.15f, 0.25f, 0.5f, 1.0f);
			else glColor4f(0.3f, 0.5f, 1.0f, 1.0f);
			
			btScalar h = fluid->getFluidHeight(index);
			btScalar g = fluid->getGroundHeight(index);

			const btScalar EPSILON = btScalar(0.03);
			//if (btFabs(h) < EPSILON) continue;
			//else if(h < -EPSILON)glColor4f(0.05f, 0.8f, 0.15f, 1.0f);
			if(h < EPSILON)glColor4f(0.05f, 0.8f, 0.15f, 1.0f);

			//btVector3 boxMin = origin + btVector3( fluid->getCellPosX(i), g, fluid->getCellPosZ(j) );
			//btVector3 boxMax = origin + btVector3( fluid->getCellPosX(i+1), g+h, fluid->getCellPosZ(j+1) );
			btVector3 boxMin = btVector3( fluid->getCellPosX(i), g, fluid->getCellPosZ(j) );
			btVector3 boxMax = btVector3( fluid->getCellPosX(i+1), g+h, fluid->getCellPosZ(j+1) );
			
			drawSolidBox(boxMin, boxMax);
		}
}
void drawHfGroundAsColumns(const btFluidHf* fluid)
{
	const btScalar MIN_HEIGHT = -100.0f;

	glDisable(GL_CULL_FACE);
	
	//const btVector3 &origin = fluid->getWorldTransform().getOrigin();
	
	for(int i = 0; i < fluid->getNumNodesX(); i++)
		for(int j = 0; j < fluid->getNumNodesZ(); j++)
		{
			if( (i % 2 && (j+1) % 2) || ((i+1) % 2 && j % 2) ) glColor4f(0.85f, 0.75f, 0.5f, 1.0f);
			else glColor4f(0.7f, 0.5f, 0.0f, 1.0f);
			
			int index = fluid->arrayIndex(i, j);
			btVector3 boxMin = btVector3( fluid->getCellPosX(i), MIN_HEIGHT, fluid->getCellPosZ(j) );
			btVector3 boxMax = btVector3( fluid->getCellPosX(i+1), fluid->getGroundHeight(index), fluid->getCellPosZ(j+1) );
			
			drawSolidBox(boxMin, boxMax);
		}
}

GLuint generateSphereList(float radius)
{
	//Sphere generation code from FLUIDS v.2
	GLuint glSphereList = glGenLists(1);
	glNewList(glSphereList, GL_COMPILE);
		glBegin(GL_TRIANGLE_STRIP);
			for(float tilt = -90.0f; tilt <= 90.0f; tilt += 10.0f) 
			{
				for(float ang = 0.f; ang <= 360.0f; ang += 30.0f) 
				{
					const float DEGREES_TO_RADIANS = 3.141592f/180.0f;
					
					float ang_radians = ang * DEGREES_TO_RADIANS;
					float tilt_radians = tilt * DEGREES_TO_RADIANS;
					float tilt1_radians = (tilt + 10.0f) * DEGREES_TO_RADIANS;
				
					float x = sin(ang_radians) * cos(tilt_radians);
					float y = cos(ang_radians) * cos(tilt_radians);
					float z = sin(tilt_radians);
					float x1 = sin(ang_radians) * cos(tilt1_radians);
					float y1 = cos(ang_radians) * cos(tilt1_radians);
					float z1 = sin(tilt1_radians);
					
					//Enabling normals causes the heightfield fluid to appear darker,
					//as the hfFluid does not specify normals.
					//glNormal3f(x, y, z);	
					glVertex3f(x*radius, y*radius, z*radius);		
					//glNormal3f(x1, y1, z1);	
					glVertex3f(x1*radius, y1*radius, z1*radius);
				}
			}
		glEnd();
	glEndList();
	
	return glSphereList;
}
inline void drawSphereWithList(GLuint glSphereList, const btVector3& position)
{	
	glPushMatrix();
		glTranslatef( position.x(), position.y(), position.z() );
		glCallList(glSphereList);
	glPopMatrix();
}

void FluidSphHfDemo_GL_ShapeDrawer::drawOpenGL(btScalar* m, const btCollisionShape* shape, const btVector3& color, 
											int debugMode, const btVector3& worldBoundsMin, const btVector3& worldBoundsMax)
{
	int shapeType = shape->getShapeType();
	
	switch(shapeType)
	{
		case HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE:
		{
			const btConvexShape* convexShape = static_cast<const btFluidHfBuoyantConvexShape*>(shape)->getConvexShape();
			GL_ShapeDrawer::drawOpenGL(m, convexShape, color, debugMode, worldBoundsMin, worldBoundsMax);
			
			break;
		}
		
		case HFFLUID_SHAPE_PROXYTYPE:
		{
			glPushMatrix(); 
			btglMultMatrix(m);
			
			GlDrawcallback drawCallback;
			drawCallback.m_wireframe = (debugMode & btIDebugDraw::DBG_DrawWireframe) != 0;
			
			const btFluidHfCollisionShape* hfFluidShape = static_cast<const btFluidHfCollisionShape*>(shape);
			const btFluidHf* fluid = hfFluidShape->m_fluid;
			
			glEnable(GL_BLEND);
			glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
			
			glColor4f(0.50f, 0.50f, 0.50f, 1.0f);
				if(m_drawHfGroundWithTriangles)fluid->forEachGroundTriangle(&drawCallback, worldBoundsMin, worldBoundsMax);
				
			glColor4f(0.3f, 0.5f, 1.0f, 0.8f);
				if(m_drawHfFluidWithTriangles)fluid->forEachSurfaceTriangle(&drawCallback, worldBoundsMin, worldBoundsMax);
				
			glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
			
			if(m_drawHfGroundAsColumns)drawHfGroundAsColumns(fluid);
			if(m_drawHfFluidAsColumns)drawHfFluidAsColumns(fluid);
			
			glDisable(GL_BLEND);
			
			glPopMatrix();
			
			break;
		}
		
		case INVALID_SHAPE_PROXYTYPE:	//	Temporarily use INVALID_SHAPE_PROXYTYPE for SPH fluids
		{
			const btFluidSphCollisionShape* sphFluidShape = static_cast<const btFluidSphCollisionShape*>(shape);
			const btFluidSph* fluid = sphFluidShape->m_owner;
		
			const btFluidSphParametersLocal& FL = fluid->getLocalParameters();
			
			static btScalar lastRadius = FL.m_particleRadius;
			static GLuint glSphereList = generateSphereList(FL.m_particleRadius);
		
			if(lastRadius != FL.m_particleRadius)
			{
				glDeleteLists(glSphereList, 1);
			
				lastRadius = FL.m_particleRadius;
				glSphereList = generateSphereList(FL.m_particleRadius);
			}
			
			glColor3f( 0.3f, 0.5f, 1.0f);
			for(int i = 0; i < fluid->numParticles(); ++i) 
				drawSphereWithList( glSphereList, fluid->getPosition(i) );
			
			break;
		}
		
		default:
			GL_ShapeDrawer::drawOpenGL(m, shape, color, debugMode, worldBoundsMin, worldBoundsMax);
			break;
	}
}

