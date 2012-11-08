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
#include "HfFluidDemo_GL_ShapeDrawer.h"

#include "LinearMath/btIDebugDraw.h"

#include "GlutStuff.h"

#include "BulletHfFluid/btHfFluidBuoyantConvexShape.h"
#include "BulletHfFluid/btHfFluid.h"
#include "BulletHfFluid/btHfFluidCollisionShape.h"

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

void HfFluidDemo_GL_ShapeDrawer::drawOpenGL(btScalar* m, const btCollisionShape* shape, const btVector3& color, 
											int debugMode, const btVector3& worldBoundsMin, const btVector3& worldBoundsMax)
{
	if( shape->getShapeType() == HFFLUID_BUOYANT_CONVEX_SHAPE_PROXYTYPE )
	{
		const btConvexShape* convexShape = static_cast<const btHfFluidBuoyantConvexShape*>(shape)->getConvexShape();
		GL_ShapeDrawer::drawOpenGL(m, convexShape, color, debugMode, worldBoundsMin, worldBoundsMax);
		
		return;
	}

	if (shape->getShapeType() == HFFLUID_SHAPE_PROXYTYPE)
	{
		glPushMatrix(); 
		btglMultMatrix(m);
		
		//glEnable(GL_BLEND);
		//glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA); 
		
		glColor4f(0.3f, 0.5f, 1.0f, 0.95f);
		
			const btHfFluidCollisionShape* hfFluidShape = static_cast<const btHfFluidCollisionShape*>(shape);
			const btHfFluid* fluid = hfFluidShape->m_fluid;
			GlDrawcallback drawCallback;
			drawCallback.m_wireframe = (debugMode & btIDebugDraw::DBG_DrawWireframe) != 0;
			fluid->forEachSurfaceTriangle(&drawCallback, worldBoundsMin, worldBoundsMax);
			
			
		glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		
		//glDisable(GL_BLEND);
		
		glPopMatrix();
		
		return;
	}
	
	
	GL_ShapeDrawer::drawOpenGL(m, shape, color, debugMode, worldBoundsMin, worldBoundsMax);
}

