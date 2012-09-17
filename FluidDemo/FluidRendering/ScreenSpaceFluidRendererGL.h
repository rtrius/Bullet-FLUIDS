/* ScreenSpaceFluidRendererGL.h
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
#ifndef SCREEN_SPACE_FLUID_RENDERER_GL_H
#define SCREEN_SPACE_FLUID_RENDERER_GL_H

#include <GL/glew.h>

#include "LinearMath/btVector3.h"
#include "LinearMath/btAlignedObjectArray.h"

#include "FrameBufferGL.h"


//ScreenSpaceFluidRendererGL constructor initializes GLEW; will fail if an OpenGL context does not exist
class ScreenSpaceFluidRendererGL
{
	GLuint m_positionVertexBuffer;
	GLfloat m_depthProjectionMatrix[16];	//Projection matrix when m_generateDepthShader is run
	
	GLuint m_generateDepthShader;
	GLuint m_blurDepthShader;
	
	GLuint m_blurThickShader;
	GLuint m_absorptionAndTransparencyShader;
	
	GLuint m_generateSurfaceShader;
	GLuint m_blitShader;
	
	FrameBufferGL m_frameBuffer;
	
	//Textures are managed(and deleted) by m_frameBuffer
	GLuint m_tempColorTexture;
	GLuint m_tempDepthTexture;
	
	GLuint m_depthTexture;
	GLuint m_blurredDepthTexturePass1;
	GLuint m_blurredDepthTexturePass2;
	
	GLuint m_thickTexture;
	GLuint m_blurredThickTexturePass1;
	GLuint m_blurredThickTexturePass2;
	GLuint m_absorptionAndTransparencyTexture;
	
	GLuint m_surfaceColorTexture;
	GLuint m_surfaceDepthTexture;
	
public:
	ScreenSpaceFluidRendererGL(int screenWidth, int screenHeight);
	~ScreenSpaceFluidRendererGL();

	void render(const btAlignedObjectArray<btVector3> &particlePositions, float sphereRadius, 
				float r, float g, float b, float absorptionR, float absorptionG, float absorptionB);
	
	void setResolution(int width, int height) { m_frameBuffer.resizeTextures(width, height); }
	
private:
	void initializeGlew();
	
	void render_stage1_generateDepthTexture(const btAlignedObjectArray<btVector3> &particlePositions, float sphereRadius);
	void render_stage2_blurDepthTexture();
	void render_stage3_generateThickTexture(const btAlignedObjectArray<btVector3> &particlePositions, float sphereRadius);
	void render_stage4_blurThickTexture();
	void render_stage5_generateAbsorptionAndTransparencyTexture(float absorptionR, float absorptionG, float absorptionB);
	void render_stage6_generateSurfaceTexture();
	
	void renderFullScreenTexture(GLuint texture2d_0, GLuint texture2d_1, GLuint texture2d_2);
	
	static GLuint compileProgram(const char *vertexShaderSource, const char *fragmentShaderSource);
};


#endif
