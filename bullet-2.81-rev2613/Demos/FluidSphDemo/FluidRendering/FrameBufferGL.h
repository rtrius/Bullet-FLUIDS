/* FrameBufferGL.h
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
#ifndef FRAME_BUFFER_GL_H
#define FRAME_BUFFER_GL_H

#include "LinearMath/btAlignedObjectArray.h"

#ifndef __APPLE__
#include <GL/glew.h>
#endif

class FrameBufferGL
{
public:
	enum FramebufferTextureType
	{
		FBT_RGBA_TEXTURE,
		FBT_DEPTH_TEXTURE,
		FBT_ALPHA_TEXTURE
	};
	
private:
	static const GLuint DEFAULT_FRAMEBUFFER = 0;	//Window manager frame buffer
	static const GLint MIPMAP_LEVEL = 0;
	
	GLuint m_frameBuffer;
	
	int m_textureWidth;
	int m_textureHeight;
	
	btAlignedObjectArray<GLuint> m_depthTextures;
	btAlignedObjectArray<GLuint> m_colorTextures;
	
public:
	void initialize(int width, int height);
	void deactivate();
	
	void attachAndSetRenderTargets(GLuint out_colorTexture, GLuint out_depthTexture);
	void detachAndUseDefaultFrameBuffer();
	
	int getWidth() const { return m_textureWidth; }
	int getHeight() const { return m_textureHeight; }
	
	void resizeTextures(int width, int height);
	
	GLuint createFrameBufferTexture(FramebufferTextureType type);	//Created textures are managed/deleted by FrameBufferGL
	
private:
	static GLuint internalCreateFrameBufferTexture(int width, int height, FramebufferTextureType type);
};

#endif
