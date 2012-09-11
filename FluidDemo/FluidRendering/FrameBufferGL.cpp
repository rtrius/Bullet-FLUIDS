/* FrameBufferGL.cpp
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
#include "FrameBufferGL.h"

#include <cstdio>

void FrameBufferGL::initialize(int width, int height)
{
	m_textureWidth = width;
	m_textureHeight = height;
	
	//
	glGenFramebuffers(1, &m_frameBuffer);
	glBindFramebuffer(GL_FRAMEBUFFER, m_frameBuffer);
	
	//
	glBindFramebuffer(GL_FRAMEBUFFER, DEFAULT_FRAMEBUFFER);
}
void FrameBufferGL::deactivate()
{
	//Detach textures
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, MIPMAP_LEVEL);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, MIPMAP_LEVEL);
	
	//
	glDeleteFramebuffers(1, &m_frameBuffer);
	
	//
	glDeleteTextures( m_depthTextures.size(), &m_depthTextures[0] );
	m_depthTextures.resize(0);
	
	glDeleteTextures( m_colorTextures.size(), &m_colorTextures[0] );
	m_colorTextures.resize(0);
	
}

void FrameBufferGL::attachAndSetRenderTargets(GLuint out_colorTexture, GLuint out_depthTexture)
{
	glBindFramebuffer(GL_FRAMEBUFFER, m_frameBuffer);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, out_colorTexture, MIPMAP_LEVEL);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, out_depthTexture, MIPMAP_LEVEL);
	
	GLenum frameBufferStatus = glCheckFramebufferStatus(GL_FRAMEBUFFER);
	if(frameBufferStatus != GL_FRAMEBUFFER_COMPLETE)
		printf("Error: OpenGL FrameBuffer is not properly configured(%d).\n", frameBufferStatus);
	
	//
	glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
	glClearDepth(1.0f);
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
}
void FrameBufferGL::detachAndUseDefaultFrameBuffer()
{
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, 0, MIPMAP_LEVEL);
	glFramebufferTexture2D(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_TEXTURE_2D, 0, MIPMAP_LEVEL);
	glBindFramebuffer(GL_FRAMEBUFFER, DEFAULT_FRAMEBUFFER);
}

void FrameBufferGL::resizeTextures(int width, int height)
{
	m_textureWidth = width;
	m_textureHeight = height;
	
	for(int i = 0; i < m_depthTextures.size(); ++i)
	{
		glBindTexture(GL_TEXTURE_2D, m_depthTextures[i]);
		glTexImage2D(GL_TEXTURE_2D, MIPMAP_LEVEL, GL_DEPTH_COMPONENT32, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	}
	
	for(int i = 0; i < m_colorTextures.size(); ++i)
	{
		glBindTexture(GL_TEXTURE_2D, m_colorTextures[i]);
		glTexImage2D(GL_TEXTURE_2D, MIPMAP_LEVEL, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
	}
	
	glBindTexture(GL_TEXTURE_2D, 0);
}

GLuint FrameBufferGL::createFrameBufferTexture(FramebufferTextureType type)
{
	GLuint textureId = internalCreateFrameBufferTexture(m_textureWidth, m_textureHeight, type);
	
	if(type == FBT_DEPTH_TEXTURE) m_depthTextures.push_back(textureId);
	else m_colorTextures.push_back(textureId);
	
	return textureId;
}
GLuint FrameBufferGL::internalCreateFrameBufferTexture(int width, int height, FramebufferTextureType type)
{
	GLuint textureId;
	glGenTextures(1, &textureId);
	glBindTexture(GL_TEXTURE_2D, textureId);
	
	if(type == FBT_DEPTH_TEXTURE)
		glTexImage2D(GL_TEXTURE_2D, MIPMAP_LEVEL, GL_DEPTH_COMPONENT32, width, height, 0, GL_DEPTH_COMPONENT, GL_FLOAT, 0);
	else
		glTexImage2D(GL_TEXTURE_2D, MIPMAP_LEVEL, GL_RGBA, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, 0);
	
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
	
	glBindTexture(GL_TEXTURE_2D, 0);
	return textureId;
}
