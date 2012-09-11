/* ScreenSpaceFluidRendererGL.cpp
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
#include "ScreenSpaceFluidRendererGL.h"

#include <cstdlib> 	//exit()
#include <cstdio>


#define STRINGIFY(A) #A

const char generateDepthVertexShader[] = STRINGIFY(
	uniform float pointRadius;  	//Point size in world space
	uniform float pointScale;  	 	//Scale to calculate size in pixels
	varying vec3 eyePosition;
	varying mat4 projectionMatrix;
	void main()
	{
		 //Calculate window-space point size
		eyePosition = ( gl_ModelViewMatrix * vec4(gl_Vertex.xyz, 1.0) ).xyz;
		float distance = length(eyePosition);
		gl_PointSize = pointRadius * (pointScale / distance);

		projectionMatrix = gl_ProjectionMatrix;
		
		gl_TexCoord[0] = gl_MultiTexCoord0;
		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);
		gl_FrontColor = gl_Color;
	}
);
const char generateDepthFragmentShader[] = STRINGIFY(
	uniform float pointRadius;  	//Point size in world space
	varying vec3 eyePosition;       //Position of sphere center in eye space
	varying mat4 projectionMatrix;
	void main()
	{
		vec3 normal;
		normal.xy = gl_TexCoord[0].xy * vec2(2.0, 2.0) - vec2(1.0, 1.0);
		float r2 = dot(normal.xy, normal.xy);
		if(r2 > 1.0) discard;
		normal.z = sqrt(1.0 - r2);
		
		vec4 pixelPosition = vec4(eyePosition + normal*pointRadius, 1.0);
		vec4 clipSpacePosition = projectionMatrix * pixelPosition;

		float depth = clipSpacePosition.z / clipSpacePosition.w;
		
		gl_FragColor = vec4(normal, 1.0);
		gl_FragDepth = depth;
	}
);

const char bilateralFilter1dVertexShader[] = STRINGIFY(
	void main()
	{
		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);
		gl_TexCoord[0] = gl_MultiTexCoord0;
	}
);
const char bilateralFilter1dFragmentShader[] = STRINGIFY(
	uniform float texelSize;
	uniform float filterRadiusPixels;
	uniform float blurScale;			//blurScale: Lower values increase blur
	uniform float blurDepthFalloff;		//blurDepthFalloff: Higher values increase 'sharpness'
	uniform vec2 blurDirection;
	uniform sampler2D depthTexture;
	void main()
	{
		float depth = texture(depthTexture, gl_TexCoord[0]).x;
		
		float sum = 0.0;
		float wsum = 0.0;
		for(float x = -filterRadiusPixels; x <= filterRadiusPixels; x += 1.0)
		{	
			float neighborDepth = texture(depthTexture, gl_TexCoord[0].xy + blurDirection*texelSize*x).x;

			//Spatial domain
			float r = x * blurScale;
			float w = exp(-r*r);

			//Range domain
			float r2 = (neighborDepth - depth) * blurDepthFalloff;
			float g = exp(-r2*r2);
			
			sum += neighborDepth * w * g;
			wsum += w * g;
		}

		if(wsum > 0.0) sum /= wsum;
		gl_FragDepth = sum;
	}
);

const char generateSurfacelVertexShader[] = STRINGIFY(
	void main()
	{
		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);
		gl_TexCoord[0] = gl_MultiTexCoord0;
		gl_FrontColor = gl_Color;
	}
);
const char generateSurfaceVertexShader[] = STRINGIFY(
	vec3 getEyePos(sampler2D depthTexture, vec2 texCoord, mat4 projectionMatrix)
	{
		float depth = texture(depthTexture, texCoord).x;
		
		vec4 unprojectedPosition = inverse(projectionMatrix) * vec4(texCoord.xy*vec2(2.0, 2.0) - vec2(1.0, 1.0), depth, 1.0);
		
		vec3 eyePosition = unprojectedPosition.xyz / unprojectedPosition.w;
		return eyePosition;
	}
	
	uniform mat4 depthProjectionMatrix;		//Projection matrix used to generate depth values
	uniform vec2 texelSize;
	uniform sampler2D depthTexture;
	void main()
	{
		const float MAX_DEPTH = 0.99;
		float depth = texture(depthTexture, gl_TexCoord[0]).x;
		if(depth >= MAX_DEPTH) 
		{
			discard;
			return;
		}
		
		//Calculate normal using texCoords, depth, and projection matrix 
		vec3 eyePosition = getEyePos(depthTexture, gl_TexCoord[0].xy, depthProjectionMatrix);

		vec3 ddx = getEyePos(depthTexture, gl_TexCoord[0].xy + vec2(texelSize.x, 0), depthProjectionMatrix) - eyePosition;
		vec3 ddx2 = eyePosition - getEyePos(depthTexture, gl_TexCoord[0].xy + vec2(-texelSize.x, 0), depthProjectionMatrix);
		if( abs(ddx.z) > abs(ddx2.z) ) ddx = ddx2;

		vec3 ddy = getEyePos(depthTexture, gl_TexCoord[0].xy + vec2(0, texelSize.y), depthProjectionMatrix) - eyePosition;
		vec3 ddy2 = eyePosition - getEyePos(depthTexture, gl_TexCoord[0].xy + vec2(0, -texelSize.y), depthProjectionMatrix);
		if( abs(ddy.z) > abs(ddy2.z) ) ddy = ddy2;
		
		vec3 normal = normalize( cross(ddx, ddy) );
		
		//
		const vec3 LIGHT_DIRECTION = vec3(0.577, 0.577, 0.577);
		const float SHININESS = 40.0;
		
		float diffuse = max( 0.0, dot(normal, LIGHT_DIRECTION) );
		
		vec3 v = normalize(-eyePosition);			//Normalized vector pointing at camera/eye
		vec3 h = normalize(LIGHT_DIRECTION + v);	//Normalized vector halfway between LIGHT_DIRECTION and v
		float specular = pow( max(0.0, dot(normal, h)), SHININESS );
		
		gl_FragColor = vec4(gl_Color.xyz * diffuse + specular, 1.0);
		gl_FragDepth = depth;
	
		//
		const int DISPLAY_DEPTH = 0;
		if(DISPLAY_DEPTH)
		{
			gl_FragColor = vec4(depth, depth, depth, 1.0);
		}
		const int DISPLAY_LINEAR_DEPTH = 0;
		if(DISPLAY_LINEAR_DEPTH)
		{
			float linearDepth = min(1.0, -eyePosition.z / 10000.0);		//10000.0 == far z value
			gl_FragColor = vec4(linearDepth, linearDepth, linearDepth, 1.0);
		}
		const int DISPLAY_NORMAL = 0;
		if(DISPLAY_NORMAL)
		{
			gl_FragColor = vec4( (normal + 1.0) * 0.5, 1.0 );
		}
	}
);

const char blitVertexShader[] = STRINGIFY(
	void main()
	{
		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);
		gl_TexCoord[0] = gl_MultiTexCoord0;
	}
);
const char blitFragmentShader[] = STRINGIFY(
	uniform sampler2D rgbaTexture;
	uniform sampler2D depthTexture;
	void main()
	{
		gl_FragColor = texture(rgbaTexture, gl_TexCoord[0]);
		gl_FragDepth = texture(depthTexture, gl_TexCoord[0]);
	}
);


ScreenSpaceFluidRendererGL::ScreenSpaceFluidRendererGL(int screenWidth, int screenHeight)
{
	initializeGlew();
	
	//
	m_generateDepthShader = compileProgram(generateDepthVertexShader, generateDepthFragmentShader);
	m_blurShader = compileProgram(bilateralFilter1dVertexShader, bilateralFilter1dFragmentShader);
	m_generateNormalShader = compileProgram(generateSurfacelVertexShader, generateSurfaceVertexShader);
	m_blitShader = compileProgram(blitVertexShader, blitFragmentShader);
	
	//
	glGenBuffers(1, &m_positionVertexBuffer);
	
	//
	m_frameBuffer.initialize(screenWidth, screenHeight);
	m_tempColorTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_RGBA_TEXTURE);
	m_surfaceDepthTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	
	m_depthTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	m_blurredDepthTexturePass1 = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	m_blurredDepthTexturePass2 = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	m_surfaceColorTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_RGBA_TEXTURE);
}
ScreenSpaceFluidRendererGL::~ScreenSpaceFluidRendererGL()
{
	glDeleteShader(m_generateDepthShader);
	glDeleteShader(m_blurShader);
	glDeleteShader(m_generateNormalShader);
	glDeleteShader(m_blitShader);
	
	//
	glDeleteBuffers(1, &m_positionVertexBuffer);
	
	//
	m_frameBuffer.deactivate();
}

void ScreenSpaceFluidRendererGL::render(const btAlignedObjectArray<btVector3> &particlePositions, float sphereRadius, float r, float g, float b)
{
	render_stage1_generateDepthTexture(particlePositions, sphereRadius);
	
	render_stage2_blurDepthTexture();
	
	glColor4f(r, g, b, 1.0f);
	render_stage3_generateSurfaceTexture();
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	
	//Blit results to the main/window frame buffer
	glUseProgram(m_blitShader);
	glUniform1i( glGetUniformLocation(m_blitShader, "rgbaTexture"), 0 );
	glUniform1i( glGetUniformLocation(m_blitShader, "depthTexture"), 1 );
		renderFullScreenTexture(m_surfaceColorTexture, m_surfaceDepthTexture);
	glUseProgram(0);
	
	//Default clear color for Bullet demos, set in DemoApplication::myinit()
	glClearColor(0.7f, 0.7f, 0.7f, 0.0f);
}
	
void ScreenSpaceFluidRendererGL::initializeGlew()
{
	GLenum errorCode = glewInit();
	if(errorCode != GLEW_OK)
	{
		printf( "GLEW error: %s(%d) \n", glewGetErrorString(errorCode), errorCode );
	}
	
	const int NUM_REQUIRED_EXTENSIONS = 5;
	const char *requiredExtensions[NUM_REQUIRED_EXTENSIONS] =
	{
		//"GL_VERSION_4_2",
		//"GL_VERSION_3_0",
		"GL_VERSION_2_0",
		"GL_ARB_multitexture",
		"GL_ARB_vertex_buffer_object",
		"GL_ARB_framebuffer_object",
		"GL_ARB_point_sprite"
	};
	
	bool areRequiredExtensionsMissing = false;
	for(int i = 0; i < NUM_REQUIRED_EXTENSIONS; ++i)
	{
		if( !glewIsSupported(requiredExtensions[i]) )
		{
			printf("Required extension: %s is missing.\n", requiredExtensions[i]);
			areRequiredExtensionsMissing = true;
		}
	}
	
	if(areRequiredExtensionsMissing) 
	{
		printf("ScreenSpaceFluidRendererGL: Required OpenGL extensions missing.\n");
		printf("Press enter to exit.\n");
		getchar();
		exit(-1);
	}
}

void ScreenSpaceFluidRendererGL::render_stage1_generateDepthTexture(const btAlignedObjectArray<btVector3> &particlePositions, float sphereRadius)
{
	btAssert( sizeof(btVector3) == 16 );

	glGetFloatv(GL_PROJECTION_MATRIX, m_depthProjectionMatrix); 	//Used to reconstruct positions from depth values

	glEnable(GL_POINT_SPRITE);
	glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);
	glUseProgram(m_generateDepthShader);

	{
		int numParticles = particlePositions.size();
		
		float screenWidth = static_cast<float>( m_frameBuffer.getWidth() );
		float screenHeight = static_cast<float>( m_frameBuffer.getHeight() );
		float lesserDistance = (screenWidth > screenHeight) ?  screenHeight : screenWidth;
		glUniform1f( glGetUniformLocation(m_generateDepthShader, "pointScale"), lesserDistance );
		glUniform1f( glGetUniformLocation(m_generateDepthShader, "pointRadius"), sphereRadius );

		glBindBuffer(GL_ARRAY_BUFFER, m_positionVertexBuffer);
		glBufferData(GL_ARRAY_BUFFER, sizeof(btVector3) * numParticles, &particlePositions[0], GL_DYNAMIC_DRAW);
		glVertexPointer(4, GL_FLOAT, 0, 0);
		glEnableClientState(GL_VERTEX_ARRAY);
		
		m_frameBuffer.attachAndSetRenderTargets(m_tempColorTexture, m_depthTexture);
			glDrawArrays(GL_POINTS, 0, numParticles);
		m_frameBuffer.detachAndUseDefaultFrameBuffer();
		
		glDisableClientState(GL_VERTEX_ARRAY);
	}
	
	glUseProgram(0);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_POINT_SPRITE_ARB);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}
void ScreenSpaceFluidRendererGL::render_stage2_blurDepthTexture()
{	
	glDepthFunc(GL_ALWAYS);
	glUseProgram(m_blurShader);
	
	//First pass blurs along the x-axis
	glUniform1f( glGetUniformLocation(m_blurShader, "texelSize"), 1.0f / static_cast<float>( m_frameBuffer.getWidth() ) );
	glUniform1f( glGetUniformLocation(m_blurShader, "filterRadiusPixels"), 64.0f );
	glUniform1f( glGetUniformLocation(m_blurShader, "blurScale"), 0.1f );
	glUniform1f( glGetUniformLocation(m_blurShader, "blurDepthFalloff"), 0.6f );
	glUniform2f( glGetUniformLocation(m_blurShader, "blurDirection"), 1.0f, 0.0f );
	glUniform1i( glGetUniformLocation(m_blurShader, "depthTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_tempColorTexture, m_blurredDepthTexturePass1);
		renderFullScreenTexture(m_depthTexture, 0);		
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	//Second pass blurs along the y-axis
	glUniform1f( glGetUniformLocation(m_blurShader, "texelSize"), 1.0f / static_cast<float>( m_frameBuffer.getHeight() ) );
	glUniform1f( glGetUniformLocation(m_blurShader, "filterRadiusPixels"), 64.0f );
	glUniform1f( glGetUniformLocation(m_blurShader, "blurScale"), 0.1f );
	glUniform1f( glGetUniformLocation(m_blurShader, "blurDepthFalloff"), 0.6f );
	glUniform2f( glGetUniformLocation(m_blurShader, "blurDirection"), 0.0f, 1.0f );
	glUniform1i( glGetUniformLocation(m_blurShader, "depthTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_tempColorTexture, m_blurredDepthTexturePass2);
		renderFullScreenTexture(m_blurredDepthTexturePass1, 0);
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	glDepthFunc(GL_LESS);
	glUseProgram(0);
}
void ScreenSpaceFluidRendererGL::render_stage3_generateSurfaceTexture()
{
	float texelSize_x = 1.0f / static_cast<float>( m_frameBuffer.getWidth() );
	float texelSize_y = 1.0f / static_cast<float>( m_frameBuffer.getHeight() );

	glUseProgram(m_generateNormalShader);
	glUniformMatrix4fv(  glGetUniformLocation(m_generateNormalShader, "depthProjectionMatrix"), 1, false, m_depthProjectionMatrix );
	glUniform2f( glGetUniformLocation(m_generateNormalShader, "texelSize"), texelSize_x, texelSize_y );
	glUniform1i( glGetUniformLocation(m_generateNormalShader, "depthTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_surfaceColorTexture, m_surfaceDepthTexture);
		//renderFullScreenTexture(m_depthTexture, 0);
		renderFullScreenTexture(m_blurredDepthTexturePass2, 0);
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	glUseProgram(0);
}

void ScreenSpaceFluidRendererGL::renderFullScreenTexture(GLuint texture2d_0, GLuint texture2d_1)
{
	//Enable states
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glEnable(GL_DEPTH_TEST);
	glEnable(GL_TEXTURE_2D);
	glEnableClientState(GL_VERTEX_ARRAY);
	glEnableClientState(GL_TEXTURE_COORD_ARRAY);
	
	//Configure matrices
	glMatrixMode(GL_MODELVIEW);	
	glPushMatrix();
	glLoadIdentity();
	
	glMatrixMode(GL_PROJECTION);	
	glPushMatrix();
	glLoadIdentity();
	glOrtho(0.0f, 1.0f, 0.0f, 1.0f, -1.0f, 1.0f);
	
	//
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, texture2d_1);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, texture2d_0);
	
	//Render
	{
		//Arranged for GL_TRIANGLE_STRIP; order: Upper Left, Upper Right, Lower Left, Lower Right
		const GLfloat SQUARE_VERTICES[4 * 2] = { 0.0f,1.0f, 1.0f,1.0f, 0.0f,0.0f, 1.0f,0.0f };
	
		glVertexPointer(2, GL_FLOAT, 0, SQUARE_VERTICES);
		glTexCoordPointer(2, GL_FLOAT, 0, SQUARE_VERTICES);
		
		glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
	}
	
	//Disable states
	glDisable(GL_BLEND);
	glDisable(GL_DEPTH_TEST);
	glDisable(GL_TEXTURE_2D);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	
	//Reset matrices
	glMatrixMode(GL_MODELVIEW);	
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	
	//
	glActiveTexture(GL_TEXTURE1);
	glBindTexture(GL_TEXTURE_2D, 0);
	glActiveTexture(GL_TEXTURE0);
	glBindTexture(GL_TEXTURE_2D, 0);
}

 GLuint ScreenSpaceFluidRendererGL::compileProgram(const char *vertexShaderSource, const char *fragmentShaderSource)
{
	GLuint vertexShader = glCreateShader(GL_VERTEX_SHADER);
	GLuint fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);

	glShaderSource(vertexShader, 1, &vertexShaderSource, 0);
	glShaderSource(fragmentShader, 1, &fragmentShaderSource, 0);
	
	glCompileShader(vertexShader);
	glCompileShader(fragmentShader);

	//
	GLuint program = glCreateProgram();

	glAttachShader(program, vertexShader);
	glAttachShader(program, fragmentShader);

	glLinkProgram(program);

	GLint success = 0;
	glGetProgramiv(program, GL_LINK_STATUS, &success);

	if(!success) 
	{
		const int MAX_STRING_LENGTH = 65536;
		btAlignedObjectArray<char> string;
		string.resize(MAX_STRING_LENGTH);
		
		glGetProgramInfoLog(program, MAX_STRING_LENGTH, 0, &string[0]);
		printf("GL Program Build Log:\n");
		printf("%s\n", &string[0]);
		
		glGetShaderInfoLog(vertexShader, MAX_STRING_LENGTH, 0, &string[0]);
		printf("Vertex Shader Build Log:\n");
		printf("%s\n", &string[0]);
		
		glGetShaderInfoLog(fragmentShader, MAX_STRING_LENGTH, 0, &string[0]);
		printf("Fragment Shader Build Log:\n");
		printf("%s\n", &string[0]);
		
	
		glDeleteProgram(program);
		program = 0;
	}
	
	//If program was compiled successfully, marks shaders for deletion(not deleted until program is deleted)
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
	
	return program;
}


