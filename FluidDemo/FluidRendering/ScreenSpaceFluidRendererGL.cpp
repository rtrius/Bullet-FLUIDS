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
		normal.xy = gl_TexCoord[0].xy*2.0 - 1.0;
		float r2 = dot(normal.xy, normal.xy);
		if(r2 > 1.0) discard;
		normal.z = sqrt(1.0 - r2);
		
		vec4 pixelPosition = vec4(eyePosition + normal*pointRadius, 1.0);
		vec4 clipSpacePosition = projectionMatrix * pixelPosition;

		float depth = clipSpacePosition.z / clipSpacePosition.w;
		
		//normal.z decreases as the distance from the sphere center increases
		float thickness = normal.z * 0.05;	
		
		gl_FragColor = vec4( vec3(1.0), thickness );
		gl_FragDepth = depth;
	}
);

const char fullScreenTextureVertexShader[] = STRINGIFY(
	void main()
	{
		gl_Position = gl_ModelViewProjectionMatrix * vec4(gl_Vertex.xyz, 1.0);
		gl_TexCoord[0] = gl_MultiTexCoord0;
		gl_FrontColor = gl_Color;
	}
);

const char bilateralFilter1dFragmentShader_depth[] = STRINGIFY(
	uniform float texelSize;
	uniform float filterRadiusPixels;
	uniform float blurScale;			//blurScale: Lower values increase blur
	uniform float blurDepthFalloff;		//blurDepthFalloff: Higher values decrease blurring between pixels of differing intensity
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
const char bilateralFilter1dFragmentShader_alpha[] = STRINGIFY(
	uniform float texelSize;
	uniform float filterRadiusPixels;
	uniform float blurScale;			//blurScale: Lower values increase blur
	uniform float blurDepthFalloff;		//blurDepthFalloff: Higher values decrease blurring between pixels of differing intensity
	uniform vec2 blurDirection;
	uniform sampler2D alphaTexture;
	void main()
	{
		float depth = texture(alphaTexture, gl_TexCoord[0]).a;
		
		float sum = 0.0;
		float wsum = 0.0;
		for(float x = -filterRadiusPixels; x <= filterRadiusPixels; x += 1.0)
		{	
			float neighborDepth = texture(alphaTexture, gl_TexCoord[0].xy + blurDirection*texelSize*x).a;

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
		gl_FragColor = vec4( vec3(1.0), sum );
	}
);

const char absorptionAndTransparencyFragmentShader[] = STRINGIFY(

	//Beer's law / absorption constants(xyz == rgb)
	//Controls the darkening of the fluid's color based on its thickness
	//For a constant k, (k > 1) == darkens faster; (k < 1) == darkens slower; (k == 0) == disable
	uniform vec3 absorption;
	
	uniform sampler2D thicknessTexture;
	
	void main()
	{
		float thickness = texture(thicknessTexture, gl_TexCoord[0]).a;
		gl_FragColor = vec4( exp(-absorption.x * thickness), exp(-absorption.y * thickness), exp(-absorption.z * thickness), thickness );
	}
);

const char generateSurfaceFragmentShader[] = STRINGIFY(
	vec3 getEyePos(sampler2D depthTexture, vec2 texCoord, mat4 projectionMatrix)
	{
		float depth = texture(depthTexture, texCoord).x;
		
		vec4 unprojectedPosition = inverse(projectionMatrix) * vec4(texCoord.xy*2.0 - 1.0, depth, 1.0);
		
		vec3 eyePosition = unprojectedPosition.xyz / unprojectedPosition.w;
		return eyePosition;
	}
	
	uniform mat4 depthProjectionMatrix;		//Projection matrix used to generate depth values
	uniform vec2 texelSize;
	uniform sampler2D depthTextureBlurred;
	uniform sampler2D depthTextureNotBlurred;
	uniform sampler2D absorptionAndTransparencyTexture;
	void main()
	{
		const float MAX_DEPTH = 0.99;
		float depthBlurred = texture(depthTextureBlurred, gl_TexCoord[0]).x;
		if(depthBlurred >= MAX_DEPTH) 
		{
			discard;
			return;
		}
		
		//Calculate normal using texCoords, depth, and projection matrix 
		vec3 eyePosition = getEyePos(depthTextureBlurred, gl_TexCoord[0].xy, depthProjectionMatrix);

		vec3 ddx = getEyePos(depthTextureBlurred, gl_TexCoord[0].xy + vec2(texelSize.x, 0), depthProjectionMatrix) - eyePosition;
		vec3 ddx2 = eyePosition - getEyePos(depthTextureBlurred, gl_TexCoord[0].xy + vec2(-texelSize.x, 0), depthProjectionMatrix);
		if( abs(ddx.z) > abs(ddx2.z) ) ddx = ddx2;

		vec3 ddy = getEyePos(depthTextureBlurred, gl_TexCoord[0].xy + vec2(0, texelSize.y), depthProjectionMatrix) - eyePosition;
		vec3 ddy2 = eyePosition - getEyePos(depthTextureBlurred, gl_TexCoord[0].xy + vec2(0, -texelSize.y), depthProjectionMatrix);
		if( abs(ddy.z) > abs(ddy2.z) ) ddy = ddy2;
		
		vec3 normal = normalize( cross(ddx, ddy) );
		
		//
		const vec3 LIGHT_DIRECTION = vec3(0.577, 0.577, 0.577);
		const float SHININESS = 40.0;
		
		float diffuse = max( 0.0, dot(normal, LIGHT_DIRECTION) );
		
		vec3 v = normalize(-eyePosition);			//Normalized vector pointing at camera/eye
		vec3 h = normalize(LIGHT_DIRECTION + v);	//Normalized vector halfway between LIGHT_DIRECTION and v
		float specular = pow( max(0.0, dot(normal, h)), SHININESS );
		
		//
		vec4 absorptionAndTransparency = texture(absorptionAndTransparencyTexture, gl_TexCoord[0]);
		
		const float MINIMUM_ALPHA = 0.80;
		vec3 color = gl_Color.xyz * absorptionAndTransparency.xyz * diffuse + specular;
		float alpha = MINIMUM_ALPHA + absorptionAndTransparency.w * (1.0 - MINIMUM_ALPHA);
		
		gl_FragColor = vec4(color, alpha);
		//gl_FragColor = vec4(gl_Color.xyz * diffuse + specular, 1.0);
		
		//Convert depth from Normalized Device Coordinates(NDC) to Window/Screen coordinates
		float depth = min( depthBlurred, texture(depthTextureNotBlurred, gl_TexCoord[0]).x );
		gl_FragDepth = (gl_DepthRange.diff*depth + gl_DepthRange.near + gl_DepthRange.far) * 0.5;
		
		//
		const int DISPLAY_FLUID = 0;
		const int DISPLAY_DEPTH = 1;
		const int DISPLAY_LINEAR_DEPTH = 2;
		const int DISPLAY_NORMAL = 3;
		const int DISPLAY_THICKNESS = 4;
		const int DISPLAY_ABSORPTION = 5;
		
		const int DISPLAY_MODE = DISPLAY_FLUID;
		switch(DISPLAY_MODE)
		{
			case DISPLAY_DEPTH:
				gl_FragColor = vec4( vec3(depthBlurred), 1.0 );
				break;
			case DISPLAY_LINEAR_DEPTH:	
				float linearDepth = min(1.0, -eyePosition.z / 10000.0);		//10000.0 == far z value
				gl_FragColor = vec4( vec3(linearDepth), 1.0 );
				break;
			case DISPLAY_NORMAL:
				gl_FragColor = vec4( (normal + 1.0) * 0.5, 1.0 );
				break;
			case DISPLAY_THICKNESS:
				float thickness = absorptionAndTransparency.a;
				gl_FragColor = vec4( vec3(thickness), 1.0 );
				break;
			case DISPLAY_ABSORPTION:
				gl_FragColor = vec4(absorptionAndTransparency.xyz, 1.0);
				break;
				
			case DISPLAY_FLUID:
			default:
				//Fluid color is set above
				break;
		}
	}
);

const char blitFragmentShader[] = STRINGIFY(
	uniform sampler2D rgbaTexture;
	uniform sampler2D depthTexture;
	void main()
	{
		gl_FragColor = texture(rgbaTexture, gl_TexCoord[0]);
		gl_FragDepth = texture(depthTexture, gl_TexCoord[0]).x;
	}
);


ScreenSpaceFluidRendererGL::ScreenSpaceFluidRendererGL(int screenWidth, int screenHeight)
{
	initializeGlew();
	
	//
	m_generateDepthShader = compileProgram(generateDepthVertexShader, generateDepthFragmentShader);
	m_blurDepthShader = compileProgram(fullScreenTextureVertexShader, bilateralFilter1dFragmentShader_depth);
	
	m_blurThickShader = compileProgram(fullScreenTextureVertexShader, bilateralFilter1dFragmentShader_alpha);
	m_absorptionAndTransparencyShader = compileProgram(fullScreenTextureVertexShader, absorptionAndTransparencyFragmentShader);
	
	m_generateSurfaceShader = compileProgram(fullScreenTextureVertexShader, generateSurfaceFragmentShader);
	m_blitShader = compileProgram(fullScreenTextureVertexShader, blitFragmentShader);
	
	//
	glGenBuffers(1, &m_positionVertexBuffer);
	
	//
	m_frameBuffer.initialize(screenWidth, screenHeight);
	m_tempColorTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_RGBA_TEXTURE);
	m_tempDepthTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	
	m_depthTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	m_blurredDepthTexturePass1 = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	m_blurredDepthTexturePass2 = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
	
	m_thickTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_ALPHA_TEXTURE);
	m_blurredThickTexturePass1 = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_ALPHA_TEXTURE);
	m_blurredThickTexturePass2 = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_ALPHA_TEXTURE);
	m_absorptionAndTransparencyTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_RGBA_TEXTURE);
	
	m_surfaceColorTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_RGBA_TEXTURE);
	m_surfaceDepthTexture = m_frameBuffer.createFrameBufferTexture(FrameBufferGL::FBT_DEPTH_TEXTURE);
}
ScreenSpaceFluidRendererGL::~ScreenSpaceFluidRendererGL()
{
	glDeleteShader(m_generateDepthShader);
	glDeleteShader(m_blurDepthShader);
	
	glDeleteShader(m_blurThickShader);
	glDeleteShader(m_absorptionAndTransparencyShader);
	
	glDeleteShader(m_generateSurfaceShader);
	glDeleteShader(m_blitShader);
	
	//
	glDeleteBuffers(1, &m_positionVertexBuffer);
	
	//
	m_frameBuffer.deactivate();
}

void ScreenSpaceFluidRendererGL::render(const btAlignedObjectArray<btVector3> &particlePositions, float sphereRadius, 
										float r, float g, float b, float absorptionR, float absorptionG, float absorptionB)
{
	btAssert( sizeof(btVector3) == 16 );

	glDepthMask(GL_TRUE);
	glEnable(GL_DEPTH_TEST);
	
	render_stage1_generateDepthTexture(particlePositions, sphereRadius);
	
	render_stage2_blurDepthTexture();
	
	render_stage3_generateThickTexture(particlePositions, sphereRadius);
	
	render_stage4_blurThickTexture();
	
	render_stage5_generateAbsorptionAndTransparencyTexture(absorptionR, absorptionG, absorptionB);
	
	glColor4f(r, g, b, 1.0f);
	render_stage6_generateSurfaceTexture();
	glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
	
	//Blit results to the main/window frame buffer
	glEnable(GL_BLEND);
	glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
	glUseProgram(m_blitShader);
	glUniform1i( glGetUniformLocation(m_blitShader, "rgbaTexture"), 0 );
	glUniform1i( glGetUniformLocation(m_blitShader, "depthTexture"), 1 );
		renderFullScreenTexture(m_surfaceColorTexture, m_surfaceDepthTexture, 0);
	glDisable(GL_BLEND);
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
	glGetFloatv(GL_PROJECTION_MATRIX, m_depthProjectionMatrix); 	//Used to reconstruct positions from depth values

	glEnable(GL_POINT_SPRITE);
	glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
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
	glDisable(GL_POINT_SPRITE);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}
void ScreenSpaceFluidRendererGL::render_stage2_blurDepthTexture()
{	
	glDepthFunc(GL_ALWAYS);
	glUseProgram(m_blurDepthShader);
	
	//First pass blurs along the x-axis
	glUniform1f( glGetUniformLocation(m_blurDepthShader, "texelSize"), 1.0f / static_cast<float>( m_frameBuffer.getWidth() ) );
	glUniform1f( glGetUniformLocation(m_blurDepthShader, "filterRadiusPixels"), 64.0f );
	glUniform1f( glGetUniformLocation(m_blurDepthShader, "blurScale"), 0.1f );
	glUniform1f( glGetUniformLocation(m_blurDepthShader, "blurDepthFalloff"), 0.6f );
	glUniform2f( glGetUniformLocation(m_blurDepthShader, "blurDirection"), 1.0f, 0.0f );
	glUniform1i( glGetUniformLocation(m_blurDepthShader, "depthTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_tempColorTexture, m_blurredDepthTexturePass1);
		renderFullScreenTexture(m_depthTexture, 0, 0);		
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	//Second pass blurs along the y-axis
	glUniform1f( glGetUniformLocation(m_blurDepthShader, "texelSize"), 1.0f / static_cast<float>( m_frameBuffer.getHeight() ) );
	glUniform1f( glGetUniformLocation(m_blurDepthShader, "filterRadiusPixels"), 64.0f );
	glUniform1f( glGetUniformLocation(m_blurDepthShader, "blurScale"), 0.1f );
	glUniform1f( glGetUniformLocation(m_blurDepthShader, "blurDepthFalloff"), 0.6f );
	glUniform2f( glGetUniformLocation(m_blurDepthShader, "blurDirection"), 0.0f, 1.0f );
	glUniform1i( glGetUniformLocation(m_blurDepthShader, "depthTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_tempColorTexture, m_blurredDepthTexturePass2);
		renderFullScreenTexture(m_blurredDepthTexturePass1, 0, 0);
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	glDepthFunc(GL_LESS);
	glUseProgram(0);
}
void ScreenSpaceFluidRendererGL::render_stage3_generateThickTexture(const btAlignedObjectArray<btVector3> &particlePositions, float sphereRadius)
{	
	glEnable(GL_POINT_SPRITE);
	glTexEnvi(GL_POINT_SPRITE, GL_COORD_REPLACE, GL_TRUE);
	glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
	glUseProgram(m_generateDepthShader);

	glDepthFunc(GL_ALWAYS);
	glEnable(GL_BLEND);
	glBlendFunc(GL_ONE, GL_ONE);
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
		
		m_frameBuffer.attachAndSetRenderTargets(m_thickTexture, m_tempDepthTexture);
			glDrawArrays(GL_POINTS, 0, numParticles);
		m_frameBuffer.detachAndUseDefaultFrameBuffer();
		
		glDisableClientState(GL_VERTEX_ARRAY);
	}
	glDepthFunc(GL_LESS);
	glDisable(GL_BLEND);
	
	glUseProgram(0);
	glDisable(GL_POINT_SPRITE);
	glBindBuffer(GL_ARRAY_BUFFER, 0);
}
void ScreenSpaceFluidRendererGL::render_stage4_blurThickTexture()
{	
	glDepthFunc(GL_ALWAYS);
	glUseProgram(m_blurThickShader);
	
	//First pass blurs along the x-axis
	glUniform1f( glGetUniformLocation(m_blurThickShader, "texelSize"), 1.0f / static_cast<float>( m_frameBuffer.getWidth() ) );
	glUniform1f( glGetUniformLocation(m_blurThickShader, "filterRadiusPixels"), 64.0f );
	glUniform1f( glGetUniformLocation(m_blurThickShader, "blurScale"), 0.1f );
	glUniform1f( glGetUniformLocation(m_blurThickShader, "blurDepthFalloff"), 0.6f );
	glUniform2f( glGetUniformLocation(m_blurThickShader, "blurDirection"), 1.0f, 0.0f );
	glUniform1i( glGetUniformLocation(m_blurThickShader, "alphaTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_blurredThickTexturePass1, m_tempDepthTexture);
		renderFullScreenTexture(m_thickTexture, 0, 0);		
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	//Second pass blurs along the y-axis
	glUniform1f( glGetUniformLocation(m_blurThickShader, "texelSize"), 1.0f / static_cast<float>( m_frameBuffer.getHeight() ) );
	glUniform1f( glGetUniformLocation(m_blurThickShader, "filterRadiusPixels"), 64.0f );
	glUniform1f( glGetUniformLocation(m_blurThickShader, "blurScale"), 0.1f );
	glUniform1f( glGetUniformLocation(m_blurThickShader, "blurDepthFalloff"), 0.6f );
	glUniform2f( glGetUniformLocation(m_blurThickShader, "blurDirection"), 0.0f, 1.0f );
	glUniform1i( glGetUniformLocation(m_blurThickShader, "alphaTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_blurredThickTexturePass2, m_tempDepthTexture);
		renderFullScreenTexture(m_blurredThickTexturePass1, 0, 0);
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	glDepthFunc(GL_LESS);
	glUseProgram(0);
}
void ScreenSpaceFluidRendererGL::render_stage5_generateAbsorptionAndTransparencyTexture(float absorptionR, float absorptionG, float absorptionB)
{
	glDepthFunc(GL_ALWAYS);
	glUseProgram(m_absorptionAndTransparencyShader);
	
	glUniform3f( glGetUniformLocation(m_absorptionAndTransparencyShader, "absorption"), absorptionR, absorptionG, absorptionB);
	glUniform1i( glGetUniformLocation(m_absorptionAndTransparencyShader, "thicknessTexture"), 0 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_absorptionAndTransparencyTexture, m_tempDepthTexture);
		renderFullScreenTexture(m_blurredThickTexturePass2, 0, 0);
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	glDepthFunc(GL_LESS);
	glUseProgram(0);
}

void ScreenSpaceFluidRendererGL::render_stage6_generateSurfaceTexture()
{
	float texelSize_x = 1.0f / static_cast<float>( m_frameBuffer.getWidth() );
	float texelSize_y = 1.0f / static_cast<float>( m_frameBuffer.getHeight() );

	glUseProgram(m_generateSurfaceShader);
	glUniformMatrix4fv( glGetUniformLocation(m_generateSurfaceShader, "depthProjectionMatrix"), 1, false, m_depthProjectionMatrix );
	glUniform2f( glGetUniformLocation(m_generateSurfaceShader, "texelSize"), texelSize_x, texelSize_y );
	glUniform1i( glGetUniformLocation(m_generateSurfaceShader, "depthTextureBlurred"), 0 );
	glUniform1i( glGetUniformLocation(m_generateSurfaceShader, "depthTextureNotBlurred"), 1 );
	glUniform1i( glGetUniformLocation(m_generateSurfaceShader, "absorptionAndTransparencyTexture"), 2 );
	
	m_frameBuffer.attachAndSetRenderTargets(m_surfaceColorTexture, m_surfaceDepthTexture);
		//renderFullScreenTexture(m_depthTexture, m_absorptionAndTransparencyTexture, 0);
		renderFullScreenTexture(m_blurredDepthTexturePass2, m_depthTexture, m_absorptionAndTransparencyTexture);
	m_frameBuffer.detachAndUseDefaultFrameBuffer();
	
	glUseProgram(0);
}

void ScreenSpaceFluidRendererGL::renderFullScreenTexture(GLuint texture2d_0, GLuint texture2d_1, GLuint texture2d_2)
{
	//Enable states
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
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, texture2d_2);
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
	glDisable(GL_TEXTURE_2D);
	glDisableClientState(GL_VERTEX_ARRAY);
	glDisableClientState(GL_TEXTURE_COORD_ARRAY);
	
	//Reset matrices
	glMatrixMode(GL_MODELVIEW);	
	glPopMatrix();
	glMatrixMode(GL_PROJECTION);
	glPopMatrix();
	
	//
	glActiveTexture(GL_TEXTURE2);
	glBindTexture(GL_TEXTURE_2D, 0);
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
		char *stringStart = &string[0];
		
		glGetProgramInfoLog(program, MAX_STRING_LENGTH, 0, stringStart);
		printf("GL Program Build Log:\n");
		printf("%s\n", stringStart);
		
		glGetShaderInfoLog(vertexShader, MAX_STRING_LENGTH, 0, stringStart);
		printf("Vertex Shader Build Log:\n");
		printf("%s\n", stringStart);
		
		glGetShaderInfoLog(fragmentShader, MAX_STRING_LENGTH, 0, stringStart);
		printf("Fragment Shader Build Log:\n");
		printf("%s\n", stringStart);
		
	
		glDeleteProgram(program);
		program = 0;
	}
	
	//If program was compiled successfully, marks shaders for deletion(not deleted until program is deleted)
	glDeleteShader(vertexShader);
	glDeleteShader(fragmentShader);
	
	return program;
}


