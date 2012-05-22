
  
    
function createDemos( demos, incdirs, linknames)
	for _, name in ipairs(demos) do
				
			project ( "App_" .. name )
			
			kind "ConsoleApp"
			targetdir ".."
	
	  	includedirs {incdirs}
	  
		configuration { "Windows" }
			defines { "GLEW_STATIC"}
	 		links { "opengl32" }
			includedirs{	"../Glut"	}
	 		libdirs {"../Glut"}
	 		files   { "../msvc/bullet.rc" }
	 		
	 		configuration {"Windows", "x32"}
				links {"glew32s","glut32"}
			configuration {"Windows", "x64"}
				links {"glew64s", "glut64"}
	
		configuration {"MacOSX"}
			--print "hello"
	 		linkoptions { "-framework Carbon -framework OpenGL -framework AGL -framework Glut" } 
		
		configuration {"not Windows", "not MacOSX"}
			links {"GL","GLU","glut"}
		
		configuration{}
	
		links { 
				linknames
			}
		
		files     { 
		 	"./" .. name .. "/*.cpp" ,
		 	"./" .. name .. "/*.h"
		 }
	end
end

-- "CharacterDemo", fixme: it includes BspDemo files

 local localdemos = {
    "BasicDemo",
    "Box2dDemo",
    "BspDemo",
    "CcdPhysicsDemo",
    "CollisionDemo",
    "CollisionInterfaceDemo",
    "ConcaveConvexcastDemo",
    "ConcaveDemo",
    "ConcaveRaycastDemo",
    "ConstraintDemo",
    "ContinuousConvexCollision",
    "ConvexHullDistance",
    "DynamicControlDemo",
    "EPAPenDepthDemo",
    "ForkLiftDemo",
    "FractureDemo",
    "GenericJointDemo",
    "GimpactTestDemo",
    "GjkConvexCastDemo",
    "HelloWorld",
    "InternalEdgeDemo",
    "MovingConcaveDemo",
    "MultiMaterialDemo",
    "RagdollDemo",
    "Raytracer",
    "SimplexDemo",
    "SliderConstraintDemo",
    "TerrainDemo",
    "UserCollisionAlgorithm",
    "VehicleDemo",
    "VoronoiFractureDemo"
  }

-- the following demos require custom include or link settings

 createDemos(localdemos,{"../src","OpenGL"},{"OpenGLSupport","LinearMath","BulletCollision","BulletDynamics"})
 
 createDemos({"ConvexDecompositionDemo"},{"../Extras/HACD","../Extras/ConvexDecomposition","../src","OpenGL"},{"OpenGLSupport","LinearMath","BulletCollision","BulletDynamics", "HACD","ConvexDecomposition"})
 
 createDemos({"SoftDemo"},{"../src","OpenGL"}, {"OpenGLSupport","LinearMath","BulletCollision","BulletDynamics", "BulletSoftBody"})
 
 createDemos({"SerializeDemo"},{"../Extras/Serialize/BulletFileLoader","../Extras/Serialize/BulletWorldImporter","../src","OpenGL"},{"OpenGLSupport","LinearMath","BulletCollision","BulletDynamics", "BulletSoftBody", "BulletFileLoader","BulletWorldImporter"})
 

------------------------------- FluidDemo --------------------------------------
function createFluidDemo(demoname, incdirs, linknames)
	
	USING_OPENCL = true
	USING_AMD_OPENCL = true			--If false, then use NVIDIA OpenCL
	
	if (USING_OPENCL) then
		if (USING_AMD_OPENCL) then
			hasCL = findOpenCL_AMD()
		else
			hasCL = findOpenCL_NVIDIA()
		end
		
		if (not hasCL) then 
			return
		end
	end
	
		project ( "App_" .. demoname )
		
		if (USING_OPENCL) then
			if (USING_AMD_OPENCL) then
				initOpenCL_AMD()
			else
				initOpenCL_NVIDIA()
			end
		end
		
		kind "ConsoleApp"
		targetdir ".."

  	includedirs { incdirs }
  
	configuration { "Windows" }
		defines { "GLEW_STATIC"}
 		links { "opengl32" }
		includedirs{	"../Glut"	}
 		libdirs {"../Glut"}
 		files   { "../msvc/bullet.rc" }
 		
 		configuration {"Windows", "x32"}
			links {"glew32s", "glut32"}
		configuration {"Windows", "x64"}
			links {"glew64s", "glut64"}

	configuration {"MacOSX"}
		--print "hello"
 		linkoptions { "-framework Carbon -framework OpenGL -framework AGL -framework Glut" } 
	
	configuration {"not Windows", "not MacOSX"}
		links {"GL", "GLU", "glut"}
	
	configuration{}

	links 
	{ 
		linknames
	}
	
	if (USING_OPENCL) then
		files    
		{ 
			"./" .. demoname .. "/*.cpp",
			"./" .. demoname .. "/*.h",
			"./" .. demoname .. "/Fluids/*.cpp",
			"./" .. demoname .. "/Fluids/*.h",
			"./" .. demoname .. "/Fluids/OpenCL_support/*.cpp",
			"./" .. demoname .. "/Fluids/OpenCL_support/*.h"
		}
	else
		files    
		{ 
			"./" .. demoname .. "/*.cpp",
			"./" .. demoname .. "/*.h",
			"./" .. demoname .. "/Fluids/*.cpp",
			"./" .. demoname .. "/Fluids/*.h",
		}
	end
	
end

include_dirs = {"../src", "OpenGL"}
linked_libs = {"OpenGLSupport", "LinearMath", "BulletCollision", "BulletDynamics"}
createFluidDemo("FluidDemo", include_dirs, linked_libs)
------------------------------- FluidDemo --------------------------------------
 
   
include "OpenGL"
 
 
