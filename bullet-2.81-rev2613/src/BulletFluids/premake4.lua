	project "BulletFluids"
		
	kind "StaticLib"
	targetdir "../../lib"
	includedirs {
		"..",
	}
	files {
		"**.cpp",
		"**.h"
	}