workspace "SCESS"
	architecture "x64"

	configurations
	{
		"Debug",
		"Release"
	}

	filter{"system:linux"}
		defines{"OS_LINUX"}

	outputdir = "%{os.host()}/%{cfg.buildcfg}-x64/"

project "QuantumSCESS"
	kind "SharedLib"
	language "C++"
	cppdialect "C++20"
	systemversion "latest"
	location ("./Projects/" .. _ACTION)
	
	targetdir("bin/" .. outputdir)
	objdir("bin-obj/" .. outputdir)
	debugdir("bin/" .. outputdir)
	
	pchheader("./src/pch.h")
	pchsource("./src/pch.cpp")

	
	files
	{
	  "./src/Quantum/**.h",
	  "./src/Quantum/**.hpp",
	  "./src/Quantum/**.cpp",
	  "./src/Utilities/**.h",
	  "./src/Utilities/**.hpp",
	  "./src/Utilities/**.cpp",
	  "./src/Resources/**.cpp",
	  "./src/Resources/**.hpp",
	  "./src/Resources/**.h"
	}
	
	defines{ "EIGEN_USE_BLAS", "EIGEN_USE_LAPACKE" }
	includedirs
	{
	  "./src/",
	  "./Libraries/GraphicsEngine-1.0.0/include/",
	   "./Libraries/glm",
	  "./Libraries/cppitertools/",
	  "./Libraries/JSON/include/",
	  "./Libraries/Eigen/",
	  "./Libraries/LibInt2/include",
	  "./Libraries/LibXC/lib/include"
	}
	
	libdirs
	{
	    "./Libraries/GraphicsEngine-1.0.0/",
	    "./Libraries/LibInt2/",
	    "./Libraries/LibXC/lib"
	}
	links {"lapack", "lapacke", "blas", "boost_filesystem", "GraphicsEngine", "libint2.a", "libxc.a" }
	
	filter {  "configurations:Debug" }
		defines { "CONFIG_DEBUG", "_DEBUG" }
		optimize "Off"
		symbols "On"
		
	filter { "configurations:Release" }
		defines { "CONFIG_RELEASE" }
		optimize "Speed"
		vectorextensions "AVX2"
		symbols "Off"
	
	
project "SCESS"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++20"
	systemversion "latest"
	location ("./Projects/" .. _ACTION)
	
	targetdir("bin/" .. outputdir)
	objdir("bin-obj/" .. outputdir)
	debugdir("bin/" .. outputdir)
	
	files
	{
	  "./src/SCESS/**.h",
	  "./src/SCESS/**.hpp",
	  "./src/SCESS/**.cpp"
	}

	includedirs
	{
	  "./src/",
	  "./Libraries/GraphicsEngine-1.0.0/include/",
	  "./Libraries/glm",
	  "./Libraries/cppitertools/",
	  "./Libraries/JSON/include/",
	  "./Libraries/Eigen/",
	  "./Libraries/LibInt2/include",
 	  "./Libraries/LibInt2/include",
	  "./Libraries/LibXC/lib/include"
	}
	
	pchheader("./src/pch.h")
	pchsource("./src/pch.cpp")
	
	libdirs
	{
		"./Libraries/LibXC/lib",
		"./Libraries/GraphicsEngine-1.0.0/"
	}
	links { "QuantumSCESS", "xc.a", "GraphicsEngine" }
		
	filter {  "configurations:Debug" }
		defines { "CONFIG_DEBUG", "_DEBUG" }
		optimize "Off"
		symbols "On"
		
	filter { "configurations:Release" }
		defines { "CONFIG_RELEASE" }
		optimize "Speed"
		vectorextensions "AVX2"
		symbols "Off"

project "ViewerSCESS"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++20"
	systemversion "latest"
	location ("./Projects/" .. _ACTION)

	targetdir("bin/" .. outputdir)
	objdir("bin-obj/" .. outputdir)
	debugdir("bin/" .. outputdir)

	
	files
	{
	  "./src/Viewer/**.h",
	  "./src/Viewer/**.hpp",
	  "./src/Viewer/**.cpp",
	  "./src/Utilities/**.h",
	  "./src/Utilities/**.hpp",
	  "./src/Utilities/**.cpp",
	  "./src/Resources/**.cpp",
	  "./src/Resources/**.hpp",
	  "./src/Resources/**.h"
	}

	includedirs
	{
	  "./src/",
	  "./Libraries/GraphicsEngine-1.0.0/include/",
	  "./Libraries/glm",
	  "./Libraries/cppitertools/",
	  "./Libraries/JSON/include/",
	  "./Libraries/Eigen/",
	  "./Libraries/LibInt2/include",
	  "./Libraries/nfd/src/include/",
  	  "./Libraries/LibXC/lib/include"
	}
	
	pchheader("./src/pch.h")
	pchsource("./src/pch.cpp")
	
	libdirs
	{
	    "./Libraries/GraphicsEngine-1.0.0/",
	    "./Libraries/nfd/build/lib/Release/x64",
	    "./Libraries/LibXC/lib"
	}
	links { "QuantumSCESS", "GraphicsEngine", "xc.a", "nfd.a", "gtk-3", "gdk-3", "glib-2.0" }
		
	filter {  "configurations:Debug" }
		--libdirs { "../../C++/GraphicsEngine/bin/linux/Debug-x64/GraphicsEngine" }
		defines { "CONFIG_DEBUG", "_DEBUG" }
		optimize "Off"
		symbols "On"
		
	filter { "configurations:Release" }
		--libdirs { "../../C++/GraphicsEngine/bin/linux/Release-x64/GraphicsEngine" }
		defines { "CONFIG_RELEASE" }
		optimize "Speed"
		vectorextensions "AVX2"
		symbols "Off"
