workspace "SCESS"
	architecture "x64"

	configurations
	{
		"Debug",
		"Release"
	}

	filter{"system:linux"}
		defines{"OS_LINUX"}

	outputdir = "%{os.host()}/%{cfg.buildcfg}-x64/%{prj.name}"
	libDir = "%{wks.location}/libs/"

project "QuantumSCESS"
	kind "SharedLib"
	language "C++"
	cppdialect "C++20"
	systemversion "latest"
	
	targetdir("bin/" .. outputdir)
	objdir("bin-obj/" .. outputdir)
	debugdir("bin/" .. outputdir)

	
	files
	{
	  "%{wks.location}/src/Quantum/**.h",
	  "%{wks.location}/src/Quantum/**.hpp",
	  "%{wks.location}/src/Quantum/**.cpp",
	  "%{wks.location}/src/Utilities/**.h",
	  "%{wks.location}/src/Utilities/**.hpp",
	  "%{wks.location}/src/Utilities/**.cpp"
	}
	
	defines{ "EIGEN_USE_BLAS", "EIGEN_USE_LAPACKE" }
	includedirs
	{
	  "%{wks.location}/src/",
	  "/usr/include",
	  "/usr/local/include",
	  "../../C++/GraphicsEngine/src",
	  "../../Libraries/json-develop/include/",
	  "../../Libraries/eigen-3.4.90/",
          "../../Libraries/",
          "../../Libraries/cppitertools",
          "/opt/etsf/include"
	}
	
	pchheader("%{wks.location}/src/pch.h")
	pchsource("%{wks.location}/src/pch.cpp")
	
	libdirs
	{
	    "../../C++/GraphicsEngine/bin/linux/Release-x64/GraphicsEngine",
	    "../../Libraries/ImGui/bin/linux/Release-x64",
	    "/opt/etsf/lib/"
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
	
	targetdir("bin/" .. outputdir)
	objdir("bin-obj/" .. outputdir)
	debugdir("bin/" .. outputdir)
	
	files
	{
	  "%{wks.location}/src/SCESS/**.h",
	  "%{wks.location}/src/SCESS/**.hpp",
	  "%{wks.location}/src/SCESS/**.cpp",
	  "%{wks.location}/src/Utilities/**.h",
	  "%{wks.location}/src/Utilities/**.hpp",
	  "%{wks.location}/src/Utilities/**.cpp",
          "/opt/etsf/include"
	}

	includedirs
	{
	  "%{wks.location}/src/",
	  "/usr/include",
	  "/usr/local/include",
	  "/opt/etsf/include",
           "../../Libraries/eigen-3.4.90/",
           "../../C++/GraphicsEngine/src",
	}
	
	pchheader("%{wks.location}/src/pch.h")
	pchsource("%{wks.location}/src/pch.cpp")
	
	libdirs
	{  
		"/opt/etsf/lib/",
		"../../C++/GraphicsEngine/bin/linux/Release-x64/GraphicsEngine"
	}
	links { "QuantumSCESS" }
		
	filter {  "configurations:Debug" }
		defines { "CONFIG_DEBUG", "_DEBUG" }
		optimize "Off"
		symbols "On"
		
	filter { "configurations:Release" }
		defines { "CONFIG_RELEASE" }
		optimize "Speed"
		vectorextensions "AVX2"
		symbols "Off"

project "Explorer"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++20"
	systemversion "latest"

	targetdir("bin/" .. outputdir)
	objdir("bin-obj/" .. outputdir)
	debugdir("bin/" .. outputdir)

	
	files
	{
	  "%{wks.location}/src/Explorer/**.h",
	  "%{wks.location}/src/Explorer/**.hpp",
	  "%{wks.location}/src/Explorer/**.cpp",
	  "%{wks.location}/resource.h",
	  "%{wks.location}/SCESS.APS",
	  "%{wks.location}/SCESS.rc",
	  "%{wks.location}/src/Utilities/**.h",
	  "%{wks.location}/src/Utilities/**.hpp",
	  "%{wks.location}/src/Utilities/**.cpp"
	}

	includedirs
	{
	  "%{wks.location}/src/",
	  "/usr/include",
	  "/usr/local/include",
          "../../C++/GraphicsEngine/src",
	  "../../Libraries/json-develop/include",
	  "../../Libraries/eigen-3.4.90",
	  "../../Libraries/glm",
          "../../Libraries/",
	  "../../Libraries/ImGui",
          "../../Libraries/nativefiledialog/src/include",
          "/opt/etsf/include"
	}
	
	pchheader("%{wks.location}/src/pch.h")
	pchsource("%{wks.location}/src/pch.cpp")
	
	libdirs
	{
	    "../../C++/GraphicsEngine/bin/linux/Release-x64/GraphicsEngine",
	    "../../Libraries/ImGui/bin/linux/Release-x64",
	    "../../Libraries/nativefiledialog/build/lib/Release/x64",
	    "/opt/etsf/lib/"
	}
	links { "QuantumSCESS", "GraphicsEngine", "glfw3", "ImGui", "GLEW.a", "GL", "m", 
		"Xrandr", "Xi", "X11", "Xxf86vm", "pthread", "nfd", "gtk-3", "z", "gdk-3", "glib-2.0" }
		
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
