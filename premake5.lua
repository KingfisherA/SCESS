workspace "SCESS"
	architecture "x64"

	configurations
	{
		"Debug",
		"Release"
	}

filter{"system:linux"}
	defines{"OS_LINUX"}
	

filter{"system:windows"}
	defines{"OS_WIN"}

outputdir = "%{os.host()}/%{cfg.buildcfg}-x64/"

libDir = "%{wks.location}/libs/"

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
	  "%{wks.location}/src/**.h",
	  "%{wks.location}/src/**.hpp",
	  "%{wks.location}/src/**.cpp",
	  "%{wks.location}/resource.h",
	  "%{wks.location}/SCESS.APS",
	  "%{wks.location}/SCESS.rc",
	}

		
	filter {"system:linux"}
		defines{ "EIGEN_USE_BLAS", "EIGEN_USE_LAPACKE" }
		includedirs
		{
		  "%{wks.location}/src/",
		  "/usr/include",
		  "/usr/local/include",
		  "/opt/etsf/include",
                  "/usr/local/include/eigen3",
                  "../../C++/GraphicsEngine/src",
		  "../../Libraries/json-develop/include",
		  "../../Libraries/eigen-3.4.90",
		  "../../Libraries/glm",
                  "../../Libraries/",
		  "../../Libraries/ImGui",
                  "../../Libraries/nativefiledialog/src/include",
                  "../../Libraries/cppitertools"
		}
		libdirs
		{
		    "../../Libraries/ImGui/bin/linux/Release-x64",
		    "../../Libraries/nativefiledialog/build/lib/Release/x64"
		}
		links {"lapack", "lapacke", "blas", "GraphicsEngine", "glfw3", "ImGui", "libboost_filesystem.a", "GLEW.a", "GL", "m", "Xrandr", "Xi", "X11", "Xxf86vm", "pthread", "/opt/etsf/lib/libxc.a", "nfd", "gtk-3", "z", "gdk-3", "glib-2.0" }


	filter { "system:linux", "configurations:Debug" }
		libdirs { "../../C++/GraphicsEngine/bin/linux/Debug-x64/GraphicsEngine" }


	filter { "system:linux", "configurations:Release" }
		libdirs { "../../C++/GraphicsEngine/bin/linux/Release-x64/GraphicsEngine" }

	filter { "system:windows" }
		defines{"GLEW_STATIC", "GENG_STATIC"}
		pchheader("pch.h")
		pchsource("%{wks.location}/src/pch.cpp")
		files {
		  "%{wks.location}/Eigen.natvis"
		}
		includedirs
		{
		  "%{wks.location}/src/",
		  "../../Libraries/json-develop/include",
		  "../../Libraries/eigen-3.4.90",
		  "../../C++/GraphicsEngine/src",
		  "../../Libraries/glm",
		  "../../Libraries/ImGui",
		  "../../Libraries/glew-2.1.0/include",
		  "../../Libraries/glfw-3.3.2/include",
		  "../../Libraries/",
		  "../../Libraries/nativefiledialog/src/include",
		  "../../Libraries/boost_1_83_0",
		  "../../Libraries/libxc-6.2.2/src",
		  "../../Libraries/cppitertools"
		}
		links
		{
          "GraphicsEngine",
		  "opengl32",
		  "ImGui",
		  "xc.lib"
		}
		linkoptions { '/NODEFAULTLIB:"LIBCMTD" /NODEFAULTLIB:"LIBCMT"' }

	filter { "configurations:Debug", "system:Windows" }
		defines { "CONFIG_DEBUG", "_DEBUG" }
		symbols "On"
		libdirs { 
		    "../../Libraries/glew-2.1.0/lib/Debug/x64/",
		    "../../Libraries/glfw-3.3.2/lib/x64/",
		    "../../Libraries/boost_1_83_0/stage/lib/",
		    "../../Libraries/ImGui/bin/windows/Debug-x64",
		    "../../Libraries/libxc-6.2.2/VS/Release",
		    "../../C++/GraphicsEngine/bin/windows/Debug-x64/GraphicsEngine",
			"../../Libraries/nativefiledialog/build/lib/Release/x64"
		     }
		links { "glew32sd", "glfw3_debug", "libboost_filesystem-vc143-mt-gd-x64-1_83.lib", "nfd.lib"}
		staticruntime "off"
		buildoptions "/MDd"

	filter { "configurations:Release", "system:Windows" }
		defines "CONFIG_RELEASE"
		optimize "On"
		libdirs {   "../../Libraries/glew-2.1.0/lib/Release/x64/",
		    "../../Libraries/glfw-3.3.2/lib/x64/",
		    "../../Libraries/boost_1_83_0/stage/lib/",
		    "../../Libraries/ImGui/bin/windows/Release-x64",
		    "../../Libraries/libxc-6.2.2/VS/Release",
		     "../../C++/GraphicsEngine/bin/windows/Release-x64/GraphicsEngine",
			 "../../Libraries/nativefiledialog/build/lib/Release/x64"
		     }
		links { "glew32s", "glfw3", "libboost_filesystem-vc143-mt-x64-1_83.lib", "nfd.lib"}
		staticruntime "off"
		buildoptions "/MD"

	filter { "configurations:Debug", "system:Linux" }
		defines { "CONFIG_DEBUG", "_DEBUG" }
		symbols "On"

	filter { "configurations:Release", "system:Linux" }
		defines { "CONFIG_RELEASE" }
		optimize "Speed"
		vectorextensions "AVX2"
		symbols "Off"
