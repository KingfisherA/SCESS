workspace "SCESS"
	architecture "x64"

	configurations
	{
		"Debug",
		"Release"
	}

--	system
--	{
--		"Windows",
--		"Linux"
--	}

publicLibDir = ""
publicLibDirC = ""

filter{"system:Windows"}
	publicLibDir = "C:/dev/Libraries/C++/"
	publicLibDirC = "C:/dev/Libraries/C/"

filter{"system:Linux"}
	publicLibDir = "/media/charlie/Data/dev/Libraries/C++/"
	publicLibDirC = "/media/charlie/Data/dev/Libraries/C/"

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

	defines { "UNI_OPENGL", "_CONSOLE", "GLEW_STATIC" }

	filter { "system:Windows" }
		pchheader("pch.h")
		pchsource("%{wks.location}/src/pch.cpp")
		postbuildcommands {
		  ("call xcopy.exe \"%{prj.location}Assets\" \"%{wks.location}" .. outputdir .. "Assets\" /s /y /i /d")
		}
		files {
		  "%{wks.location}/Eigen.natvis"
		}



		filter "system:windows"
		defines
		{
		  "PLATFORM_WIN"
		}

		links
		{
		  "opengl32",
		  "ImGui"
		}

		linkoptions { '/NODEFAULTLIB:"LIBCMTD" /NODEFAULTLIB:"LIBCMT"' }
		
	filter {"system:Linux"}
		includedirs
		{
		  "%{wks.location}/src/",
		  "/usr/include",
		  "/usr/local/include",
		  publicLibDir .. "json-develop/include",
		  publicLibDir .. "eigen-3.4.0",
		  publicLibDir .. "glm",
		  publicLibDir .. "ImGui"
		}
		buildoptions{"$(pkg-config --cflags --libs gtk+-3.0)"}


	filter { "configurations:Debug", "system:Windows" }
		defines { "CONFIG_DEBUG", "_DEBUG" }
		symbols "On"
		libdirs { publicLibDir .. "glew-2.1.0/lib/Debug/x64/",
		  publicLibDir .. "glfw-3.3.2/lib/x64/",
		  publicLibDir .. "boost_1_83_0/stage/lib/",
		  publicLibDir .. "ImGui/bin/ImGui/Debug-x64/ImGui",
		  publicLibDirC .. "libxc-6.2.2/VS/Release" }
		links { "glew32sd", "glfw3_debug", "libboost_filesystem-vc143-mt-gd-x64-1_83.lib", "ImGui", "xc.lib", publicLibDirC .. "nativefiledialog/build/lib/Release/x64/nfd.lib" }
		staticruntime "off"
		buildoptions "/MDd"

	filter { "configurations:Release", "system:Windows" }
		defines "CONFIG_RELEASE"
		optimize "On"
		libdirs { publicLibDir .. "glew-2.1.0/lib/Release/x64/",
		  publicLibDir .. "glfw-3.3.2/lib/x64/",
		  publicLibDir .. "boost_1_83_0/stage/lib/",
		  publicLibDir .. "ImGui/bin/ImGui/Release-x64/ImGui",
		  publicLibDirC .. "libxc-6.2.2/VS/Release" }
		links { "glew32s", "glfw3", "libboost_filesystem-vc143-mt-x64-1_83.lib", "ImGui", "xc.lib", publicLibDirC .. "nativefiledialog/build/lib/Release/x64/nfd.lib" }
		staticruntime "off"
		buildoptions "/MD"

	filter { "configurations:Debug", "system:Linux" }
		defines { "CONFIG_DEBUG", "_DEBUG" }
		symbols "On"
		links { "libglfw3.a", "/usr/local/lib/libImGui.a", "/usr/local/lib/libboost_filesystem.a", "/usr/lib64/libGLEW.a", "GL", "m", "Xrandr", "Xi", "X11","Xxf86vm", "pthread", "/usr/local/lib/nfd.a", "/usr/local/lib/libnfd.a", "/opt/etsf/lib/libxc.a" }
		staticruntime "off"

	filter { "configurations:Release", "system:Linux" }
		defines { "CONFIG_RELEASE" }
		optimize "On"
		links { "libglfw3.a", "/usr/local/lib/libImGui.a", "/usr/local/lib/libboost_filesystem.a", "/usr/lib64/libGLEW.a","/usr/lib64/libGLEW.so", "GL", "m", "Xrandr", "Xi", "X11","Xxf86vm", "pthread", "/usr/local/lib/nfd.a","/usr/local/lib/libnfd.a", "/opt/etsf/lib/libxc.a" }
		staticruntime "off"
