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
	
	filter { "system:windows" }
		defines{"GLEW_STATIC", "_CONSOLE",}
		pchheader("pch.h")
		pchsource("%{wks.location}/src/pch.cpp")
		postbuildcommands {
		  ("call xcopy.exe \"%{prj.location}Assets\" \"%{wks.location}" .. outputdir .. "Assets\" /s /y /i /d")
		}
		files {
		  "%{wks.location}/Eigen.natvis"
		}
        publicLibDir = "/media/charlie/Data/dev/Libraries/"
		includedirs
		{
		  "%{wks.location}/src/",
		  "D:/dev/Libraries/json-develop/include",
		  "D:/dev/Libraries/eigen-3.4.0",
		  "D:/dev/Libraries/glm",
		  "D:/dev/Libraries/ImGui",
		  "D:/dev/Libraries/glew-2.1.0/include",
		  "D:/dev/Libraries/glfw-3.3.2/include",
		  "D:/dev/Libraries/",
		  "D:/dev/Libraries/nativefiledialog/src/include",
		  "D:/dev/Libraries/boost_1_83_0",
		  "D:/dev/Libraries/libxc-6.2.2/src"
		}
		links
		{
		  "opengl32",
		  "ImGui"
		}
		linkoptions { '/NODEFAULTLIB:"LIBCMTD" /NODEFAULTLIB:"LIBCMT"' }
		
	filter {"system:linux"}
		includedirs
		{
		  "%{wks.location}/src/",
		  "/usr/include",
		  "/usr/local/include",
		  "/opt/etsf/include",
          "/usr/local/include/eigen3",
		  "/media/charlie/Data/dev/Libraries/json-develop/include",
		  "/media/charlie/Data/dev/Libraries/eigen-3.4.0",
		  "/media/charlie/Data/dev/Libraries/glm",
          "/media/charlie/Data/dev/Libraries/",
		  "/media/charlie/Data/dev/Libraries/ImGui",
          "/media/charlie/Data/dev/Libraries/nativefiledialog/src/include"
		}
        libdirs
        {
            "/media/charlie/Data/dev/Libraries/ImGui/bin/linux/Release-x64",
            "/media/charlie/Data/dev/Libraries/nativefiledialog/build/lib/Release/x64"
        }
		links { "glfw3", "ImGui", "boost_filesystem", "GLEW.a", "GL", "m", "Xrandr", "Xi", "X11", "Xxf86vm", "pthread", "/opt/etsf/lib/libxc.a", "nfd", "gtk-3", "z", "gdk-3", "glib-2.0" }
		--linkoptions { "`pkg-config --libs gtk+-3.0`" }


	filter { "configurations:Debug", "system:Windows" }
		defines { "CONFIG_DEBUG", "_DEBUG" }
		symbols "On"
		libdirs { 
		    "D:/dev/Libraries/glew-2.1.0/lib/Debug/x64/",
		    "D:/dev/Libraries/glfw-3.3.2/lib/x64/",
		    "D:/dev/Libraries/boost_1_83_0/stage/lib/",
		    "D:/dev/Libraries/ImGui/bin/windows/Debug-x64",
		    "D:/dev/Libraries/libxc-6.2.2/VS/Release" }
		links { "glew32sd", "glfw3_debug", "libboost_filesystem-vc143-mt-gd-x64-1_83.lib", "ImGui", "xc.lib",   "nativefiledialog/build/lib/Release/x64/nfd.lib" }
		staticruntime "off"
		buildoptions "/MDd"

	filter { "configurations:Release", "system:Windows" }
		defines "CONFIG_RELEASE"
		optimize "On"
		libdirs {   "D:/dev/Libraries/glew-2.1.0/lib/Release/x64/",
		    "D:/dev/Libraries/glfw-3.3.2/lib/x64/",
		    "D:/dev/Libraries/boost_1_83_0/stage/lib/",
		    "D:/dev/Libraries/ImGui/bin/windows/Release-x64",
		    "D:/dev/Libraries/libxc-6.2.2/VS/Release" }
		links { "glew32s", "glfw3", "libboost_filesystem-vc143-mt-x64-1_83.lib", "ImGui", "xc.lib",   "nativefiledialog/build/lib/Release/x64/nfd.lib" }
		staticruntime "off"
		buildoptions "/MD"

	filter { "configurations:Debug", "system:Linux" }
		defines { "CONFIG_DEBUG", "_DEBUG" }
		symbols "On"

	filter { "configurations:Release", "system:Linux" }
		defines { "CONFIG_RELEASE" }
		optimize "On"
		symbols "Off"
