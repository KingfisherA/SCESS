workspace "SCESS"
    architecture "x64"

    configurations
    {
    "Debug",
    "Release"
    }

outputdir = "%{prj.name}/%{cfg.buildcfg}-x64/"
publicLibDir = "C:/dev/Libraries/C++/"
publicLibDirC = "C:/dev/Libraries/C/"
libDir = "%{wks.location}/libs/"

project "SCESS"

     kind "ConsoleApp"
     language "C++"

     targetdir ("bin/" .. outputdir)
     objdir ("bin-obj/" .. outputdir)
     debugdir ("bin/" .. outputdir)

	 pchheader("pch.h")
	 pchsource("%{wks.location}/src/pch.cpp")

    postbuildcommands {
      ("call xcopy.exe \"%{prj.location}Assets\" \"%{wks.location}" .. outputdir .. "Assets\" /s /y /i /d")
    }

     files
     {
          "%{wks.location}/src/**.h",
          "%{wks.location}/src/**.hpp",
          "%{wks.location}/src/**.cpp",
          "%{wks.location}/resource.h",
          "%{wks.location}/SCESS.APS",
          "%{wks.location}/SCESS.rc",
          "%{wks.location}/Eigen.natvis"
     }

     includedirs
     {
        "%{wks.location}/src/",
        publicLibDir .. "glew-2.1.0/include/",
        publicLibDir .. "glfw-3.3.2/include/",
        publicLibDir .. "glm/",
		publicLibDir .. "json-develop/include",
		publicLibDir .. "eigen-3.4.0",
		publicLibDir .. "boost_1_83_0/",
        publicLibDir .. "ImGui/",
		publicLibDirC .. "libxc-6.2.2/src",
		publicLibDirC .. "nativefiledialog/src/include"
     }

     defines {"UNI_OPENGL", "_CONSOLE", "GLEW_STATIC"}

     filter "system:windows"
        cppdialect "C++20"
        staticruntime "On"
        systemversion "latest"

        defines 
        {
          "PLATFORM_WIN"
        }

         links
         {
           "opengl32",
           "ImGui"
         }

        linkoptions {'/NODEFAULTLIB:"LIBCMTD" /NODEFAULTLIB:"LIBCMT"'}

    filter "configurations:Debug"
        defines {"CONFIG_DEBUG", "_DEBUG"}
        symbols "On"
        libdirs {publicLibDir .. "glew-2.1.0/lib/Debug/x64/",
         publicLibDir .. "glfw-3.3.2/lib/x64/",
		 publicLibDir .. "boost_1_83_0/stage/lib/",
		 publicLibDir .. "ImGui/bin/ImGui/Debug-x64/ImGui",
		 publicLibDirC .. "libxc-6.2.2/VS/Release"}
      links {"glew32sd", "glfw3_debug", "libboost_filesystem-vc143-mt-gd-x64-1_83.lib", "ImGui", "xc.lib", publicLibDirC .. "nativefiledialog/build/lib/Release/x64/nfd.lib"}
      staticruntime "off"
      buildoptions "/MDd"

     filter "configurations:Release"
        defines "CONFIG_RELEASE"
        optimize "On"
        libdirs {publicLibDir .. "glew-2.1.0/lib/Release/x64/",
        publicLibDir .. "glfw-3.3.2/lib/x64/",
		publicLibDir .. "boost_1_83_0/stage/lib/",
		publicLibDir .. "ImGui/bin/ImGui/Release-x64/ImGui",
		publicLibDirC .. "libxc-6.2.2/VS/Release"}
        links {"glew32s", "glfw3", "libboost_filesystem-vc143-mt-x64-1_83.lib", "ImGui", "xc.lib", publicLibDirC .. "nativefiledialog/build/lib/Release/x64/nfd.lib"}
        staticruntime "off"
        buildoptions "/MD"
	