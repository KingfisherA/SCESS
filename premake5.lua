workspace "SCESS"
    architecture "x64"

    configurations
    {
    "Debug",
    "Release"
    }

outputdir = "%{prj.name}/%{cfg.buildcfg}-x64/"
publicLibDir = "C:/dev/Libraries/"
libDir = "%{wks.location}/libs/"

project "SCESS"

     kind "ConsoleApp"
     language "C++"

     targetdir ("bin/" .. outputdir)
     objdir ("bin-obj/" .. outputdir)
     debugdir ("bin/" .. outputdir)

     pchheader ("pch.h")
     pchsource ("%{wks.location}pch.cpp")

     --postbuildcommands  {
       -- ("call xcopy.exe \"%{prj.location}shaders\\\" \"%{wks.location}bin\\%{prj.name}\\Debug-x64\\assets\\\" /s /y /i /d")
     --}

     files
     {
          "%{wks.location}/src/**.h",
          "%{wks.location}/src/**.hpp",
          "%{wks.location}/src/**.cpp"
     }

     includedirs
     {
        "%{wks.location}/src/",
        publicLibDir .. "glew-2.1.0/include/",
        publicLibDir .. "glfw-3.3.2/include/",
        publicLibDir .. "glm/",
		libDir .. "json-develop/include",
		libDir .. "eigen-3.4.0",
		publicLibDir .. "boost_1_83_0/",
        publicLibDir .. "ImGui/"
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
         publicLibDir .. "glfw-3.3.2/lib/x64/"}
      links {"glew32sd", "glfw3_debug"}
      staticruntime "off"
      buildoptions "/MDd"

     filter "configurations:Release"
        defines "CONFIG_RELEASE"
        optimize "On"
        libdirs {publicLibDir .. "glew-2.1.0/lib/Release/x64/",
        publicLibDir .. "glfw-3.3.2/lib/x64/"}
        links {"glew32s", "glfw3"}
        staticruntime "off"
        buildoptions "/MD"
--
--
include (publicLibDir .. "ImGui/")