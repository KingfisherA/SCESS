workspace "ImPlot"
	architecture "x64"

	configurations
	{
		"Debug",
		"Release"
	}
	
	outputdir = "%{os.host()}/%{cfg.buildcfg}-x64/"

project "ImPlot"
	kind "StaticLib"
	language "C++"
	cppdialect "C++20"
	systemversion "latest"
	
	targetdir("../ImPlot/bin/" .. outputdir)
	objdir("../ImPlot/bin-obj/" .. outputdir)
	debugdir("../ImPlot/bin/" .. outputdir)

	files
	{
	  "../ImPlot/**.cpp",
	  "../ImPlot/**.hpp",
	  "./**.h"
	}
	
	includedirs
	{
	  "../ImPlot/",
	  "../GraphicsEngine-1.0.0/include/ImGui"
	}
	
	filter {  "configurations:Debug" }
		defines { "CONFIG_DEBUG", "_DEBUG" }
		optimize "Off"
		symbols "On"
		
	filter { "configurations:Release" }
		defines { "CONFIG_RELEASE" }
		optimize "Speed"
		symbols "Off"
	
