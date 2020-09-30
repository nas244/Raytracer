workspace "RayTracer"
	architecture "x64"

	configurations
	{
		"Debug",
		"Release"
	}
outputdir = "%{cfg.buildcfg}-%{cfg.system}-%{cfg.architecture}"

IncludeDir = {}
IncludeDir["glm"] = "RayTracer/lib/glm"
IncludeDir["CImg"] = "RayTracer/lib/CImg"
IncludeDir["tinyObj"] = "RayTracer/lib/tinyObj"

project "RayTracer"
	location "TestProgram"
	kind "ConsoleApp"
	language "C++"
	cppdialect "C++17"

	targetdir ("bin/" .. outputdir .. "/%{prj.name}")
	objdir ("bin-int/" .. outputdir .. "/%{prj.name}")

	files
	{
		"%{prj.name}/src/**.h",
		"%{prj.name}/src/**.cpp"
	}

	includedirs
	{
		"%{IncludeDir.glm}",
		"%{IncludeDir.CImg}",
		"%{IncludeDir.tinyObj}"
	}

