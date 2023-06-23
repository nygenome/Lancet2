# https://en.wikichip.org/wiki/intel/cpuid, https://en.wikichip.org/wiki/amd/cpuid
# Currently only the x86 CPUs available on Google Cloud Compute are auto-detected
# https://cloud.google.com/compute/docs/cpu-platforms#x86_processors
# https://gcc.gnu.org/onlinedocs/gcc/x86-Options.html
# For all other CPUs, -march=native -mtune=native is used to let compiler decide

function(DetectHostCPUInfo)

	if (NOT EXISTS "/proc/cpuinfo")
		message(FATAL_ERROR "Could not find /proc/cpuinfo")
	endif ()

	file(READ "/proc/cpuinfo" _cpuinfo)
	string(REGEX REPLACE ".*vendor_id[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" CPU_VENDOR_ID "${_cpuinfo}")
	string(REGEX REPLACE ".*cpu family[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" CPU_FAMILY "${_cpuinfo}")
	string(REGEX REPLACE ".*model[ \t]*:[ \t]+([a-zA-Z0-9_-]+).*" "\\1" CPU_MODEL "${_cpuinfo}")

	set(CPU_VENDOR_ID ${CPU_VENDOR_ID} PARENT_SCOPE)
	set(CPU_FAMILY ${CPU_FAMILY} PARENT_SCOPE)
	set(CPU_MODEL ${CPU_MODEL} PARENT_SCOPE)

endfunction(DetectHostCPUInfo)

function(DetectHostCpuArch)

	DetectHostCPUInfo()

	# Default fallback arch to use
	set(DETECTED_CPU_ARCH "native")
	set(LANCET_CPU_NAME "Unknown")

	if (CPU_VENDOR_ID STREQUAL "GenuineIntel" AND CPU_FAMILY EQUAL 6)
		if (CPU_MODEL EQUAL 143)
			# C3 instances & Sapphire Rapids
			set(DETECTED_CPU_ARCH "sapphirerapids")
			set(LANCET_CPU_NAME "Intel SapphireRapids")
		elseif (CPU_MODEL EQUAL 106 OR CPU_MODEL EQUAL 108)
			# N2, M3 instances & Ice Lake
			set(DETECTED_CPU_ARCH "icelake-server")
			set(LANCET_CPU_NAME "Intel IceLake")
		elseif (CPU_MODEL EQUAL 85)
			# N1, N2, C2, M2 instances & Skylake (or) Cascade Lake
			set(DETECTED_CPU_ARCH "skylake-avx512")
			set(LANCET_CPU_NAME "Intel Skylake/CascadeLake")
		elseif (CPU_MODEL EQUAL 79 OR CPU_MODEL EQUAL 86)
			# N1 instances & Broadwell
			set(DETECTED_CPU_ARCH "broadwell")
			set(LANCET_CPU_NAME "Intel Broadwell")
		elseif (CPU_MODEL EQUAL 63)
			# N1 instances & Haswell
			set(DETECTED_CPU_ARCH "haswell")
			set(LANCET_CPU_NAME "Intel Haswell")
		elseif (CPU_MODEL EQUAL 62)
			# N1 instances & Ivy Bridge
			set(DETECTED_CPU_ARCH "ivybridge")
			set(LANCET_CPU_NAME "Intel IvyBridge")
		elseif (CPU_MODEL EQUAL 45)
			# N1 instances & Sandy Bridge
			set(DETECTED_CPU_ARCH "sandybridge")
			set(LANCET_CPU_NAME "Intel SandyBridge")
		endif ()
	endif ()

	if (CPU_VENDOR_ID STREQUAL "AuthenticAMD")
		if (CPU_FAMILY EQUAL 25 AND CPU_MODEL EQUAL 1)
			# N2D, T2D, C2D instances & Milan
			set(DETECTED_CPU_ARCH "znver3")
			set(LANCET_CPU_NAME "AMD Milan")
		elseif (CPU_FAMILY EQUAL 23 AND CPU_MODEL EQUAL 49)
			# N2D instances & Rome
			set(DETECTED_CPU_ARCH "znver2")
			set(LANCET_CPU_NAME "AMD Rome")
		endif ()
	endif ()

	set(LANCET_BUILD_ARCH ${DETECTED_CPU_ARCH} CACHE STRING "Host CPU Architecture" FORCE)
	message(STATUS "Detected CPU information: ${LANCET_CPU_NAME} Family ${CPU_FAMILY} Model ${CPU_MODEL}")

	if (LANCET_CPU_NAME STREQUAL "Unknown")
		message(STATUS "Could not detect CPU, letting compiler decide best native optimizations for host machine")
	endif ()

endfunction(DetectHostCpuArch)
