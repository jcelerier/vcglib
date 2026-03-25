cmake_minimum_required(VERSION 3.13)
project(trimesh_fast_winding_number)

if (VCG_HEADER_ONLY)
	set(SOURCES
		trimesh_fast_winding_number.cpp
		${VCG_INCLUDE_DIRS}/wrap/ply/plylib.cpp)
endif()

add_executable(trimesh_fast_winding_number
	${SOURCES})

target_link_libraries(
	trimesh_fast_winding_number
	PUBLIC
		vcglib
	)

