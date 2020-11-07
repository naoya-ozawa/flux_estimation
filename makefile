all:	flux.cpp
	`root-config --cxx --cflags` -o flux flux.cpp `root-config --glibs`

mbflux:	mbflux.cpp
	`root-config --cxx --cflags` -o mbflux mbflux.cpp `root-config --glibs`

rdep:	distance_dependence.cpp
	`root-config --cxx --cflags` -o rdep distance_dependence.cpp `root-config --glibs`

material:	material_thickness.cpp
	`root-config --cxx --cflags` -o material material_thickness.cpp `root-config --glibs`
