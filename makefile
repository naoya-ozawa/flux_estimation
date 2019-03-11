all:	flux.cpp
	`root-config --cxx --cflags` -o flux flux.cpp `root-config --glibs`

mbflux:	mbflux.cpp
	`root-config --cxx --cflags` -o mbflux mbflux.cpp `root-config --glibs`
