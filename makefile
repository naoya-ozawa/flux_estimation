all:	flux.cpp
	`root-config --cxx --cflags` -o flux flux.cpp `root-config --glibs`
