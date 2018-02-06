WORKDIR = `pwd`

CC = gcc
CXX = g++
AR = ar
LD = g++
WINDRES = windres

INC = 
CFLAGS = -std=gnu++11 -Wall 
RESINC = 
LIBDIR = 
LIB = 
LDFLAGS = 

INC_DEBUG = $(INC) -I/home/jco35/gsl/include -I/home/jco35/Downloads/Healpix_3.31/include -I/home/jco35/cwru-pat-cosmic_topology/T3-correlations 
CFLAGS_DEBUG = $(CFLAGS) -g
RESINC_DEBUG = $(RESINC)
RCFLAGS_DEBUG = $(RCFLAGS)
LIBDIR_DEBUG = $(LIBDIR) -L/home/jco35/Downloads/cfitsio -L/home/jco35/Downloads/Healpix_3.31/lib -L/home/jco35/gsl/lib 
LIB_DEBUG = $(LIB)/home/jco35/Downloads/Healpix_3.31/lib/libchealpix.a /home/jco35/Downloads/cfitsio/libcfitsio.a /home/jco35/gsl/lib/libgsl.a
LDFLAGS_DEBUG = $(LDFLAGS)
OBJDIR_DEBUG = obj/Debug
DEP_DEBUG = 
OUT_DEBUG = bin/Debug/T3_correlations

INC_RELEASE = $(INC) -I/home/jco35/gsl/include -I/home/jco35/Downloads/Healpix_3.31/include -I/home/jco35/cwru-pat-cosmic_topology/T3-correlations 
CFLAGS_RELEASE = $(CFLAGS) -g
RESINC_RELEASE = $(RESINC)
RCFLAGS_RELEASE = $(RCFLAGS)
LIBDIR_RELEASE = $(LIBDIR) -L/home/jco35/Downloads/cfitsio -L/home/jco35/Downloads/Healpix_3.31/lib -L/home/jco35/gsl/lib 
LIB_RELEASE = $(LIB)/home/jco35/Downloads/Healpix_3.31/lib/libchealpix.a /home/jco35/Downloads/cfitsio/libcfitsio.a /home/jco35/gsl/lib/libgsl.a
LDFLAGS_RELEASE = $(LDFLAGS) -s
OBJDIR_RELEASE = obj/Release
DEP_RELEASE = 
OUT_RELEASE = bin/Release/T3_correlations

OBJ_DEBUG = $(OBJDIR_DEBUG)/HPCC_Spherical.o $(OBJDIR_DEBUG)/correlation_matrices.o $(OBJDIR_DEBUG)/eigenbasis_lin_T3.o $(OBJDIR_DEBUG)/main.o $(OBJDIR_DEBUG)/readToVec.o $(OBJDIR_DEBUG)/rotAngles.o

OBJ_RELEASE = $(OBJDIR_RELEASE)/HPCC_Spherical.o $(OBJDIR_RELEASE)/correlation_matrices.o $(OBJDIR_RELEASE)/eigenbasis_lin_T3.o $(OBJDIR_RELEASE)/main.o $(OBJDIR_RELEASE)/readToVec.o $(OBJDIR_RELEASE)/rotAngles.o

all: debug release

clean: clean_debug clean_release

before_debug: 
	test -d bin/Debug || mkdir -p bin/Debug
	test -d $(OBJDIR_DEBUG) || mkdir -p $(OBJDIR_DEBUG)

after_debug: 

debug: before_debug out_debug after_debug

out_debug: before_debug $(OBJ_DEBUG) $(DEP_DEBUG)
	$(LD) $(LIBDIR_DEBUG) -o $(OUT_DEBUG) $(OBJ_DEBUG)  $(LDFLAGS_DEBUG) $(LIB_DEBUG)

$(OBJDIR_DEBUG)/HPCC_Spherical.o: HPCC_Spherical.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c HPCC_Spherical.cpp -o $(OBJDIR_DEBUG)/HPCC_Spherical.o

$(OBJDIR_DEBUG)/correlation_matrices.o: correlation_matrices.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c correlation_matrices.cpp -o $(OBJDIR_DEBUG)/correlation_matrices.o

$(OBJDIR_DEBUG)/eigenbasis_lin_T3.o: eigenbasis_lin_T3.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c eigenbasis_lin_T3.cpp -o $(OBJDIR_DEBUG)/eigenbasis_lin_T3.o

$(OBJDIR_DEBUG)/main.o: main.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c main.cpp -o $(OBJDIR_DEBUG)/main.o

$(OBJDIR_DEBUG)/readToVec.o: readToVec.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c readToVec.cpp -o $(OBJDIR_DEBUG)/readToVec.o

$(OBJDIR_DEBUG)/rotAngles.o: rotAngles.cpp
	$(CXX) $(CFLAGS_DEBUG) $(INC_DEBUG) -c rotAngles.cpp -o $(OBJDIR_DEBUG)/rotAngles.o

$(OBJDIR_DEBUG)/T3-correlations: $(OBJ_DEBUG)
	$(CXX) $(CFLAGS_DEBUG) -o T3-correlations $(OBJ_DEBUG) $(LIB_DEBUG)

clean_debug: 
	rm -f $(OBJ_DEBUG) $(OUT_DEBUG)
	rm -rf bin/Debug
	rm -rf $(OBJDIR_DEBUG)

before_release: 
	test -d bin/Release || mkdir -p bin/Release
	test -d $(OBJDIR_RELEASE) || mkdir -p $(OBJDIR_RELEASE)

after_release: 

release: before_release out_release after_release

out_release: before_release $(OBJ_RELEASE) $(DEP_RELEASE)
	$(LD) $(LIBDIR_RELEASE) -o $(OUT_RELEASE) $(OBJ_RELEASE) $(LDFLAGS_RELEASE) $(LIB_RELEASE)

$(OBJDIR_RELEASE)/HPCC_Spherical.o: HPCC_Spherical.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c HPCC_Spherical.cpp -o $(OBJDIR_RELEASE)/HPCC_Spherical.o

$(OBJDIR_RELEASE)/correlation_matrices.o: correlation_matrices.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c correlation_matrices.cpp -o $(OBJDIR_RELEASE)/correlation_matrices.o

$(OBJDIR_RELEASE)/eigenbasis_lin_T3.o: eigenbasis_lin_T3.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c eigenbasis_lin_T3.cpp -o $(OBJDIR_RELEASE)/eigenbasis_lin_T3.o

$(OBJDIR_RELEASE)/main.o: main.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c main.cpp -o $(OBJDIR_RELEASE)/main.o

$(OBJDIR_RELEASE)/readToVec.o: readToVec.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c readToVec.cpp -o $(OBJDIR_RELEASE)/readToVec.o

$(OBJDIR_RELEASE)/rotAngles.o: rotAngles.cpp
	$(CXX) $(CFLAGS_RELEASE) $(INC_RELEASE) -c rotAngles.cpp -o $(OBJDIR_RELEASE)/rotAngles.o

$(OBJDIR_RELEASE)/T3-correlations: $(OBJ_RELEASE)
	$(CXX) $(CFLAGS_RELEASE) -o T3-correlations $(OBJ_RELEASE) $(LIB_RELEASE)

clean_release: 
	rm -f $(OBJ_RELEASE) $(OUT_RELEASE)
	rm -rf bin/Release
	rm -rf $(OBJDIR_RELEASE)

.PHONY: before_debug after_debug clean_debug before_release after_release clean_release
