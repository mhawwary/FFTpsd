# Specifying some environmental variables for Linux, note this can be done in the shell prompt
COMP	= GCC
CODE	= RELEASE
OPENMP	= NO
PLATFORM = LINUX

# Specifing Standard Variables:
ifeq ($(PLATFORM), LINUX)
	CXX = g++ #-pedantic-errors # c++ gcc compiler
	LDLFLAGS =  # linker flags
        BINNAME = fftpsd
else
	CXX = x86_64-w64-mingw32-g++ # c++ gcc compiler
	LDLFLAGS = -static-libgcc -static-libstdc++ # linker flags
        BINNAME = fftpsd.exe
endif
CXXFLAGS= -std=gnu++11 # C++ compiler flags
CPPFLAGS=  # c/c++ preprocessor flags

OPTS	= # optimization flags and other options

# Includes
OPTS	+= -I include

ifeq ($(CODE),RELEASE)
    ifeq ($(COMP),GCC)
        OPTS += -O3  
	endif
	
	ifeq ($(COMP),INTEL)
	    OPTS   += -xHOST -fast
	endif	
endif

ifeq ($(OPENMP),YES)	
	OPTS	+=  -lgomp -fopenmp 
endif

ifeq ($(CODE),DEBUG)
	ifeq ($(COMP),GCC)
		OPTS	+= -fbuiltin -g -Wall #-Werror
	endif
	
	ifeq ($(COMP),INTEL)
		OPTS	+= -xHOST -fast
	endif	
endif


# Source

SRC	= src/
OBJ	= ./obj/
BIN	= ./bin/
INC	= include/      
OUT_DIR = $(OBJ) $(BIN)

vpath %.cpp $(SRC)
vpath %.c $(SRC)
vpath %.o $(OBJ)
vpath %.h $(INC) $(SRC)
vpath %.hpp $(INC) $(SRC)

#source files:
src_filename = FFTpsd
OBJo = $(OBJ)$(src_filename).o 
src_cpp = $(src_filename).cpp

# Objects
OBJS	=  $(OBJo) $(OBJ)fft_driver.o $(OBJ)fft.o $(OBJ)fft_param.o $(OBJ)general_tools.o $(OBJ)string_to_type.o  # objects 
INCLS	= 

# Compile
.PHONY: default help clean directories

default: all
help:	
	@echo 'help'

all: directories $(BINNAME)

directories: ${OUT_DIR}

${OUT_DIR}:
	mkdir -p ${OUT_DIR}

$(BINNAME): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OPTS) $(LDLFLAGS) -o $(BIN)$@ $+

$(OBJ)%.o : %.cpp
	$(CXX) $(CXXFLAGS) $(OPTS) -c -o $@ $<
$(OBJ)%.o : %.c 
	$(CXX) $(CXXFLAGS) $(OPTS) -c -o $@ $<

$(OBJo): $(src_cpp)
$(OBJ)fft_driver.o: fft_driver.cpp 
$(OBJ)fft.o: fft.cpp
$(OBJ)fft_param.o: fft_param.cpp
$(OBJ)general_tools.o: general_tools.cpp
$(OBJ)string_to_type.o: string_to_type.c

clean:
	rm -f ./$(OBJ)*.o ./$(BIN)*.exe 
	rm -f *.exe *.o ./$(BIN)fftpsd* fftpsd*
	@echo  removing all object and executable files
