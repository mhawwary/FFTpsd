# Specifying some environmental variables for Linux, note this can be done in the shell prompt
COMP	= GCC
CODE	= RELEASE
OPENMP	= NO

# Specifing Standard Variables:
CXX	= g++ -std=gnu++11 #-pedantic-errors # c++ gcc compiler
CXXFLAGS=       # C++ compiler flags
LDLFLAGS=	# linker flags
CPPFLAGS=	# c/c++ preprocessor flags

OPTS	= 	# optimization flags and other options

# Includes
OPTS	+= -I include

ifeq ($(CODE),RELEASE)
	ifeq ($(COMP),GCC)
		OPTS	+= -O3 
	endif
	
	ifeq ($(COMP),INTEL)
		OPTS	+= -xHOST -fast
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
OBJ	= obj/
BIN	= bin/
INC	= include/      

vpath %.cpp src
vpath %.c src
vpath %.o   obj
vpath %.h include src
vpath %.hpp include src

#source files:
src_filename = FFTpsd
OBJo = $(OBJ)$(src_filename).o 
src_cpp = $(src_filename).cpp

# Objects
OBJS	=  $(OBJo) $(OBJ)fft_driver.o $(OBJ)fft.o $(OBJ)fft_param.o $(OBJ)general_tools.o $(OBJ)string_to_type.o  # objects 
INCLS	= 

# Compile
.PHONY: default help clean

default: all
help:	
	@echo 'help'

all: fftpsd
fftpsd: $(OBJS)
	$(CXX) $(OPTS) -o $(BIN)$@ $+

$(OBJ)%.o : %.cpp 
	$(CXX) $(OPTS) -c -o $@ $<
$(OBJ)%.o : %.c 
	$(CXX) $(OPTS) -c -o $@ $<

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

