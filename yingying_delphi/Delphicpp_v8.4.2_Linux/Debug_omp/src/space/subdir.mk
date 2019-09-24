################################################################################
# Automatically-generated file. Do not edit!
################################################################################

#Added an external fft page
# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/space/space_mpi.cpp \
../src/space/grid_space.cpp \
../src/space/space_crgarr.cpp \
../src/space/space_cube.cpp \
../src/space/space_eps.cpp \
../src/space/space_indver.cpp \
../src/space/space_msrf.cpp \
../src/space/space_run.cpp \
../src/space/space_sas.cpp \
../src/space/space_sclbp.cpp \
../src/space/space_setout.cpp \
../src/space/space_validateInput.cpp \
../src/space/space_vwtms.cpp \
../src/space/space_vwtms_inc.cpp \
../src/space/space_setgaussian.cpp

OBJS += \
./src/space/space_mpi.o \
./src/space/grid_space.o \
./src/space/space_crgarr.o \
./src/space/space_cube.o \
./src/space/space_eps.o \
./src/space/space_indver.o \
./src/space/space_msrf.o \
./src/space/space_run.o \
./src/space/space_sas.o \
./src/space/space_sclbp.o \
./src/space/space_setout.o \
./src/space/space_validateInput.o \
./src/space/space_vwtms.o \
./src/space/space_vwtms_inc.o \
./src/space/space_setgaussian.o

CPP_DEPS += \
./src/space/space_mpi.d \
./src/space/grid_space.d \
./src/space/space_crgarr.d \
./src/space/space_cube.d \
./src/space/space_eps.d \
./src/space/space_indver.d \
./src/space/space_msrf.d \
./src/space/space_run.d \
./src/space/space_sas.d \
./src/space/space_sclbp.d \
./src/space/space_setout.d \
./src/space/space_validateInput.d \
./src/space/space_vwtms.d \
./src/space/space_vwtms_inc.d \
./src/space/space_setgaussian.d


# Each subdirectory must supply rules for building sources it contributes
src/space/%.o: ../src/space/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) -g -c -Wall -std=c++11 -fopenmp -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


