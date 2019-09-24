################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/delphicpp/space/grid_space.cpp \
../src/delphicpp/space/space_crgarr.cpp \
../src/delphicpp/space/space_cube.cpp \
../src/delphicpp/space/space_eps.cpp \
../src/delphicpp/space/space_indver.cpp \
../src/delphicpp/space/space_msrf.cpp \
../src/delphicpp/space/space_run.cpp \
../src/delphicpp/space/space_sas.cpp \
../src/delphicpp/space/space_sclbp.cpp \
../src/delphicpp/space/space_setout.cpp \
../src/delphicpp/space/space_validateInput.cpp \
../src/delphicpp/space/space_vwtms.cpp \
../src/delphicpp/space/space_vwtms_inc.cpp \
../src/delphicpp/space/space_setgaussian.cpp 

OBJS += \
./src/delphicpp/space/grid_space.o \
./src/delphicpp/space/space_crgarr.o \
./src/delphicpp/space/space_cube.o \
./src/delphicpp/space/space_eps.o \
./src/delphicpp/space/space_indver.o \
./src/delphicpp/space/space_msrf.o \
./src/delphicpp/space/space_run.o \
./src/delphicpp/space/space_sas.o \
./src/delphicpp/space/space_sclbp.o \
./src/delphicpp/space/space_setout.o \
./src/delphicpp/space/space_validateInput.o \
./src/delphicpp/space/space_vwtms.o \
./src/delphicpp/space/space_vwtms_inc.o \
./src/delphicpp/space/space_setgaussian.o

CPP_DEPS += \
./src/delphicpp/space/grid_space.d \
./src/delphicpp/space/space_crgarr.d \
./src/delphicpp/space/space_cube.d \
./src/delphicpp/space/space_eps.d \
./src/delphicpp/space/space_indver.d \
./src/delphicpp/space/space_msrf.d \
./src/delphicpp/space/space_run.d \
./src/delphicpp/space/space_sas.d \
./src/delphicpp/space/space_sclbp.d \
./src/delphicpp/space/space_setout.d \
./src/delphicpp/space/space_validateInput.d \
./src/delphicpp/space/space_vwtms.d \
./src/delphicpp/space/space_vwtms_inc.d \
./src/delphicpp/space/space_setgaussian.d


# Each subdirectory must supply rules for building sources it contributes
src/delphicpp/space/%.o: ../src/delphicpp/space/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 --std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


