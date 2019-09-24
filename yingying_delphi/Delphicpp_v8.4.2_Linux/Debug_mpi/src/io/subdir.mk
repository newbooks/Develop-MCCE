################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/io/io_epsmap.cpp \
../src/io/io_force.cpp \
../src/io/io_frc.cpp \
../src/io/io_misc.cpp \
../src/io/io_pdb.cpp 

OBJS += \
./src/io/io_epsmap.o \
./src/io/io_force.o \
./src/io/io_frc.o \
./src/io/io_misc.o \
./src/io/io_pdb.o 

CPP_DEPS += \
./src/io/io_epsmap.d \
./src/io/io_force.d \
./src/io/io_frc.d \
./src/io/io_misc.d \
./src/io/io_pdb.d 


# Each subdirectory must supply rules for building sources it contributes
src/io/%.o: ../src/io/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS)  -c -Wall -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


