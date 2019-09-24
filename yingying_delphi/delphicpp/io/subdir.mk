################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/delphicpp/io/io_epsmap.cpp \
../src/delphicpp/io/io_force.cpp \
../src/delphicpp/io/io_frc.cpp \
../src/delphicpp/io/io_misc.cpp \
../src/delphicpp/io/io_pdb.cpp 

OBJS += \
./src/delphicpp/io/io_epsmap.o \
./src/delphicpp/io/io_force.o \
./src/delphicpp/io/io_frc.o \
./src/delphicpp/io/io_misc.o \
./src/delphicpp/io/io_pdb.o 

CPP_DEPS += \
./src/delphicpp/io/io_epsmap.d \
./src/delphicpp/io/io_force.d \
./src/delphicpp/io/io_frc.d \
./src/delphicpp/io/io_misc.d \
./src/delphicpp/io/io_pdb.d 


# Each subdirectory must supply rules for building sources it contributes
src/delphicpp/io/%.o: ../src/delphicpp/io/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 --std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


