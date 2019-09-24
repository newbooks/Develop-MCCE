################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/interface/exceptions.cpp \
../src/interface/interface_datacontainer.cpp \
../src/interface/interface_datamarshal.cpp 

OBJS += \
./src/interface/exceptions.o \
./src/interface/interface_datacontainer.o \
./src/interface/interface_datamarshal.o 

CPP_DEPS += \
./src/interface/exceptions.d \
./src/interface/interface_datacontainer.d \
./src/interface/interface_datamarshal.d 


# Each subdirectory must supply rules for building sources it contributes
src/interface/%.o: ../src/interface/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) -g -c -Wall -std=c++11 -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


