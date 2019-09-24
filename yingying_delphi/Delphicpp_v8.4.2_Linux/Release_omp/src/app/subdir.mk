################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/app/app_delphi.cpp 

OBJS += \
./src/app/app_delphi.o 

CPP_DEPS += \
./src/app/app_delphi.d 


# Each subdirectory must supply rules for building sources it contributes
src/app/%.o: ../src/app/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) -O3 -c -Wall -std=c++11 -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


