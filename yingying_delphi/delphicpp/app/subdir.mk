################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/delphicpp/app/app_mcce.cpp 

OBJS += \
./src/delphicpp/app/app_mcce.o 

CPP_DEPS += \
./src/delphicpp/app/app_mcce.d 


# Each subdirectory must supply rules for building sources it contributes
src/delphicpp/app/%.o: ../src/delphicpp/app/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 --std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


