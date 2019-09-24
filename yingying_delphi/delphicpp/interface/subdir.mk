################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/delphicpp/interface/exceptions.cpp \
../src/delphicpp/interface/interface_datacontainer.cpp \
../src/delphicpp/interface/interface_datamarshal.cpp 

OBJS += \
./src/delphicpp/interface/exceptions.o \
./src/delphicpp/interface/interface_datacontainer.o \
./src/delphicpp/interface/interface_datamarshal.o 

CPP_DEPS += \
./src/delphicpp/interface/exceptions.d \
./src/delphicpp/interface/interface_datacontainer.d \
./src/delphicpp/interface/interface_datamarshal.d 


# Each subdirectory must supply rules for building sources it contributes
src/delphicpp/interface/%.o: ../src/delphicpp/interface/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 --std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


