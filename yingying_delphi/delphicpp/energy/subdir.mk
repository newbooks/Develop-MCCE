################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/delphicpp/energy/energy_clb.cpp \
../src/delphicpp/energy/energy_clbmedia.cpp \
../src/delphicpp/energy/energy_clbnonl.cpp \
../src/delphicpp/energy/energy_clbtotal.cpp \
../src/delphicpp/energy/energy_nl.cpp \
../src/delphicpp/energy/energy_react.cpp \
../src/delphicpp/energy/energy_run.cpp 

OBJS += \
./src/delphicpp/energy/energy_clb.o \
./src/delphicpp/energy/energy_clbmedia.o \
./src/delphicpp/energy/energy_clbnonl.o \
./src/delphicpp/energy/energy_clbtotal.o \
./src/delphicpp/energy/energy_nl.o \
./src/delphicpp/energy/energy_react.o \
./src/delphicpp/energy/energy_run.o 

CPP_DEPS += \
./src/delphicpp/energy/energy_clb.d \
./src/delphicpp/energy/energy_clbmedia.d \
./src/delphicpp/energy/energy_clbnonl.d \
./src/delphicpp/energy/energy_clbtotal.d \
./src/delphicpp/energy/energy_nl.d \
./src/delphicpp/energy/energy_react.d \
./src/delphicpp/energy/energy_run.d 


# Each subdirectory must supply rules for building sources it contributes
src/delphicpp/energy/%.o: ../src/delphicpp/energy/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 --std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


