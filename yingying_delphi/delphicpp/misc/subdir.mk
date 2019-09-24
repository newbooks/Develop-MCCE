################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/delphicpp/misc/misc_grid_opts.cpp \
../src/delphicpp/misc/misc_interpl.cpp \
../src/delphicpp/misc/misc_timer.cpp 

OBJS += \
./src/delphicpp/misc/misc_grid_opts.o \
./src/delphicpp/misc/misc_interpl.o \
./src/delphicpp/misc/misc_timer.o 

CPP_DEPS += \
./src/delphicpp/misc/misc_grid_opts.d \
./src/delphicpp/misc/misc_interpl.d \
./src/delphicpp/misc/misc_timer.d 


# Each subdirectory must supply rules for building sources it contributes
src/delphicpp/misc/%.o: ../src/delphicpp/misc/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 --std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


