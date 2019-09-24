################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/delphicpp/delphi/delphi_data_flash.cpp \
../src/delphicpp/delphi/delphi_data_reset.cpp \
../src/delphicpp/delphi/delphi_data_setMap.cpp \
../src/delphicpp/delphi/delphi_data_showMap.cpp \
../src/delphicpp/delphi/delphi_datamarshal_getFunction.cpp \
../src/delphicpp/delphi/delphi_datamarshal_getStatement.cpp \
../src/delphicpp/delphi/delphi_datamarshal_setDefault.cpp \
../src/delphicpp/delphi/delphi_datamarshal_showParameters.cpp \
../src/delphicpp/delphi/delphi_datamarshal_updateParameters.cpp 

OBJS += \
./src/delphicpp/delphi/delphi_data_flash.o \
./src/delphicpp/delphi/delphi_data_reset.o \
./src/delphicpp/delphi/delphi_data_setMap.o \
./src/delphicpp/delphi/delphi_data_showMap.o \
./src/delphicpp/delphi/delphi_datamarshal_getFunction.o \
./src/delphicpp/delphi/delphi_datamarshal_getStatement.o \
./src/delphicpp/delphi/delphi_datamarshal_setDefault.o \
./src/delphicpp/delphi/delphi_datamarshal_showParameters.o \
./src/delphicpp/delphi/delphi_datamarshal_updateParameters.o 

CPP_DEPS += \
./src/delphicpp/delphi/delphi_data_flash.d \
./src/delphicpp/delphi/delphi_data_reset.d \
./src/delphicpp/delphi/delphi_data_setMap.d \
./src/delphicpp/delphi/delphi_data_showMap.d \
./src/delphicpp/delphi/delphi_datamarshal_getFunction.d \
./src/delphicpp/delphi/delphi_datamarshal_getStatement.d \
./src/delphicpp/delphi/delphi_datamarshal_setDefault.d \
./src/delphicpp/delphi/delphi_datamarshal_showParameters.d \
./src/delphicpp/delphi/delphi_datamarshal_updateParameters.d 


# Each subdirectory must supply rules for building sources it contributes
src/delphicpp/delphi/%.o: ../src/delphicpp/delphi/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 --std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


