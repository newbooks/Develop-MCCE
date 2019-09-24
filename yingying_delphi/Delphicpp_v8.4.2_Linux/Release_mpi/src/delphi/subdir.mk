################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/delphi/delphi_dplt_updateParameters.cpp \
../src/delphi/delphi_data_flash.cpp \
../src/delphi/delphi_data_reset.cpp \
../src/delphi/delphi_data_setMap.cpp \
../src/delphi/delphi_data_showMap.cpp \
../src/delphi/delphi_datamarshal_getFunction.cpp \
../src/delphi/delphi_datamarshal_getStatement.cpp \
../src/delphi/delphi_datamarshal_setDefault.cpp \
../src/delphi/delphi_datamarshal_showParameters.cpp \
../src/delphi/delphi_datamarshal_updateParameters.cpp \

OBJS += \
./src/delphi/delphi_dplt_updateParameters.o \
./src/delphi/delphi_data_flash.o \
./src/delphi/delphi_data_reset.o \
./src/delphi/delphi_data_setMap.o \
./src/delphi/delphi_data_showMap.o \
./src/delphi/delphi_datamarshal_getFunction.o \
./src/delphi/delphi_datamarshal_getStatement.o \
./src/delphi/delphi_datamarshal_setDefault.o \
./src/delphi/delphi_datamarshal_showParameters.o \
./src/delphi/delphi_datamarshal_updateParameters.o \

CPP_DEPS += \
./src/delphi/delphi_dplt_updateParameters.d \
./src/delphi/delphi_data_flash.d \
./src/delphi/delphi_data_reset.d \
./src/delphi/delphi_data_setMap.d \
./src/delphi/delphi_data_showMap.d \
./src/delphi/delphi_datamarshal_getFunction.d \
./src/delphi/delphi_datamarshal_getStatement.d \
./src/delphi/delphi_datamarshal_setDefault.d \
./src/delphi/delphi_datamarshal_showParameters.d \
./src/delphi/delphi_datamarshal_updateParameters.d \


# Each subdirectory must supply rules for building sources it contributes
src/delphi/%.o: ../src/delphi/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) $(CXXFLAGS) -O3 -c -Wall -std=c++11 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


