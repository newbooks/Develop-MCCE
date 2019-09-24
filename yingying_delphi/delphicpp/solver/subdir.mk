################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/delphicpp/solver/solver_bndy_isCoulombBndy.cpp \
../src/delphicpp/solver/solver_bndy_isDipolarBndy.cpp \
../src/delphicpp/solver/solver_bndy_isFocusBndy.cpp \
../src/delphicpp/solver/solver_bndy_setBndy.cpp \
../src/delphicpp/solver/solver_fastSOR_initOddEvenItr.cpp \
../src/delphicpp/solver/solver_fastSOR_itit.cpp \
../src/delphicpp/solver/solver_fastSOR_itrEvenPoints.cpp \
../src/delphicpp/solver/solver_fastSOR_itrOddPoints.cpp \
../src/delphicpp/solver/solver_fastSOR_mkdbsf.cpp \
../src/delphicpp/solver/solver_fastSOR_nitit.cpp \
../src/delphicpp/solver/solver_fastSOR_postItr.cpp \
../src/delphicpp/solver/solver_fastSOR_relfac.cpp \
../src/delphicpp/solver/solver_fastSOR_run.cpp \
../src/delphicpp/solver/solver_fastSOR_setcrg.cpp \
../src/delphicpp/solver/solver_fastSOR_validateInput.cpp 

OBJS += \
./src/delphicpp/solver/solver_bndy_isCoulombBndy.o \
./src/delphicpp/solver/solver_bndy_isDipolarBndy.o \
./src/delphicpp/solver/solver_bndy_isFocusBndy.o \
./src/delphicpp/solver/solver_bndy_setBndy.o \
./src/delphicpp/solver/solver_fastSOR_initOddEvenItr.o \
./src/delphicpp/solver/solver_fastSOR_itit.o \
./src/delphicpp/solver/solver_fastSOR_itrEvenPoints.o \
./src/delphicpp/solver/solver_fastSOR_itrOddPoints.o \
./src/delphicpp/solver/solver_fastSOR_mkdbsf.o \
./src/delphicpp/solver/solver_fastSOR_nitit.o \
./src/delphicpp/solver/solver_fastSOR_postItr.o \
./src/delphicpp/solver/solver_fastSOR_relfac.o \
./src/delphicpp/solver/solver_fastSOR_run.o \
./src/delphicpp/solver/solver_fastSOR_setcrg.o \
./src/delphicpp/solver/solver_fastSOR_validateInput.o 

CPP_DEPS += \
./src/delphicpp/solver/solver_bndy_isCoulombBndy.d \
./src/delphicpp/solver/solver_bndy_isDipolarBndy.d \
./src/delphicpp/solver/solver_bndy_isFocusBndy.d \
./src/delphicpp/solver/solver_bndy_setBndy.d \
./src/delphicpp/solver/solver_fastSOR_initOddEvenItr.d \
./src/delphicpp/solver/solver_fastSOR_itit.d \
./src/delphicpp/solver/solver_fastSOR_itrEvenPoints.d \
./src/delphicpp/solver/solver_fastSOR_itrOddPoints.d \
./src/delphicpp/solver/solver_fastSOR_mkdbsf.d \
./src/delphicpp/solver/solver_fastSOR_nitit.d \
./src/delphicpp/solver/solver_fastSOR_postItr.d \
./src/delphicpp/solver/solver_fastSOR_relfac.d \
./src/delphicpp/solver/solver_fastSOR_run.d \
./src/delphicpp/solver/solver_fastSOR_setcrg.d \
./src/delphicpp/solver/solver_fastSOR_validateInput.d 


# Each subdirectory must supply rules for building sources it contributes
src/delphicpp/solver/%.o: ../src/delphicpp/solver/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O3 -g3 -Wall -c -fmessage-length=0 --std=c++0x -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


