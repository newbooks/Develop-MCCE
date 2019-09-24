################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/solver/solver_bndy_isCoulombBndy.cpp \
../src/solver/solver_bndy_isDipolarBndy.cpp \
../src/solver/solver_bndy_isFocusBndy.cpp \
../src/solver/solver_bndy_setBndy.cpp \
../src/solver/solver_fastSOR_initOddEvenItr.cpp \
../src/solver/solver_fastSOR_itit.cpp \
../src/solver/solver_fastSOR_itrEvenPoints.cpp \
../src/solver/solver_fastSOR_itrOddPoints.cpp \
../src/solver/solver_fastSOR_mkdbsf.cpp \
../src/solver/solver_fastSOR_nitit.cpp \
../src/solver/solver_fastSOR_postItr.cpp \
../src/solver/solver_fastSOR_relfac.cpp \
../src/solver/solver_fastSOR_run.cpp \
../src/solver/solver_fastSOR_setcrg.cpp \
../src/solver/solver_fastSOR_validateInput.cpp 

OBJS += \
./src/solver/solver_bndy_isCoulombBndy.o \
./src/solver/solver_bndy_isDipolarBndy.o \
./src/solver/solver_bndy_isFocusBndy.o \
./src/solver/solver_bndy_setBndy.o \
./src/solver/solver_fastSOR_initOddEvenItr.o \
./src/solver/solver_fastSOR_itit.o \
./src/solver/solver_fastSOR_itrEvenPoints.o \
./src/solver/solver_fastSOR_itrOddPoints.o \
./src/solver/solver_fastSOR_mkdbsf.o \
./src/solver/solver_fastSOR_nitit.o \
./src/solver/solver_fastSOR_postItr.o \
./src/solver/solver_fastSOR_relfac.o \
./src/solver/solver_fastSOR_run.o \
./src/solver/solver_fastSOR_setcrg.o \
./src/solver/solver_fastSOR_validateInput.o 

CPP_DEPS += \
./src/solver/solver_bndy_isCoulombBndy.d \
./src/solver/solver_bndy_isDipolarBndy.d \
./src/solver/solver_bndy_isFocusBndy.d \
./src/solver/solver_bndy_setBndy.d \
./src/solver/solver_fastSOR_initOddEvenItr.d \
./src/solver/solver_fastSOR_itit.d \
./src/solver/solver_fastSOR_itrEvenPoints.d \
./src/solver/solver_fastSOR_itrOddPoints.d \
./src/solver/solver_fastSOR_mkdbsf.d \
./src/solver/solver_fastSOR_nitit.d \
./src/solver/solver_fastSOR_postItr.d \
./src/solver/solver_fastSOR_relfac.d \
./src/solver/solver_fastSOR_run.d \
./src/solver/solver_fastSOR_setcrg.d \
./src/solver/solver_fastSOR_validateInput.d 


# Each subdirectory must supply rules for building sources it contributes
src/solver/%.o: ../src/solver/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	$(CXX) -O3 -c -Wall -std=c++11 -o"$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


