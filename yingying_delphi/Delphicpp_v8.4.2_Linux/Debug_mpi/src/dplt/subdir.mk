################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 

## For DPLT
DARTDIR=../src/dplt/ddm/dart-impl
DDMIDIR=../src/dplt/ddm
DPLTDIR=../src/dplt

## The Source Code to Make
SOURCES= $(DARTDIR)/array  \
		$(DARTDIR)/dart_config  \
		$(DARTDIR)/dart_initialization  \
		$(DARTDIR)/dart_locality  \
		$(DARTDIR)/dart_mem  \
		$(DARTDIR)/dart_synchronization  \
		$(DARTDIR)/dart_team_private  \
		$(DARTDIR)/host_topology  \
		$(DARTDIR)/locality  \
		$(DARTDIR)/papi  \
		$(DARTDIR)/unit_locality  \
		$(DARTDIR)/dart_communication  \
		$(DARTDIR)/dart_globmem  \
		$(DARTDIR)/dart_io_hdf5  \
		$(DARTDIR)/dart_locality_priv  \
		$(DARTDIR)/dart_segment  \
		$(DARTDIR)/dart_team_group  \
		$(DARTDIR)/domain_locality  \
		$(DARTDIR)/hwinfo  \
		$(DARTDIR)/logging  \
		$(DARTDIR)/string  \
		$(DDMIDIR)/BenchmarkParams  \
		$(DDMIDIR)/Distribution  \
		$(DDMIDIR)/Init  \
		$(DDMIDIR)/LocalityDomain  \
		$(DDMIDIR)/Logging  \
		$(DDMIDIR)/StackTrace  \
		$(DDMIDIR)/SUMMA  \
		$(DDMIDIR)/TeamLocality  \
		$(DDMIDIR)/TimestampClockPosix  \
		$(DDMIDIR)/TimestampPAPI  \
		$(DDMIDIR)/TypeInfo  \
		$(DDMIDIR)/Config  \
		$(DDMIDIR)/GlobPtr  \
		$(DDMIDIR)/Locality  \
		$(DDMIDIR)/LocalityJSONPrinter  \
		$(DDMIDIR)/Math  \
		$(DDMIDIR)/StreamConversion  \
		$(DDMIDIR)/Team  \
		$(DDMIDIR)/Timer  \
		$(DDMIDIR)/TimestampCounterPosix  \
		$(DDMIDIR)/Trace  \
		
CPP_SRCS += \
#../src/dplt/dplt.cpp \
../src/dplt/any.cpp 

OBJS += \
#./src/dplt/dplt.o \
./src/dplt/any.o 

CPP_DEPS += \
#./src/dplt/dplt.d \
./src/dplt/any.d 

## The Objects to Generate
DPLTOBJECTS=$(addsuffix .o, $(SOURCES))


# Each subdirectory must supply rules for building sources it contributes
#src/dplt/%.o: ../src/dplt/%.cpp
#	@echo 'Building file: $<'
#	@echo 'Invoking: GCC C++ Compiler'
#	$(CXX) $(CXXFLAGS) $(CXXFLAGS) -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o"$@" "$<"
#	@echo 'Finished building: $<'
#	@echo ' '


