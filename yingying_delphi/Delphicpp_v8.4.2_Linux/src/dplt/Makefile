# The Compiler to Use
CC=mpicc 
CXX=mpicxx

## The Target Executable to Make
##TARGET=NArray-Any

## The Compiler Options to Add (if use OpenMP parallel, plese add -fopenmp flag)
CFLAGS=-O3 -DNDEBUG -std=c99 -mtune=native -march=native -mfpmath=sse -DDDM_RELEASE -DDDM_ENABLE_ASSERTIONS -DDART_ENABLE_ASSERTIONS -DDART -DDART_MPI_IMPL_ID="mpich" -DDART_MPI_IMPL_IS_MPICH -DDART_MPI_ENABLE_UNIFIED_MEMORY_MODEL -DDART_MPI_DISABLE_SHARED_WINDOWS

CXXFLAGS= -std=c++11

## The Libraries to Add
LIB=

## The Header Search to Add
INCLUDE=

## The Linker Options to Add
LDFLAGS=

## The Source Code Directories
DARTDIR=./ddm/dart-impl
DDMIDIR=./ddm

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


## The Objects to Generate
OBJECTS=$(addsuffix .o, $(SOURCES))

## Default rule executed
all: dplt
        @true


## Rule for making the actual target

dplt: $(OBJECTS)
        @echo ' '
        @echo 'Building target: $@'
        @echo 'Invoking: C++ Linker'
        $(CXX) dplt.cpp $(CXXFLAGS) $(LDFLAGS) $(LIB) $(INCLUDE) $(OBJECTS) -o $@
        @echo 'Finished building target: $@'
        @echo ' '
        
NArray-Any: $(OBJECTS)
        @echo ' '
        @echo 'Building target: $@'
        @echo 'Invoking: C++ Linker'
        $(CXX) NArray-Any.cpp $(CFLAGS) $(LDFLAGS) $(LIB) $(INCLUDE) $(OBJECTS) -o $@
        @echo 'Finished building target: $@'
        @echo ' '


## Clean Rule
clean:
        @echo "Clean up .o"
        @-rm -f $(TARGET) $(OBJECTS)
