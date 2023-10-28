# Source files in the current directory
SRC_FILES := $(wildcard *.cu)
# Object files derived from source files
OBJ := $(SRC_FILES:.cu=.o)
# The name of the output binary
TARGET := app

# Compiler and compiler options
NVCC := nvcc
NVCCFLAGS := -I.

# Rule to build the target binary
$(TARGET): $(OBJ)
	$(NVCC) $^ -o $@

# Rule to build object files from corresponding source files
%.o: %.cu
	$(NVCC) -x cu $(NVCCFLAGS) -dc $< -o $@

.PHONY: clean
clean:
	rm -f $(OBJ) $(TARGET)

