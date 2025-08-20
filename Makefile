# Compiler and flags
CXX = g++
# Added -MMD and -MP for automatic dependency generation
CXXFLAGS = -Wall -g -MMD -MP

# Source files
SRCS = main.cpp \
       encoder.cpp

# Object files (derived from source files)
OBJS = $(SRCS:.cpp=.o)

# Dependency files (derived from source files)
DEPS = $(SRCS:.cpp=.d)

# Executable name
TARGET = my_program

# Default target: build the executable
all: $(TARGET)

# Rule to link object files into the executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) $(OBJS) -o $@

# Rule to compile C++ source files into object files
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# Include the generated dependency files
# The hyphen ignores errors if the files don't exist yet
-include $(DEPS)

# Clean rule: remove object files, dependency files, and the executable
clean:
	rm -f $(OBJS) $(DEPS) $(TARGET)