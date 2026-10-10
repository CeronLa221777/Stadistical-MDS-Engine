# 1. OS Detection
ifeq ($(OS),Windows_NT)
    TARGET = sim3D.exe
    CLEAN_CMD = del /Q /F *.o $(TARGET)
    PYTHON_CMD = python
    RUN_CMD = .\$(TARGET)
else
    TARGET = sim3D.x
    CLEAN_CMD = rm -f $(OBJS) $(TARGET)
    PYTHON_CMD = python3
    RUN_CMD = ./$(TARGET)
endif

# 2. Variables
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall -fopenmp
SRCS = simulation3D.cpp verlet.cpp observables.cpp
OBJS = $(SRCS:.cpp=.o)

# 3. Default target
all: $(TARGET)

# 4. Build executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# 5. Compile objects
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 6. Run simulation
run: $(TARGET)
	$(RUN_CMD)

# 7. Generate plots
plot:
	$(PYTHON_CMD) PlotterMaster.py

# 8. Clean
clean:
	$(CLEAN_CMD)