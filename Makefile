# 1. Variables (Makes it easy to change things later)
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall
TARGET = sim3D.x
SRCS = simulation3D.cpp verlet.cpp observables.cpp
OBJS = $(SRCS:.cpp=.o)

# 2. The default rule (What happens when you just type 'make')
all: $(TARGET)

# 3. How to build the final executable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# 4. How to compile each individual .cpp file into an .o (object) file
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 5. A cleanup rule to delete compiled files so you can start fresh
clean:
	rm -f $(OBJS) $(TARGET)
