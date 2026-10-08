# 1. Detección de Sistema Operativo (Magia multiplataforma)
ifeq ($(OS),Windows_NT)
    TARGET = sim3D.exe
    CLEAN_CMD = del /Q /F *.o $(TARGET)
else
    TARGET = sim3D.x
    CLEAN_CMD = rm -f $(OBJS) $(TARGET)
endif

# 2. Variables
CXX = g++
CXXFLAGS = -std=c++17 -O3 -Wall
SRCS = simulation3D.cpp verlet.cpp observables.cpp
OBJS = $(SRCS:.cpp=.o)

# 3. Regla por defecto
all: $(TARGET)

# 4. Construir el ejecutable
$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

# 5. Compilar cada .cpp a .o
%.o: %.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

# 6. Limpieza multiplataforma
clean:
	$(CLEAN_CMD)