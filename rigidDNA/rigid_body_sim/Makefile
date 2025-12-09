CXX = g++
CXXFLAGS = -O3 -std=c++11
TARGET = RigidBodySim
SRC = RigidBodySim.cpp

all: $(TARGET)

$(TARGET): $(SRC)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(SRC)

clean:
	rm -f $(TARGET)
