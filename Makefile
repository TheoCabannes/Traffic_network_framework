CXX = g++
SRC = graph_cuckoo_2.cpp dijkstra.cpp road_network.cpp contraction_hierarchy.cpp test.cpp 
OPENMP = -fopenmp
OBJS = $(SRC:%.cpp=%.o)

LIBS = -Wall -lm
TARGET = test

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LIBS) $(OPENMP)

$(OBJS): $(SRC)
	$(CXX) -c $(SRC) $(OPENMP)

clean:
	rm -f $(OBJS) $(TARGET)

