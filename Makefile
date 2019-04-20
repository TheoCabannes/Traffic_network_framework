OBJS = frank_wolfe.o parse_csv.o dijkstra.o

LIBS = -Wall -g
TARGET = frank_wolfe

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LIBS)

$(OBJS): frank_wolfe.cpp parse_csv.cpp dijkstra.cpp
	$(CXX) -c frank_wolfe.cpp parse_csv.cpp dijkstra.cpp

clean:
	rm -f $(OBJS) $(TARGET)
