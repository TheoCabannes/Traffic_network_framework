OBJS = frank_wolfe.o parse_csv.o

LIBS = -Wall -g
TARGET = frank_wolfe

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(OBJS) -o $(TARGET) $(LIBS)

$(OBJS): frank_wolfe.cpp parse_csv.cpp
	$(CXX) -c frank_wolfe.cpp parse_csv.cpp

clean:
	rm -f $(OBJS) $(TARGET)
