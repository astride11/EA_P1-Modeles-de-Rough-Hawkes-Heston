CXX = g++
CXXFLAGS = -std=c++11 -Wall -Wextra -O2
TARGET = test_lifted
OBJS = Option.o Heston.o Lifted.o test_lifted.o

all: $(TARGET)

$(TARGET): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $(TARGET) $(OBJS)

Option.o: Option.cpp Option.hpp
	$(CXX) $(CXXFLAGS) -c Option.cpp

Heston.o: Heston.cpp Heston.hpp Option.hpp
	$(CXX) $(CXXFLAGS) -c Heston.cpp

Lifted.o: Lifted.cpp Lifted.hpp Option.hpp
	$(CXX) $(CXXFLAGS) -c Lifted.cpp

test_lifted.o: test_lifted.cpp Lifted.hpp Heston.hpp
	$(CXX) $(CXXFLAGS) -c test_lifted.cpp

clean:
	rm -f $(OBJS) $(TARGET)

.PHONY: all clean