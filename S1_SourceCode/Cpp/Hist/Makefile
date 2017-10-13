TARGET ?= bin/aligner
VERSION ?= 0.2.1
CXX ?= g++
OBJECTS := $(patsubst src/%.cpp, %.o, $(wildcard src/*.cpp))

CXXFLAGS ?= -fopenmp -Wall -Wextra -std=c++11
ifeq ($(debug), yes)
	CXXFLAGS +=  -ggdb
else
	CXXFLAGS += -O3 -march=native
endif
LDFLAGS ?= -lm
default: $(TARGET)

$(TARGET): $(OBJECTS)
	mkdir -p bin
	$(CXX) $(CXXFLAGS) $(OBJECTS) $(LDFLAGS) -o $@

%.o: src/%.cpp
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) $(TARGET) $(OBJECTS)

dist: clean
	$(shell cd .. && tar czf $(TARGET)-$(VERSION).tar.gz $(TARGET) && cp $(TARGET)-$(VERSION).tar.gz $(TARGET))
.PHONY: clean default dist
