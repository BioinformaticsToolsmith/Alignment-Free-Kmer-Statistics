
CXX ?= g++
NAME ?= alignment-tool
VERSION ?= 0.9.2
CXXFLAGS ?= -fopenmp -DOPENMP -O3 -march=native -DHIST -g -std=c++11

SOURCES := $(shell find ./src -name '*.cpp')
OBJECTS = $(SOURCES:./src/%.cpp=bin/%.o)

all: bin/$(NAME)

bin/$(NAME): $(OBJECTS)
	$(CXX) $(CXXFLAGS) $^ -o $@ $(LDFLAGS)

bin/%.o: src/%.cpp | bin/align bin/nonltr bin/exception bin/utility
	$(CXX) $(CXXFLAGS) -c $< -o $@
bin:
	mkdir bin
bin/nonltr: | bin
	mkdir bin/nonltr
bin/exception: | bin
	mkdir bin/exception
bin/align: | bin
	mkdir bin/align
bin/utility: | bin
	mkdir bin/utility
clean:
	$(RM) bin/$(NAME)
	$(RM) $(OBJECTS)

dist:
	$(MAKE) clean
	$(shell cd .. && tar czf $(NAME)-v$(VERSION).tar.gz $(NAME) && mv $(NAME)-v$(VERSION).tar.gz $(NAME))
.PHONY: all clean dist
