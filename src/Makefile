CXX=g++
CXXFLAGS=-Wall -g -fopenmp
LIBS=-fopenmp

SRCS=main.cpp

OBJS=$(SRCS:.cpp=.o)

MAIN=main

HEADERS=helpers.hpp implementation.hpp

.PHONY: depend clean
all: $(MAIN)
	@echo Compiled the executable $(MAIN)

$(MAIN): $(OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MAIN) $(OBJS) $(LFLAGS) $(LIBS)

%.o: %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	$(RM) *.o *~ $(MAIN)

depend: $(SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE THIS LINE -- make depend needs it

