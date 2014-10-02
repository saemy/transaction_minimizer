# Should debugging be enabled?
DEBUG ?= 0

OBJECTS = main.o transaction_minimizer.o
EXEC = transaction_minimizer
CXXFLAGS = -Wall
CXX = g++
RM = rm -f

ifeq ($(DEBUG), 0)
    CXXFLAGS += -O3
else
    CXXFLAGS += -O0 -g -DDEBUG=1
endif

$(EXEC): $(OBJECTS)
	$(CXX) $(FLAGS) -o $@ $(OBJECTS)

all: $(EXEC)

clean:
	$(RM) $(OBJECTS) $(EXEC)
