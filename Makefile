# Should debugging be enabled?
DEBUG = 0

OBJECTS = main.o transaction_minimizer.o
EXEC = transaction_minimizer
CXXFLAGS = -Wall -O3 -DDEBUG=$(DEBUG)
CC = g++
RM = rm 

all: $(OBJECTS)
	$(CC) $(FLAGS) -o $(EXEC) $(OBJECTS)

clean:
	$(RM) -f $(OBJECTS) $(EXEC)
