CC = g++
ROOT = `root-config --cflags --glibs --ldflags`
SRC=$(shell ls ./src/*.cpp)
OBJ=$(SRC:.cpp=.o)
.PHONY: clean
.PHONY: clean-all
all: analyze 
	cat include/woa.txt
clean:
	rm -f ./src/*.o
clean-all:
	rm -f ./src/*.o analyze
%.o: %.cpp
	$(CC) -c $(ROOT) -o $@ $<
analyze: $(OBJ) 
	$(CC) $(ROOT) -o $@ $^ 
