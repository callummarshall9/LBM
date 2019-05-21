CC=g++ -std=c++11 -pedantic
CFLAGS=-Wall
HEADERS=headers/*.hpp
SRC=source/*.cpp
OUTPUT=main
hellomake:
	$(CC) $(CFLAGS) -I$(shell pwd) main.cpp $(SRC) -o $(OUTPUT)
