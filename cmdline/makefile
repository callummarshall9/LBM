CC=g++
CFLAGS=-g -Wall -std=c++11 -pedantic
HEADERS=headers/*.hpp
SRC=source/*.cpp
OUTPUT=main

hellomake:
	$(CC) $(CFLAGS) -I$(shell pwd) main.cpp $(SRC) -o $(OUTPUT)
	sh make_buildnum.sh
