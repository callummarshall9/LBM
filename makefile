CC=g++
CFLAGS=-g -Wall -std=c++11 -pedantic -pthread
HEADERS=headers/*.hpp
SRC=source/*.cpp
OUTPUT=main
LIBS=`pkg-config gtkmm-3.0 --cflags --libs`

hellomake:
	$(CC) $(CFLAGS) -I$(shell pwd) main.cpp $(SRC) -o $(OUTPUT) $(LIBS)
	sh make_buildnum.sh
