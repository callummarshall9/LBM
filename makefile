CC=g++
CFLAGS=-g -Wall -std=c++11 -pedantic -pthread
LBM_HEADERS=LBM/headers/*.hpp
LBM_SRC=LBM/source/*.cpp
CMD_OUTPUT=main_console
GUI_HEADERS=gui/headers/*.hpp
GUI_SRC=gui/source/*.cpp
GUI_OUTPUT=main_gui
LIBS=`pkg-config gtkmm-3.0 --cflags --libs`

hellomake:
	$(CC) $(CFLAGS) -I$(shell pwd) cmdline/main.cpp $(LBM_SRC) -o $(CMD_OUTPUT) $(LIBS)
	sh make_buildnum.sh

gui:
	$(CC) $(CFLAGS) -I$(shell pwd) gui/main.cpp $(GUI_SRC) $(LBM_SRC) -o $(GUI_OUTPUT) $(LIBS)
	sh make_buildnum.sh
.PHONY: hellomake gui
clean:
	rm main_gui main_console
