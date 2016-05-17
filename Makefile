
SRC1 = code_circle.c 
SRC2 = code_semicircle.c 
TARGET1 = circle_meshgen
TARGET2 = semicircle_meshgen
HOME = /home/neal
BIN = bin
CC = gcc
FLAGS = -lm

$(TARGET1) : $(SRC1)
	$(CC) $(FLAGS) $(SRC1) -o $(TARGET1)

$(TARGET2) : $(SRC2)
	$(CC) $(FLAGS) $(SRC2) -o $(TARGET2)

install :
	-mv $(TARGET1) $(HOME)/$(BIN)
	-mv $(TARGET2) $(HOME)/$(BIN)

all : $(TARGET1) $(TARGET2)

clean :
	-rm $(TARGET1) $(TARGET2) *.o core
 
