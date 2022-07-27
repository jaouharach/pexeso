##
CC=gcc
LDFLAGS := -lm -lgsl -lgslcblas -lmcheck
CFLAGS_INC := 
CFLAGS := -g -O2
#

## directories
SRC_DIR :=src/
BIN_DIR :=bin/
OBJ_DIR :=obj/
#

## files
SRCS := $(wildcard $(SRC_DIR)*.c)
PRGS := $(patsubst %.c,%,$(SRCS))
OBJS := $(patsubst $(SRC_DIR)%,$(OBJ_DIR)%.o,$(PRGS))
BIN	 :=$(BIN_DIR)pexeso
#

all: $(OBJS)
	$(CC) $^  $(CFLAGS) $(LDFLAGS) -o $(BIN)
	$(info Binary file $(BIN) has been generated.)

$(OBJ_DIR)%.o: $(SRC_DIR)%.c
	$(CC) -o $@ -c $<

clean:
	rm -f $(OBJ_DIR)*.o $(BIN)