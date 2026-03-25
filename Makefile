CC = gcc

CFLAGS = -Wall -O2 -std=c11

LDFLAGS = -lm

SRCS = BH.c

OBJS = $(SRCS:.c=.o)

EXEC = bh

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CC) $(CFLAGS) -o $@ $(OBJS) $(LDFLAGS)

%.o: %.c
	$(CC) $(CFLAGS) -c $<

clean:
	rm -f $(OBJS) $(EXEC)

run: $(EXEC)
	./$(EXEC) > output.txt

rebuild: clean all
