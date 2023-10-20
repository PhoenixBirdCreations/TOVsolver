CC = gcc
CFLAGS = -Wall -Werror
LIBS = -lm

SRCS = intMain.c EOS_cold.c read_par_tov.c read_par_eos.c TOV_library.c
OBJS = $(SRCS:.c=.o)
EXECUTABLE = solver

.PHONY: all clean

all: $(EXECUTABLE)

$(EXECUTABLE): $(OBJS)
	$(CC) $(CFLAGS) $(OBJS) $(LIBS) -o $@

%.o: %.c
	$(CC) $(CFLAGS) -c $< -o $@

clean:
	rm -f $(EXECUTABLE) $(OBJS)
