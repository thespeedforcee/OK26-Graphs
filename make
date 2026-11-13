CC=gcc
CFLAGS=-O3 -Wall -Wextra
BINS=generator algorytm_zachlanny algorytm_zachlanny_ulosowiony metaheurystyka_wyzarzanie weryfikator rozszerzanie_grafu
all: $(BINS)
%: %.c integral.h
	$(CC) $(CFLAGS) -o $@ $< -lm
clean:
	rm -f $(BINS)