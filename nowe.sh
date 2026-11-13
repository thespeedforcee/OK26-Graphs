#!/bin/bash
# Wypisuje grafy z pliku "znalezione" (graph6), ktorych NIE MA w liscie "znane"
# (porownanie z dokladnoscia do izomorfizmu: obie listy przechodza przez labelg).
# Uzycie: ./nowe.sh znalezione.g6 grafy18L.graph6
a=$(mktemp); b=$(mktemp)
sed 's/\r//' "$1" | labelg -q | sort -u > "$a"
sed 's/\r//' "$2" | labelg -q | sort -u > "$b"
echo "znalezione (rozne): $(wc -l < "$a"), znane: $(wc -l < "$b")" >&2
comm -23 "$a" "$b"
rm -f "$a" "$b"