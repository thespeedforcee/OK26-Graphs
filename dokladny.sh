#!/bin/bash
# Sciezka DOKLADNA: nauty (geng) generuje nieizomorficzne grafy spojne G(n,k),
# nasze sito (weryfikator) sprawdza kazdy z nich. Dziala przez limit_s sekund
# na losowo wybranych fragmentach przestrzeni (geng res/mod), potem raportuje:
# ile grafow przejrzano, ile znaleziono i jakie.
#
# Uzycie: ./dokladny.sh n k limit_s [mod=1048576] [workers=1] [plik_wynikow]

n=$1; k=$2; limit=${3:-60}; mod=${4:-1048576}; workers=${5:-1}
out=${6:-dokladny_${n}_${k}.g6}
[ -z "$n" ] || [ -z "$k" ] && { echo "Uzycie: $0 n k limit_s [mod] [workers] [plik]" >&2; exit 1; }
tmp=$(mktemp -d); end=$(( $(date +%s) + limit ))

worker() {
    local id=$1 total=0 found=0
    while [ $(date +%s) -lt $end ]; do
        local left=$(( end - $(date +%s) )); [ $left -gt 60 ] && left=60
        local res=$(( (RANDOM * 32768 + RANDOM) % mod ))
        local part
        part=$(timeout $left geng -c -q $n $k:$k $res/$mod 2>/dev/null | ./weryfikator -g 2>"$tmp/err.$id" | tee -a "$out" | wc -l)
        local c; c=$(grep -o 'Zweryfikowano [0-9]*' "$tmp/err.$id" | grep -o '[0-9]*$')
        total=$(( total + ${c:-0} )); found=$(( found + part ))
        echo "[worker $id] fragment $res/$mod: razem przejrzano=$total, znaleziono=$found" >&2
    done
    echo "$total $found" > "$tmp/w.$id"
}

for i in $(seq 1 $workers); do worker $i & done
wait
cat "$tmp"/w.* | awk '{t+=$1; f+=$2} END {printf "PODSUMOWANIE (metoda dokladna): przejrzano=%d grafow, znaleziono=%d\n", t, f}'
echo "Znalezione grafy (graph6): $out"; rm -rf "$tmp"