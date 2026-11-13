#!/bin/bash
# Postac kanoniczna grafow (nauty): zamienia kazdy graph6 na kanoniczny graph6 i usuwa
# duplikaty izomorficzne. Uzycie: ./kanon.sh < wynik.g6 > unikalne.g6
sed 's/\r//' | labelg -q | sort -u