#!/usr/bin/env sh
cut -d "$DELIM" -f 1,2,4 | sed 's/>//g' | awk -F "$DELIM" '{$3=$3*100;  print $1"'$DELIM'"$2"'$DELIM'"$3; print $2"'$DELIM'"$1"'$DELIM'"$3} ' | uniq | sort -n | column -t -s "$DELIM" -o "$DELIM"
