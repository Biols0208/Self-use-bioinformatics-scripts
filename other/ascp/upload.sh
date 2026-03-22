#!/usr/bin/env bash

D="subasp@upload.ncbi.nlm.nih.gov:uploads/XXXX/YYYYY"   # YYYYY can be customized. XXXX from NCBI

A="/soft/.aspera/sdk/ascp"
K="/software/ascp/aspera.openssh"                       # aspera.openssh from NCBI
R="120m"; S=10
for F in "$@"; do
  [ -f "$F" ] || continue
  until "$A" -i "$K" -QT -l"$R" -k1 "$F" "$D"; do
    sleep "$S"
  done
done
