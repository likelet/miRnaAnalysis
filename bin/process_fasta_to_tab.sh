#!/bin/sh
file=$1
  seqkit fx2tab ${file} | awk '{OFS="\t"}{print $1,"9606",$2}'