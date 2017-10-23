#!/usr/bin/env bash

#if [ ! -e "published_genomes.html" ]; then
curl http://genome.jgi.doe.gov/fungi/fungi.info.html > published_genomes.html
#fi

cat published_genomes.html \
 | grep -B4 "/pubmed/" \
 | grep "href" \
 | awk '{
     if (NR % 2 == 1){
       split($0,a,"<|>")
       printf "%s\t", a[5]
     } else {
       split($0,a,"/|\"")
       printf "%s\n", a[6]
     }
   }' \
  > published_genomes

cut -f1 published_genomes > published_genomes_names
