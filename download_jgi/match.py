#!/usr/bin/env python

import csv
import gff2gff3 as g2g3
import gzip
import json

###############################################################################

def loadFastaLengths(fastaFile):

  F = {'': 0}

  current_seq = ""
  with (gzip.open(fastaFile, "r") if fastaFile[-2:] == "gz" else open(fastaFile, "r")) as fd:
    for line in fd:
      line = line.strip()
      if len(line) == 0:
        continue
      #fi
      if line[0] == '>':
        current_seq = line[1:].split(' ')[0]
        F[current_seq] = 0
      else:
        F[current_seq] = F[current_seq] + len(line.strip())
      #fi
  #ewith
  return F
#edef

###############################################################################

def verify_match(G, F):

  errors = 0

  for e in G:
    if e.seqname not in F:
      errors += 1
    else:
      seqlen = F[e.seqname]
      if (int(e.start) > seqlen) or (int(e.end) > seqlen):
        errors += 1
      #fi
    #fi 
  #efor
  return errors
#edef

###############################################################################
if __name__ == "__main__":

  # Read the data
  F = {}
  with open("download_list_cleaned", "r") as ifd:
    reader = csv.reader(ifd, delimiter="\t", quotechar='"')
    for row in reader:
      (longName, shortName, timeStamp, fileName) = row
      if (longName, shortName) not in F:
        F[(longName, shortName)] = [row]
      else:
        F[(longName, shortName)] = F[(longName, shortName)]+ [row]
      #fi
    #efor
  #ewith


  filePairs = []
  noMatch   = []

  for (longName, shortName) in F.keys():
    gffs = [ f for f in F[(longName, shortName)] if "gff3" in f[-1] ]
    fas  = [ f for f in F[(longName, shortName)] if "gff3" not in f[-1] ]
 
    if (len(gffs) == 1) and (len(fas) == 1):
      fasta = loadFastaLengths(fas[0][-1])
      genes = g2g3.readGFF3(gffs[0][-1])
      errors = verify_match(genes, fasta)
      if errors == 0:
        print("Match verified for %s!" % shortName)
        filePairs +=[(longName, shortName, fas[0][-1], gffs[0][-1])]
    else:
      fastas = sorted([ f + [loadFastaLengths(f[-1])] for f in fas], key=lambda x: int(x[2]))[::-1]
      genes  = sorted([ f + [g2g3.readGFF3(f[-1])] for f in gffs], key=lambda x: int(x[2]))[::-1]

      matched = False
      for (longName, shortName, timestamp_fa, fileName_fa, fasta_data) in fastas:
        if matched:
          break
        for (longName, shortName, timestamp_gff, fileName_gff, gff_data) in genes:
          if "_mito" in fileName_gff:
            continue
          #fi
          errors = verify_match(gff_data, fasta_data)
          if errors == 0:
            print("Match verified for %s!" % shortName)
            filePairs += [(longName, shortName, fileName_fa, fileName_gff)]
            matched = True
            break
          #fi
        #efor
      #efor
      noMatch += [(longName, shortName)]
    #fi
  #efor

  print("\n#####################################################")
  print("Matches found for %d genomes." % len(filePairs))
  print("Data is found in '%s'." % "data.json")
  print("No match found for %d genomes:" % len(noMatch))
  for (i, (nm_l, nm_s)) in enumerate(noMatch):
    print("%d: %s, %s" % (i+1, nm_l, nm_s))
  #efor

  J = { "%d_%s" % (i, shortName) : { "name": longName, "shortname": shortName, "genome": fasta, "gff": gff } for (i, (longName, shortName, fasta, gff)) in enumerate(filePairs) }

  with open("data.json", "w") as ofd:
    res = json.dump(J, ofd, sort_keys=True, indent=4, separators=(',', ': '))
  #ewith

#fi
