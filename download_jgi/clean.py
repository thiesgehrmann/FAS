#!/usr/bin/env python

import csv
import gff2gff3 as g2g3

###############################################################################
# From https://stackoverflow.com/a/3431835

import hashlib

def hash_bytestr_iter(bytesiter, hasher, ashexstr=False):
    for block in bytesiter:
        hasher.update(block)
    return (hasher.hexdigest() if ashexstr else hasher.digest())

def file_as_blockiter(afile, blocksize=65536, maxBlocks=None):
    noMaxBlocks = True if maxBlocks is None else False;
      
    i = 0
    with afile:
        block = afile.read(blocksize)
        while (len(block) > 0) and (noMaxBlocks or (i < maxBlocks)):
            yield block
            i += 1
            block = afile.read(blocksize)


def hashFile(fname, maxBlocks=None):
  return hash_bytestr_iter(file_as_blockiter(open(fname, 'rb'), maxBlocks=maxBlocks), hashlib.sha256())
#edef

#data/Laetiporus_sulphureus_var._sulphureus_v1.0/Laesu1/Laesu1.filtered_proteins.FilteredModels1.gff3.gz

def isGzipped(filename):
  with open(filename, "rb") as ifd:
    bytes = ifd.read(2)
    bytes = [ str(b) for b in bytes ]
    if (len(bytes) > 1) and (bytes[0] == "\x1f") and (bytes[1] == "\x8b") :
      return True
    else:
      return False

###############################################################################

checkSums = {}

with open("download_list", "r") as ifd, open("download_list_cleaned", "w") as ofd:
  reader = csv.reader(ifd, delimiter="\t", quotechar='"')
  for row in reader:
    (longName, shortName, timeStamp, fileName) = row

      # Skip the header line
    if longName[0] == '#':
      continue
    #fi

    if (fileName[-2:] == "gz") and not(isGzipped(fileName)):
      print("File not in GZIP format... Removing %s" % fileName)
      continue
    #fi

      # Remove duplicate files from the list, probably 4 blocks are enough to
      # satistfy this 
    checksum = hashFile(fileName, maxBlocks=4)
    if checksum in checkSums:
      print("Detected duplicate: %s  (exists as %s)" % (fileName, checkSums[checksum]))
      continue
    #fi
    checkSums[checksum] = fileName

    if ("gff3" not in fileName) and ("gff" in fileName):
      G = g2g3.readGTF2(fileName)

      if len(G.keys()) == 0:
        print("Detected Empty GFF file. Skipping %s" % (fileName))
        continue
      #fi

      g2g3.writeGTF2asGFF3(G, fileName + ".gff3.gz")
      print("Converted "+ fileName + " to GFF3")
      ofd.write("%s\t%s\t%s\t%s\n" % (longName, shortName, timeStamp, fileName + ".gff3.gz"))
    elif "gff" not in fileName:
      ofd.write("%s\t%s\t%s\t%s\n" % (longName, shortName, timeStamp, fileName))
    #fi
  #efor
#ewith

