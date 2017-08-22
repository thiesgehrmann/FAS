#!/usr/bin/env python

import csv
import shlex
import sys
import gzip

class GTF2Entry(object):

  seqname = ""
  source  = ""
  feature = ""
  start   = ""
  end     = ""
  score   = ""
  strand  = ""
  frame   = ""
  attrRaw = ""
  attr    = {}

  def __parseAttr__(self):
    self.attr = dict([tuple(shlex.split(x.strip())[0:2]) for x in self.attrRaw.split(";") ])
  #edef

  def __init__(self, fileRow):
    (self.seqname, self.source, self.feature, self.start, self.end, self.score, self.strand, self.frame, self.attrRaw) = fileRow[:9]
    self.__parseAttr__()
  #edef

  def getAttr(self, field):
    if field in self.attr:
      return self.attr[field]
    else:
      return ""
    #fi

#eclass

##############################################################################

def readGTF2(filename):

  G = {}

  with ( gzip.open(filename, "r") if filename[-2:] == "gz" else open(filename, "r")) as gffFile:
    gffReader = csv.reader(gffFile, delimiter="\t", quotechar='"')
    for row in gffReader:
      if len(row) < 9:
        continue
      #fi
      if row[0].strip()[0] == '#':
        continue
      #fi

      entry = GTF2Entry(row)
      entryName = entry.getAttr("name")
      if entryName == "":
        entryName = entry.getAttr("gene") # Try for mitocondrial genes
      #fi
      if not(entryName in G):
        G[entryName] = [entry]
      else:
        G[entryName] = G[entryName] + [entry]
      #fi
    #efor
  #ewith
  return G
#edef

###############################################################################

class GFF3Entry(GTF2Entry):
  def __parseAttr__(self):
    self.attr = dict([tuple(x.strip().split('=')[0:2]) for x in self.attrRaw.split(";") ])
  #edef
#eclass

###############################################################################

def readGFF3(filename):
  G = []
  with ( gzip.open(filename, "r") if filename[-2:] == "gz" else open(filename, "r")) as gffFile:
    gffReader = csv.reader(gffFile, delimiter="\t", quotechar='"')
    for row in gffReader:
      if len(row) != 9:
        continue
      #fi
      entry = GFF3Entry(row)
      G.append(entry)
      #fi
    #efor
  #ewith
  return G
#edef

###############################################################################

def writeGFF3(GFF3, gff3File):
  with ( gzip.open(gff3File, "wb") if gff3File[-2:] == "gz" else open(gff3File, "w")) as fd:
    fd.write("##gff-version 3\n")
    for e in GFF3:
      fd.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (e.seqname, e.source, e.feature, e.start, e.end, e.score, e.strand, e.frame, ';'.join([ "%s=%s" % (k,v) for (k,v) in e.attr ])))
    #efor
  #ewith
  fd.close()
#edef

###############################################################################

def writeGTF2asGFF3(GTF2, gff3File):
  with ( gzip.open(gff3File, "wb") if gff3File[-2:] == "gz" else open(gff3File, "w")) as fd:
    fd.write("##gff-version 3\n")
    for gene in GTF2.keys():
      entries = GTF2[gene]
      seqname  = entries[0].seqname
      source = entries[0].source
      start  = min([ int(e.start) for e in entries])
      end    = max([ int(e.end) for e in entries])
      strand = entries[0].strand
      proteinIDs    = [ p for p in [ e.getAttr("proteinId") for e in entries] if len(p) > 0]
      transcriptIDs = [ t for t in [ e.getAttr("transcriptId") for e in entries] if len(t) > 0]

      geneID       = gene
      proteinID    = "p" + geneID if len(proteinIDs) == 0 else proteinIDs[0]
      transcriptID = "t" + geneID if len(transcriptIDs) == 0 else transcriptIDs[0]
  
        # Print gene information
      fd.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seqname, source, "Gene", start, end, ".", strand, ".", "ID=%s;Name=%s;proteinId=%s" % (geneID,geneID,proteinID)))
      fd.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seqname, source, "mRNA", start, end, ".", strand, ".", "ID=%s;Parent=%s" % (transcriptID, geneID)))
      for (i,e) in enumerate(sorted(entries, key=(lambda el: int(el.start) if strand == '+' else -int(el.start)))):
        if e.feature.lower() == "gene":
          continue
        #fi
        fd.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (seqname, source, e.feature, e.start, e.end, e.score, strand, e.frame, "ID=%s.%s.%d;Parent=%s" % (transcriptID, e.feature, i, transcriptID) ))
      #efor 
    #efor
  fd.close()
#edef

if __name__ == "__main__":
  gtf2 = sys.argv[1]
  gff3 = sys.argv[2]

  G = readGTF2(gtf2)
  writeGTF2asGFF3(G, gff3)
