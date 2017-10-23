import json
import utils
import csv

import random
import inspect, os
__INSTALL_DIR__ = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

configfile: "config.json"
dconfig = json.load(open(config["genomes"], "r"))
OUTDIR = config["outdir"]
__RUN_DIR__ = "%s/run" % OUTDIR

###############################################################################

GENOMES = dconfig["data"].keys()

###############################################################################
# BLAST                                                                       #
###############################################################################

__BLAST_OUTDIR__ = "%s/blast" % __RUN_DIR__
__ANALYSIS_OUTDIR__ = "%s/analysis" % __RUN_DIR__
__ITOL_OUTDIR__ = "%s/itol" % __RUN_DIR__


__BLASTFIELDS__  = "qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qlen"

###############################################################################
# BLAST                                                                       #
###############################################################################


rule make_blastdb:
  input:
    prots = lambda wildcards: "%s/%s" % (dconfig["dataprefix"], dconfig["data"][wildcards.genome]["aa"])
  output:
    db = "%s/blastdb.{genome}.db" % __BLAST_OUTDIR__
  params:
    blastfields = __BLASTFIELDS__
  shell: """
    mkdir -p `dirname {output.db}`
    zcat {input.prots} > {output.db}.input
    makeblastdb -dbtype prot -in {output.db}.input -out {output.db}
    touch {output.db}
  """

rule run_query_blast:
  input:
    query = config["prot_query"],
    db    = lambda wildcards: "%s/blastdb.%s.db" % (__BLAST_OUTDIR__, wildcards.genome)
  output:
    hits = "%s/blastoutput.{genome}.tsv" % __BLAST_OUTDIR__
  threads: 1
  params:
    blastfields = __BLASTFIELDS__
  shell: """
    blastp -query {input.query} -db {input.db} -outfmt "6 {params.blastfields}" -out {output.hits} -num_threads {threads}
  """

rule all_blasts:
  input:
    hits = expand("%s/blastoutput.{genome}.tsv" % (__BLAST_OUTDIR__), genome=GENOMES)

#http://www.metagenomics.wiki/tools/blast/blastn-output-format-6
# 1.	 qseqid	 query (e.g., gene) sequence id
# 2.	 sseqid	 subject (e.g., reference genome) sequence id
# 3.	 pident	 percentage of identical matches
# 4.	 length	 alignment length
# 5.	 mismatch	 number of mismatches
# 6.	 gapopen	 number of gap openings
# 7.	 qstart	 start of alignment in query
# 8.	 qend	 end of alignment in query
# 9.	 sstart	 start of alignment in subject
# 10.	 send	 end of alignment in subject
# 11.	 evalue	 expect value
# 12.	 bitscore	 bit score
# 13.    slen
# 14.    qlen

###############################################################################
# Reconstruct the fragments                                                   #
###############################################################################

rule filter_hits:
  input:
    hits    = lambda wildcards: "%s/blastoutput.%s.tsv" % (__BLAST_OUTDIR__, wildcards.genome),
    queries = config["prot_query"]
  output:
    hits = "%s/{genome}/filtered_hits.tsv" % __ANALYSIS_OUTDIR__
  params:
    blastfields = __BLASTFIELDS__
  run:
    queries = utils.loadFasta(input.queries)
    hits    = utils.readBlastFile(input.hits, params.blastfields)

    # Filter by e-value, threshold
    hits = [ h for h in hits if h.evalue < config["evalue_threshold"] ]

    # As we probably use multiple queries, we need to assign hit genes to one query or the other (otherwise we count them twice)
    # To do this, I select the query for which we have the highest coverage originating from a single subject.
    # For example:
    #  X -> Y -------||||||||--
    #  X -> Y -||||||----------
    #  X -> Z ---|||||||||||---
    # In this case, I will select Y as the best query, because I cover more of it with X than I cover Y
    filtered_hits = []
    for subject in set([ h.sseqid for h in hits ]):
      subject_hits = [ h for h in hits if h.sseqid == subject ]
      subject_queries = set([ h.qseqid for h in subject_hits ])
      if len(subject_queries) > 1:
        hits_per_query = { query: [h for h in subject_hits if h.qseqid == query] for query in subject_queries }
        # Pick hits for the query for which we have the largest coverage
        best_query_for_subject = max(subject_queries, key=lambda q: utils.region_covered([ (h.qstart, h.qend) for h in hits_per_query[q] ], len(queries[q])).mean())
        subject_query_hits = hits_per_query[best_query_for_subject]
        print("%s, %s, %s" % (wildcards.genome, subject, best_query_for_subject))
        filtered_hits += subject_query_hits
      else:
        filtered_hits += subject_hits
      #fi
    #efor

    # We may have a repeat which can confuse our results later.
    # This manifests itself as TWO overlapping hits in the blast results.
    # my solution at this moment is to take the hit which covers the largest area.

    hits = filtered_hits
    hits_sorted_by_length = sorted(hits, key=lambda x: x.qend - x.qstart, reverse=True)
    nr_hits = []
    for h_i in hits_sorted_by_length:
      overlapping = [ h_j for h_j in nr_hits if utils.blast_hits_overlap(h_i, h_j) ]
      if len(overlapping) == 0:
        nr_hits.append(h_i)
      else:
        h_j = overlapping[0]
        print("%s(%d-%d) was overlapping with %s(%d-%d)!" % (h_i.sseqid, h_i.sstart, h_i.send, h_j.sseqid, h_j.sstart, h_j.send))
      #fi
    #efor

    utils.writeBlastFile(nr_hits, output.hits)

rule all_filter_hits:
  input:
    filtered = expand("%s/{genome}/filtered_hits.tsv" % __ANALYSIS_OUTDIR__, genome=GENOMES)    

rule reconstruct_queries:
  input:
    hits = lambda wildcards: "%s/%s/filtered_hits.tsv" % (__ANALYSIS_OUTDIR__, wildcards.genome)
  output:
    reconstructed = "%s/{genome}/reconstructed_list.tsv" % __ANALYSIS_OUTDIR__
  params:
    rule_outdir = "%s/{genome}/reconstructed" % __ANALYSIS_OUTDIR__,
    blastfields = __BLASTFIELDS__
  run:
    hits    = utils.readBlastFile(input.hits, params.blastfields)

    # LOOP OVER HITS SORTED BY LENGTH!
    # Add the next hits to the set, if they are not overlapping
    # If they do overlap, then do not add to set
    # Finally, add the group to a list
    # Start again with the remaining hits
    all_reconstructed = []
    queries = set([ h.qseqid for h in hits])

    for query in queries:
      query_reconstructed = []
      # Select only the hits for this query
      parts = [ h for h in hits if h.qseqid == query ]
      # group the hits per subject (if we select a subject, we should select all the hits it has...)
      # We already guarantee that those are not overlapping, based on a previous filtering step
      parts = [ [ h for h in parts if h.sseqid == s ] for s in set([j.sseqid for j in parts ]) ]
      # Sort these groups by the size of their hits
      parts = sorted(parts, key=lambda x: (sum([ h.length for h in x]), x[0].sseqid), reverse=True)
      remaining = []
      while True:
        print("#####################")
        print("Reconstruction %d(%d)" % (len(all_reconstructed), len(parts)))
        for H in parts:
          if not(utils.blast_hits_overlap_group(query_reconstructed + H)):
            print("Successfully added to reconstruction:")
            print(H)
            query_reconstructed.extend(H)
          else:
            print("This hitset could not be added to current reconstruction:")
            print(H)
            remaining.append(H)
          #fi
        #efor
        parts = [ h for h in remaining ]
        remaining = []

        print("Final reconstruction!")
        for h in sorted(query_reconstructed, key=lambda x: x.qstart):
          print(h)
        #efor
        print("##################")
        all_reconstructed.append(query_reconstructed)
        query_reconstructed = []

        if len(parts) == 0:
          break
        #fi
      #ewhile
    #efor

    utils.mkdir_p(params.rule_outdir)
    with open(output.reconstructed, "w") as ofd_list:
      for (i,r) in enumerate(all_reconstructed):
        utils.writeBlastFile(r, "%s/reconstructed_%d.tsv" % (params.rule_outdir, i+1))
        ofd_list.write("%s\t%s\t%d\t%s\t%s\n" % ("%s.%s.%d"% (wildcards.genome, r[0].qseqid, i+1), wildcards.genome, i+1, r[0].qseqid, "%s/reconstructed_%d.tsv" % (params.rule_outdir, i+1)))
      #efor
    #ewith

rule all_reconstruct_queries:
  input:
    r = expand("%s/{genome}/reconstructed_list.tsv" % __ANALYSIS_OUTDIR__, genome=GENOMES)

rule reconstructed_masked:
  input:
    reconstructed = lambda wildcards: "%s/%s/reconstructed_list.tsv" % (__ANALYSIS_OUTDIR__, wildcards.genome),
    queries       = config["prot_query"]
  output:
    masked = "%s/{genome}/masked_reconstructions.fa" % __ANALYSIS_OUTDIR__
  params:
    blastfields = __BLASTFIELDS__
  run:
    queries = utils.loadFasta(input.queries)

    with open(input.reconstructed, "r") as ifd:
      with open(output.masked, "w") as ofd:
        reader = csv.reader(ifd, delimiter="\t")
        for (r_id, genome, r_num, query, blastfile) in reader:
          hits = utils.readBlastFile(blastfile, params.blastfields)
          covered = utils.region_covered([ (h.qstart, h.qend) for h in hits], len(queries[query]))
          masked_seq = ''.join([ 'N' if m == 0 else b for (b,m) in  zip(queries[query], covered) ])
          ofd.write(">%s\n%s\n" % (r_id, masked_seq))
        #efor
      #ewith
    #ewith
    
rule reconstructed_masked_combined:
  input:
    masked = expand("%s/{genome}/masked_reconstructions.fa" % __ANALYSIS_OUTDIR__, genome=GENOMES)
  output:
    masked = "%s/masked_reconstructions.fa" % __ANALYSIS_OUTDIR__
  shell:"""
    cat {input.masked} > {output.masked}
  """

rule reconstructions_types_blast:
  input:
    masked_seqs = rules.reconstructed_masked_combined.output.masked,
    fas12       = config["fas12"]
  output:
    hits = "%s/types/hits.tsv" % __ANALYSIS_OUTDIR__
  threads: 20
  params:
    blastfields = __BLASTFIELDS__
  shell: """
    makeblastdb -dbtype prot -in {input.masked_seqs} -out {input.masked_seqs}.db
    blastp -query {input.fas12} -db {input.masked_seqs}.db -outfmt "6 {params.blastfields}" -out {output.hits} -num_threads {threads}
  """

rule reconstructed_summary:
  input:
    reconstructed = lambda wildcards: "%s/%s/reconstructed_list.tsv" % (__ANALYSIS_OUTDIR__, wildcards.genome),
    recon_types   = rules.reconstructions_types_blast.output.hits,
    queries = config["prot_query"]
  output:
    summary = "%s/{genome}/reconstruction_summary" % __ANALYSIS_OUTDIR__
  params:
    blastfields = __BLASTFIELDS__
  run: 
    gffs    = {}
    queries = utils.loadFasta(input.queries)
    type_hits = utils.indexListBy(utils.readBlastFile(input.recon_types, params.blastfields), lambda x: x.sseqid)

    with open(input.reconstructed, "r") as ifd:
      with open(output.summary, "w") as ofd:
        ofd.write("#reconstruction_id\tgenome\treconstruction_number\tquery\treconstruction_file\tcoverage\tnpieces\tngenes\ttypes\ttandem_flag\n")
        reader = csv.reader(ifd, delimiter="\t")
        for (r_id, genome, r_num, query, blastfile) in reader:
          if genome not in gffs:
            gffs[genome] = utils.readGFF3File(dconfig["dataprefix"] + '/' + dconfig["data"][genome]["gff"])
          hits = utils.readBlastFile(blastfile, params.blastfields)

          typeQueries = set([ (h.qseqid, h.qlen) for h in type_hits.get(r_id, [])])
          typeQueriesCovered = [utils.region_covered([(h.qstart, h.qend) for h in type_hits.get(r_id, []) if h.qseqid == tqid], tqlen).sum() for (tqid,tqlen) in typeQueries]
          covered = sum(typeQueriesCovered) / float(sum([tqlen for (tqid,tqlen) in typeQueries])) * 100

          pieces  = len(hits)
          genes   = list(set([h.sseqid.split("|")[1] for h in hits]))
          ngenes  = len(genes)
          types   = sorted(list(set([ h.qseqid for h in type_hits.get(r_id, []) ])))
          types   = ','.join(types)
          neighboring = any([ (gffs[genome].areTandem(h_i, h_j) and gffs[genome].areSameStrand(h_i, h_j)) for (i, h_i) in enumerate(genes[:-1]) for (j, h_j) in enumerate(genes[i+1:]) ]) if ngenes > 1 else False
          tandemFlag  = "X" if neighboring else ""
          ofd.write("%s\t%s\t%s\t%s\t%s\t%f\t%d\t%d\t%s\t%s\n" % (r_id, genome, r_num, query, blastfile, covered, pieces, ngenes, types, tandemFlag))
        #efor
      #ewith
    #ewith

rule reconstructed_summaries:
  input:
    summary = expand("%s/{genome}/reconstruction_summary" % __ANALYSIS_OUTDIR__, genome=GENOMES)
  output:
    summary = "%s/reconstruction_summary" % __ANALYSIS_OUTDIR__
  shell: """
    cat {input.summary} | sort -r | uniq  > {output.summary}
  """

rule coverage_summary:
  input:
    reconstructed = rules.reconstructed_summaries.output.summary
  output:
    coverage_summary = "%s/coverage_summary" % __ANALYSIS_OUTDIR__
  run:
    R = { genome: { "FAS1": 0, "FAS2": 0, "FAS1,FAS2" : 0 } for genome in GENOMES}
    with open(input.reconstructed, "r") as ifd:
      reader = csv.reader(ifd, delimiter="\t")
      for row in reader:
        if row[0][0] == '#':
          continue
        #fi
        (r_id, genome, r_num, query, blastfile, covered, pieces, ngenes, types, tandemFlag) = row
        R[genome][types] += 1
      #efor
    #ewith

    with open(output.coverage_summary, "w") as ofd:
      ofd.write("%s\t%s\t%s\t%s\n" % ("genome", "fas1", "fas2", "fas1,fas2"))
      for genome in GENOMES:
        ofd.write("%s\t%d\t%d\t%d\n" % (genome, R[genome]["FAS1"], R[genome]["FAS2"], R[genome]["FAS1,FAS2"]))
      #efor
    #ewith

rule itol_summaries:
  input:
    reconstructed = rules.reconstructed_summaries.output.summary,
    template      = "%s/raw_itol_colored_gradients.txt" % __INSTALL_DIR__
  output:
    itolMin = "%s/itol_min.txt" % __ITOL_OUTDIR__,
    itolMax = "%s/itol_max.txt" % __ITOL_OUTDIR__,
    itolCopies = "%s/itol_copies.txt" % __ITOL_OUTDIR__
  run:
    R = utils.readColumnFile(input.reconstructed, "reconstruction_id genome reconstruction_number query reconstruction_file coverage npieces ngenes types tandem_flag", skip=1, types="str str int str str float int int str str")
    
    RG = utils.indexListBy([r for r in R if (r.types == "FAS1,FAS2") and (r.coverage > config["reconstruction_min_coverage"]) ], lambda x: x.genome)
    maxFASnGenesPerGenome = [ (genome, max(RG[genome], key=lambda x: x.ngenes).ngenes if genome in RG else 0) for genome in RG ]
    minFASnGenesPerGenome = [ (genome, min(RG[genome], key=lambda x: x.ngenes).ngenes if genome in RG else 0) for genome in RG ]
    
    with open(input.template, "r") as ifd:
      with open(output.itolMax, "w") as ofd:
        ofd.write(ifd.read() % ( "MaxFasPieces", "\n".join([ "%s %d" % x for x in maxFASnGenesPerGenome ])))
      #ewith
    #ewith

    with open(input.template, "r") as ifd:
      with open(output.itolMin, "w") as ofd:
        ofd.write(ifd.read() % ( "MinFasPieces", "\n".join([ "%s %d" % x for x in minFASnGenesPerGenome ])))
      #ewith
    #ewith

    nFAScopies = [ (genome, len(RG[genome])) for genome in RG ]
    with open(input.template, "r") as ifd:
      with open(output.itolCopies, "w") as ofd:
        ofd.write(ifd.read() % ( "FasCopies", "\n".join([ "%s %d" % x for x in nFAScopies ])))
      #ewith
    #ewith


rule all:
  input:
    reconstructions = rules.reconstructed_summaries.output.summary,
    coverages       = rules.coverage_summary.output.coverage_summary,
    itol            = rules.itol_summaries.output
