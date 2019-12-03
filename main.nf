#!/usr/bin/env nextflow

/*
 * Copyright (c) 2017-2018, Centre for Genomic Regulation (CRG) and the authors.
 *
 *   This file is part of 'XXXXXX'.
 *
 *   XXXXXX is free software: you can redistribute it and/or modify
 *   it under the terms of the GNU General Public License as published by
 *   the Free Software Foundation, either version 3 of the License, or
 *   (at your option) any later version.
 *
 *   XXXXXX is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *   GNU General Public License for more details.
 *
 *   You should have received a copy of the GNU General Public License
 *   along with XXXXXX.  If not, see <http://www.gnu.org/licenses/>.
 */

/* 
 * Main XXX pipeline script
 *
 * @authors
 * Edgar Garriga

 */

/*
 * defaults parameter definitions
 */

// input sequences to align in fasta format
params.seqs = "$baseDir/test/small/ANK"

// input reference sequences aligned in 
params.refs = "$baseDir/test/ref/ANK"

// input guide trees in Newick format. Or `false` to generate trees
params.trees = false

// which alignment methods to run
params.align_method = "CLUSTALO" //,MAFFT-FFTNS1"  //,MAFFT-FFTNS1,MAFFT-GINSI,PROBCONS,UPP"

// which tree methods to run if `trees` == `false`
params.tree_method = "CLUSTALO" //,MAFFT_PARTTREE"  //,MAFFT-FFTNS1,MAFFT_PARTTREE"

// generate regressive alignments ?
params.regressive_align = true

// create standard alignments ?
params.standard_align = false

// create default alignments ? 
params.default_align = false

// evaluate alignments ?
params.evaluate = true

// bucket sizes for regressive algorithm
params.buckets= '1000'

// output directory
params.output = "$baseDir/results"


log.info """\
         PIPELINE  ~  version 0.1"
         ======================================="
         Input sequences (FASTA)                        : ${params.seqs}
         Input references (Aligned FASTA)               : ${params.refs}
         Input trees (NEWICK)                           : ${params.trees}
         Output directory (DIRECTORY)                   : ${params.output}
         Alignment methods                              : ${params.align_method}
         Tree methods                                   : ${params.tree_method}
         Generate default alignments                    : ${params.default_align}
         Generate standard alignments                   : ${params.standard_align}
         Generate regressive alignments (DPA)           : ${params.regressive_align}
         Bucket Sizes for regressive alignments         : ${params.buckets}
         Perform evaluation? Requires reference         : ${params.evaluate}
         Output directory (DIRECTORY)                   : ${params.output}
         """
         .stripIndent()


// Channels containing sequences
if ( params.seqs ) {
  Channel
  .fromPath(params.seqs)
  .map { item -> [ item.baseName, item] }
  .view()
  .into { seqsCh; seqs2 }
}

if ( params.refs ) {
  Channel
  .fromPath(params.refs)
  .map { item -> [ item.baseName, item] }
  .set { refs }
}
// Channels for user provided trees or empty channel if trees are to be generated [OPTIONAL]
if ( params.trees ) {
  Channel
    .fromPath(params.trees)
    .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item] }
    .set { trees }
}
else { 
  Channel
    .empty()
    .set { trees }
}

tree_methods = params.tree_method
align_methods = params.align_method

process reformatSeqs {
    tag "${id}"
    publishDir "${params.output}/seqs", mode: 'copy', overwrite: true

    input:
      set val(id), file(seqs) from seqsCh

    output:
     set val(id), file("${id}.fa") into reformatSeqsOut

    script:
    // sed -r '/^\s*$/d'
    """
    sed -r "/^\\s*\$/d" ${seqs} > ${id}.fa 
    """
}

process guide_trees {

    tag "${id}.${tree_method}"
    publishDir "${params.output}/guide_trees", mode: 'copy', overwrite: true
   
    input:

     set val(id), \
         file(seqs) \
         from reformatSeqsOut
     each tree_method from tree_methods.tokenize(',') 

   output:
     set val(id), \
       val(tree_method), \
       file("${id}.${tree_method}.dnd") \
       into treesGenerated

   when:
     !params.trees

   script:
     template "tree/generate_tree_${tree_method}.sh"
}
treesGenerated
  .mix ( trees )
  .combine ( seqs2, by:0 )
  .into { seqsAndTreesForStandardAlignment; seqsAndTreesForRegressiveAlignment }
process regressive_alignment {

    tag "${id}.${align_method}.DPA.${bucket_size}.${tree_method}"
    publishDir "${params.output}/alignments", mode: 'copy', overwrite: true

    input:
      set val(id), \
        val(tree_method), \
        file(guide_tree), \
        file(seqs) \
        from seqsAndTreesForRegressiveAlignment

      each bucket_size from params.buckets.tokenize(',')
       
      each align_method from align_methods.tokenize(',')   

    output:
      set val(id), \
        val("${align_method}"), \
        val(tree_method), \
        val("dpa_align"), \
        val(bucket_size), \
        file("*.aln") \
        into regressive_alignments

    when:
      params.regressive_align

    script:
       template "dpa_align/dpa_align_${align_method}.sh"
}

refs
  .cross ( regressive_alignments )
  .map { it -> [it[0][0], it[1][1], it[1][2], it[1][3], it[1][4], it[1][5], it[0][1]] }
  .view()
  .set { toEvaluate }

process eval {
    
    tag "${id}.${align_method}.${tree_method}.${align_type}.${bucket_size}"
    publishDir "${params.output}/individual_scores", mode: 'copy', overwrite: true

    input:
      set val(id), \
          val(align_method), \
          val(tree_method), \
          val(align_type), \
          val(bucket_size), \
          file(test_alignment), \
          file(ref_alignment) \
          from toEvaluate

    output:
      set val(id), \
          val(tree_method), \
          val(align_method), \
          val(align_type), \
          val(bucket_size), \
          file("*.sp"), \
          file("*.tc"), \
          file("*.col") \
          into scores

    when:
      params.evaluate

     script:
     """
       ## Sum-of-Pairs Score ##
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode sp \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.${tree_method}.sp"
       
       ## Total Column Score ##	
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode tc \
            | grep -v "seq1" | grep -v '*' | \
            awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.${tree_method}.tc"
       ## Column Score ##
       t_coffee -other_pg aln_compare \
             -al1 ${ref_alignment} \
             -al2 ${test_alignment} \
            -compare_mode column \
            | grep -v "seq1" | grep -v '*' | \
              awk '{ print \$4}' ORS="\t" \
            > "${id}.${align_type}.${bucket_size}.${align_method}.${tree_method}.col"
    """
}