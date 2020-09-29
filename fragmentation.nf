#!/usr/bin/env nextflow

// input files
params.msa = "/users/cn/sjin/projects/homoplasy/slaveTree_results/results_1000_CLUSTALO/alignments/{seatoxin,hip,scorptoxin,cyt3,rnasemam,bowman,toxin,ghf11,TNF,sti,Stap_Strp_toxin,profilin,ricin,ghf22,ChtBD,ins,trfl,slectin,phoslip,ltn}.*.aln" //,il8,az,kringle,cryst,DEATH,cah,mmp,rub,ghf10,tgfb,sodcu,KAS,DMRL_synthase,tms,GEL,kunitz,Sulfotransfer,mofe,Ald_Xan_dh_2,ghf5,phc,aadh,annexin,serpin,cytb,asp,oxidored_q6,hpr,hormone_rec,hr,tim,glob,ace,cys,ghf1,sodfe,peroxidase,uce,flav,HMG_box,OTCace,msb,icd,proteasome,cyclo,LIM,HLH,ldh,subt,int,lyase_1,gpdh,egf,blm}.*.aln"
// params.msa = "/users/cn/sjin/projects/homoplasy/slaveTree_results/results_1000_CLUSTALO/alignments/{il8,az,kringle,cryst,DEATH,cah,mmp,rub,ghf10,tgfb,sodcu,KAS,DMRL_synthase,tms,GEL,kunitz,Sulfotransfer,mofe,Ald_Xan_dh_2,ghf5,phc,aadh,annexin,serpin,cytb,asp,oxidored_q6,hpr,hormone_rec,hr,tim,glob,ace,cys,ghf1,sodfe,peroxidase,uce,flav,HMG_box,OTCace,msb,icd,proteasome,cyclo,LIM,HLH,ldh,subt,int,lyase_1,gpdh,egf,blm}.*.aln"
// params.msa = "/users/cn/sjin/projects/homoplasy/slaveTree_results/results_1000_CLUSTALO/alignments/seatoxin.*.aln"
params.trees = "/users/cn/sjin/projects/homoplasy/trees_nolength/*.{codnd,dpparttreednd0,fastaparttreednd,fftns1dnd,mafftdnd,parttreednd0,FAMSA}.dnd"
// params.trees = "/users/cn/sjin/projects/homoplasy/trees_nolength/*.FAMSA.dnd"
params.post = 'none'
params.bjorn = '' // '-bjorn'
params.checkexist = true

// output directory
params.output = "/users/cn/sjin/projects/homoplasy/slaveTree_results/results_50_CLUSTALO/fragmentation"
// params.output = "/users/cn/sjin/test_pastml"


log.info """\
         P A S T M L - version 1.9.24"
         ======================================="
         Input alignments (FASTA)                       : ${params.msa}
         Input trees (NEWICK)                           : ${params.trees}
         Output directory (DIRECTORY)                   : ${params.output}
         """
         .stripIndent()


// channel from msa
Channel
  .fromPath(params.msa)
  .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[5], item.baseName, item ] }   // family name, tree method, id, path to msa
  .set { ch_msa }

// channel containing trees
Channel
  .fromPath(params.trees)
  .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item ] }    // family name, tree method, path to tree
  .set { ch_trees }

// merge channels
ch_msa
  .join(ch_trees, by:[0,1])
  .set { ch_input }  // family name, tree method, id, path to msa, path to tree



/* DETERMINE FRAGMENTATION SCORE */
process fragmentation {

    tag "${id}"
    publishDir "${params.output}", mode: 'copy', overwrite: true

    input:

     set val(fam), \
         val(tree_method), \
         val(id), \
         file(msa), \
         file(tree) \
         from ch_input
      val(post) from params.post
      val(bjorn) from params.bjorn

    output:
     set val(id), \
         file("${id}.fragmentation") \
         into ch_output

    when:
      if (params.checkexist){!(file("${params.output}/${id}.fragmentation").exists())}else{true}

    script:
     """
     python /users/cn/sjin/projects/homoplasy/nf_homoplasty/bin/fragmentation/fragmentation.py -msa $msa -tree $tree -post $post $bjorn > ${id}.fragmentation

     """
}
