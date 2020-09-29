#!/usr/bin/env nextflow

// input files
params.msa = "/users/cn/sjin/projects/homoplasy/slaveTree_results/results_1000_CLUSTALO/alignments/{seatoxin,hip,scorptoxin,cyt3,rnasemam,bowman,toxin,ghf11,TNF,sti,Stap_Strp_toxin,profilin,ricin,ghf22,ChtBD,ins,trfl,slectin,phoslip,ltn,il8,az,kringle,cryst,DEATH,cah,mmp,rub,ghf10,tgfb,sodcu,KAS,DMRL_synthase,tms,GEL,kunitz,Sulfotransfer,mofe,Ald_Xan_dh_2,ghf5,phc,aadh,annexin,serpin,cytb,asp,oxidored_q6,hpr,hormone_rec,hr,tim,glob,ace,cys,ghf1,sodfe,peroxidase,uce,flav,HMG_box,OTCace,msb,icd,proteasome,cyclo,LIM,HLH,ldh,subt,int,lyase_1,gpdh,egf,blm}.*.aln"
// params.msa = "/users/cn/sjin/projects/homoplasy/slaveTree_results/results_1000_CLUSTALO/alignments/seatoxin.*.aln"
params.trees = "/users/cn/sjin/projects/homoplasy/trees_nolength/*.{codnd,dpparttreednd0,fastaparttreednd,fftns1dnd,mafftdnd,parttreednd0,FAMSA}.dnd"
// params.trees = "/users/cn/sjin/projects/homoplasy/trees_nolength/*.FAMSA.dnd"
params.regtrimN = "1000"
params.regtrim = false
params.checkexist = true

// output directory
params.output = "/users/cn/sjin/projects/homoplasy/slaveTree_results/results_50_CLUSTALO/pastml"
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
  .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[5], item ] }   // family name, tree method, path to msa
  .into { ch_msa; ch_msa2 }
  // .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[4], item ] }   // family name, tree method, path to msa

// channel containing trees
Channel
  .fromPath(params.trees)
  .map { item -> [ item.baseName.tokenize('.')[0], item.baseName.tokenize('.')[1], item ] }    // family name, tree method, path to tree
  .into { ch_trees; ch_trees2 }

// merge channels
ch_msa
  .join(ch_trees, by:[0,1])
  .set { ch_input }  // family name, tree method, path to msa, path to tree



/* PRUNE MSA: KEEP N INFORMATIVE SEQUENCES */
process regtrim_msa {

    container 'suzannejin/tcoffee:slave'

    tag "${id}.${tree_method}"
    publishDir "${params.output}/${id}_${tree_method}", mode: 'copy', overwrite: true

    label 'process_low'

    input:

     set val(id), \
         val(tree_method), \
         file(msa), \
         file(tree) \
         from ch_input
      each regtrimN from params.regtrimN.tokenize(',')

    output:
     set val(id), \
         val(tree_method), \
         file("msa") \
         into ch_msapruned

    when:
      params.regtrim

    script:
     """
     t_coffee -other_pg seq_reformat -in $msa -in2 $tree -action +regtrim $regtrimN > msa    

     """
}

// merge msa with trees
if (params.regtrim) {
  ch_msapruned
    .join(ch_trees2, by:[0,1])
    .set{ ch_toPrepastml }
} else {
  ch_msa2
    .join(ch_trees2, by:[0,1])
    .set{ ch_toPrepastml }
}


/* PREPARE THE INPUT DATA FOR PASTML
*
* 1. Prune tree
* 2. Get annotation table in csv, with gap/ungapped positions signed as 0/1
*/ 
process prepastml {

    tag "${id}.${tree_method}"
    publishDir "${params.output}/${id}_${tree_method}", mode: 'copy', overwrite: true

    label 'process_low'

    input:

     set val(id), \
         val(tree_method), \
         file(msa), \
         file(tree) \
         from ch_toPrepastml
      each regtrimN from params.regtrimN.tokenize(',')

    output:
     set val(id), \
         val(tree_method), \
         file("states"), \
         file("tree"), \
         file("msa") \
         into ch_toPastml

    when:
      if (params.checkexist){!(file("${params.output}/${id}_${tree_method}/steps").exists())}else{true}

    script:
     """
     # prune tree
     python $baseDir/bin/pastml/prune_tree_from_msa.py -msa $msa -tree $tree -out tree

     # get annotation table
     python $baseDir/bin/pastml/get_states.py -msa $msa -out states

     mv $msa msa
     
     """
}


/* RUN PASTML 
*
* 1. Determine ancestral states using a max parsimony method -> ACCTRAN
* We are using a MP approach, instead of ML, since we want to save computation time and also some of our trees don't have branch length for ML.
* 2. Count total state transitions across columns.
*/
process pastml {

    tag "${id}.${tree_method}"
    publishDir "${params.output}/${id}_${tree_method}", mode: 'copy', overwrite: true

    label 'process_high'

    input:
     set val(id), \
         val(tree_method), \
         file(states), \
         file(tree), \
         file(msa) \
         from ch_toPastml
    each regtrimN from params.regtrimN.tokenize(',') 

    output:
     set val(id), \
         val(tree_method), \
         file("steps*"), \
         file("hpastml"), \
         file("out"), \
         file("*.tab"), \
         file("named.tree*") \
         into ch_output

    when:
      if (params.checkexist){!(file("${params.output}/${id}_${tree_method}/steps").exists())}else{true}

    script:
     """
     # run pastml
     pastml -t $tree -d $states -s , --prediction_method ACCTRAN --work_dir . -o out
     
     # count state transitions
     # write the counts to files: steps, steps2, steps3
     python /nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/bin/pastml/count_steps.py -dir . -msa $msa -bycol -bycolrow

     # count homoplasic steps
     python /nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty/bin/pastml/count_homoplasy.py -tree named.tree_tree -action homo > hpastml

     """
}

