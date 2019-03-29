
/*
Dependencies 
miranda 
targetScan
java 
R 
vendiagram



*/

version="0.0.1"
// Nextflow Version check
println LikeletUtils.print_green("This workflow requires Nextflow version 0.26 or greater -- You are running version ")+ LikeletUtils.print_red(nextflow.version)



params.help = null
if (params.help) {
    this.helpMessage()
    exit 0
}




inputfasta=file(params.fasta)
microRNA_fasta_db=file(params.microRNA_fasta_db)
if( !microRNA_fasta_db.exists() ) exit 1, LikeletUtils.print_red("microRNA_fasta_db  not found: ${params.microRNA_fasta_db}")

targetscan_lib=file(params.targetScanMiRNALib)
if( !targetscan_lib.exists() ) exit 1, LikeletUtils.print_red("targetScanMiRNALib  not found: ${params.targetScanMiRNALib}")



process predictByMiranda {
   publishDir path:"Result", pattern: "miranda_result.tsv",  mode: "move", overwrite: true

  input:
        file inputfasta
        file microRNA_fasta_db
  output:
        file "miranda_result.tsv" 
        file "miranda.tsv" into miranda_overlap

  shell:
  '''
    miranda  !{microRNA_fasta_db} !{inputfasta} > miranda.result.txt

    grep -A 1 "Scores for this hit:" miranda.result.txt | sort -k 2 | grep '>' |sed 's/>//' | sed '1i\\mirna Target  Score Energy-Kcal/Mol Query-Aln(start-end) Subjetct-Al(Start-End) Al-Len Subject-Identity Query-Identity' > miranda_result.tsv
    
    grep -A 1 "Scores for this hit:" miranda.result.txt | sort -k 2 | grep '>' |sed 's/>//'  | awk '{print $2"\t"$1}' - > miranda.tsv
    rm miranda.result.txt

'''
  
}

process predictByTargetScan {
  publishDir path:"Result", pattern:"RNA22_predict_result.tsv",  mode: "move", overwrite: true

  input:
        file inputfasta
        file targetscan_lib
  output:
        file "targetScan_predict_result.tsv"
        file "targetScan.tsv" into targetScan_overlap
  script:
  """
  sh ${baseDir}/bin/process_fasta_to_tab.sh ${inputfasta} > process.fa
  targetscan_70.pl ${targetscan_lib} process.fa  targetScan_predict_result.tsv
  rm process.fa 
  cut -f 1,2 targetScan_predict_result.tsv |  sed '1d' -  > targetScan.tsv  
  """
}

// RNA22 need internet connect 
process predictByRNA22 {
  publishDir path:"Result", pattern:"RNA22_predict_result.tsv", mode: "move", overwrite: true

  input:
        file inputfasta
        file microRNA_fasta_db
  output:
        file "RNA22_predict_result.tsv"
        file "RNA22.tsv" into RNA22_overlap

  script:
  """
  ln -s ${baseDir}/bin/RNA22v2.class .
  sed 's/myMirInputFile.txt/${microRNA_fasta_db}/g' ${baseDir}/bin/Parameters.properties > Parameters.properties
  sed -i 's/myTranscriptInputFile.txt/${inputfasta}/g' Parameters.properties
  sed -i 's/output.txt/RNA22_predict_result.tsv/g' Parameters.properties 
  java RNA22v2
   awk '{print $2"\t"$1}'  RNA22_predict_result.tsv  > RNA22.tsv
  """
}



process get_overlap_result {
    publishDir path:"Result/overlap", mode: "move", overwrite: true
  input:
        file miranda from miranda_overlap
        file  targetscan from targetScan_overlap
        file rna22 from RNA22_overlap
  output:
        file "*.pdf"
        file "*.txt"
  script:
  """
        Rscript ${baseDir}/bin/Venn_plot.R
      

  """
}




// Log information

/*
Working completed message
 */
workflow.onComplete {
    println LikeletUtils.print_green("=================================================")
    println LikeletUtils.print_green("Cheers! miRNA prediction task run Complete!")
    println LikeletUtils.print_green("=================================================")
    //email information
    if (params.mail) {
        recipient = params.mail
        def subject = 'My miRNAprediction-SYSUCC execution'

       def  msg = """\
            miRNAprediction-SYSUCC execution summary
            ---------------------------
            Your command line: ${workflow.commandLine}
            Completed at: ${workflow.complete}
            Duration    : ${workflow.duration}
            Success     : ${workflow.success}
            workDir     : ${workflow.workDir}
            exit status : ${workflow.exitStatus}
            Error report: ${workflow.errorReport ?: '-'}
        
            """.stripIndent()

        sendMail(to: recipient, subject: subject, body: msg)
    }
}
workflow.onError {

    println LikeletUtils.print_yellow("Oops... Pipeline execution stopped with the following message: ")+LikeletUtils.print_red(workflow.errorMessage)
}







def minimalInformationMessage() {

    println LikeletUtils.print_yellow("=====================================")
    println "\n"
    // Minimal information message
    println LikeletUtils.print_green("-------------------------------------------------------------")
    println LikeletUtils.print_green("                       Checking Parameters                   ")
    println LikeletUtils.print_green("-------------------------------------------------------------")
    checkAnalysis("\tFastq file extension:           ",params.fasta)
    println LikeletUtils.print_green("-------------------------------------------------------------")
}


def helpMessage(){

    println ''
    println LikeletUtils.print_purple('------------------------------------------------------------------------')
    println "miRNAanalysis_SYSUCC:  v$version"
    println LikeletUtils.print_purple('------------------------------------------------------------------------')
    println ''
    println LikeletUtils.print_yellow('Usage: ')
    println LikeletUtils.print_yellow('    The typical command for running the pipeline is as follows (we do not recommend users passing configuration parameters through command line, please modify the config.file instead):\n') 
            LikeletUtils.print_purple('       Nextflow run miRNAprediction/main.nf ') +

            LikeletUtils.print_yellow('    General arguments:             Input and output setting') 
            print_parmeter('--fasta <path> ','Path to input data, current path default') 
           



    println '------------------------------------------------------------------------'
    println LikeletUtils.print_yellow('Contact information: zhaoqi@sysucc.org.cn')
    println LikeletUtils.print_yellow('Copyright (c) 2013-2019 Sun Yat-sen University Cancer Center.')
    println '------------------------------------------------------------------------'

}

def print_parameter(content, parameter){
    println LikeletUtils.print_cyan(LikeletUtils,addstringToalign(content, 30))+LikeletUtils.print_green(parameter)
}

def checkAnalysis(software,param){
if(param) println LikeletUtils.print_yellow(software)+LikeletUtils.print_green(param)
}