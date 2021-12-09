#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// help message
if (params.help) { exit 0, helpMSG() }

// error codes
if (params.profile) {
    exit 1, "--profile is WRONG use -profile" }
if (params.workdir) {
    exit 1, "--workdir is WRONG use -w" }

// warnings
if ( workflow.profile == 'standard' ) { 
    println " "
    println "\033[0;33mWARNING:NO EXECUTION PROFILE SELECTED, using [-profile local,conda]" }

println " "
println "\u001B[32mProfile: $workflow.profile\033[0m"
println " "
println "\033[2mCurrent User: $workflow.userName"
println "Nextflow-version: $nextflow.version"
println "Starting time: $workflow.start"
println "Workdir location:"
println "  $workflow.workDir"
println "Launchdir location:"
println "  $workflow.launchDir"
if ( workflow.profile.contains('singularity') ) {
    println "Singularity cache directory:"
    println "  $params.singularity_cache_dir"
}
if ( workflow.profile.contains('conda') || workflow.profile.contains('standard') ) { 
    println "Conda cache directory:"
    println "  $params.conda_cache_dir"
}
println "Configuration files:"
println "  $workflow.configFiles"
println "Cmd line:"
println "  $workflow.commandLine\u001B[0m"
if (workflow.repository != null){ println "\033[2mGit info: $workflow.repository - $workflow.revision [$workflow.commitId]\u001B[0m" }
println " "
if (workflow.profile.contains('standard') || workflow.profile.contains('local')) {
    println "\033[2mCPUs to use: $params.cores, maximal CPUs to use: $params.max_cores\u001B[0m"
    println " "
}

if ( !params.cache_dir ) {
    exit 1, "input missing, use [--cache_dir]"
}

// input
seq_ch = Channel.fromFilePairs( "$params.cache_dir/$params.seq_dir/*/*.[a,b]seq" ){ file -> file.getParent().getName() }
// globs, returns a key and a tuple of two flies with matching pattern
// emits [parent_dir, [xyz.aseq, xyz.bseq]]
// see https://www.nextflow.io/docs/latest/channel.html#fromfilepairs

process align {
    label 'emboss'
    publishDir "${params.cache_dir}/${params.algn_dir}/${parent}/", mode: params.publish_dir_mode
    
    input:
    tuple val(parent), path(aseq_bseq)

    output:
    tuple val(parent), path("*.algn")

    script:
    """
    stretcher -asequence ${aseq_bseq[0]} -bsequence ${aseq_bseq[1]} -outfile ${aseq_bseq[0].baseName}.algn -auto -aformat markx3
    """
    stub:
    """
    touch ${aseq_bseq[0].baseName}.algn
    """
}

process diff {
    label 'python'
    publishDir "${params.cache_dir}/${params.diff_dir}/${parent}/", mode: params.publish_dir_mode
    
    input:
    tuple val(parent), path(algn)

    output:
    path("${algn.baseName}.diff")

    script:
    """
    markx3_to_diff.py ${algn} ${algn.baseName}.diff
    """
    stub:
    """
    touch ${algn.baseName}.diff
    """
}

workflow {
    align(seq_ch)
    diff(align.out)
}

def helpMSG() {
    c_green = "\033[0;32m";
    c_reset = "\033[0m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_dim = "\033[2m";
    log.info """
    ${c_yellow}Input:${c_reset}
    --cache_dir             Location of the cache. REQUIRED [default: $params.cache_dir]

    --algn_dir              Name of the cache algn subdirectory [default: $params.algn_dir]
    --diff_dir              Name of the cache diff subdirectory [default: $params.diff_dir]
    --seq_dir               Name of the cache seq subdirectory [default: $params.seq_dir]

    ${c_yellow}Computing options:${c_reset}
    --cores                  Max cores per process for local use [default: $params.cores]
    --max_cores              Max cores used on the machine for local use [default: $params.max_cores]
    --memory                 Max memory in GB for local use [default: $params.memory]

    ${c_yellow}Caching:${c_reset}
    --conda_cache_dir          Location for storing the conda environments [default: $params.conda_cache_dir]
    --singularity_cache_dir    Location for storing the singularity images [default: $params.singularity_cache_dir]
    --publish_dir_mode       'copy' or 'symlink' the results [default: $params.publish_dir_mode]

    ${c_dim}Nextflow options:
    -with-tower              Activate monitoring via Nextflow Tower (needs TOWER_ACCESS_TOKEN set).
    -with-report rep.html    CPU / RAM usage (may cause errors).
    -with-dag chart.html     Generates a flowchart for the process tree.
    -with-timeline time.html Timeline (may cause errors).${c_reset}

    ${c_yellow}Execution/Engine profiles:${c_reset}
    The pipeline supports profiles to run via different ${c_green}Executers${c_reset} and ${c_blue}Engines${c_reset} e.g.: -profile ${c_green}local${c_reset},${c_blue}conda${c_reset}
    
    ${c_green}Executer${c_reset} (choose one):
      local
      slurm
    
    ${c_blue}Engines${c_reset} (choose one):
      conda
      docker
      singularity
    
    Per default: -profile local,conda is executed. 
    ${c_reset}
    """.stripIndent()
}