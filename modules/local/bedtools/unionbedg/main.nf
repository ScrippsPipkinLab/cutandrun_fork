process BEDTOOLS_UNIONBEDG_MEAN {
    tag "$meta.id"
    label 'process_medium'

    conda "bioconda::bedtools=2.31.0"
    container "${ workflow.containerEngine == 'singularity' && !task.ext.singularity_pull_docker_container ?
        'https://depot.galaxyproject.org/singularity/bedtools:2.31.0--hf5e1c6e_2' :
        'biocontainers/bedtools:2.31.0--hf5e1c6e_2' }"

    input:
    tuple val(meta), path(bedgraphs)

    output:
    tuple val(meta), path("*.bedGraph"), emit: bedgraph
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def prefix  = task.ext.prefix ?: "${meta.id}"
    def bg_list = bedgraphs instanceof List ? bedgraphs : [bedgraphs]
    def n       = bg_list.size()
    def files   = bg_list.join(' ')
    """
    # Lay the already spike-in-scaled replicate bedGraphs side by side over a
    # common interval set (missing coverage filled with 0) and take the mean of
    # the per-replicate signal columns. Averaging (not summing) keeps the pooled
    # track in the same dynamic range as a single replicate and preserves the
    # per-cell spike-in scale.
    bedtools unionbedg -i $files \\
        | awk -v n=$n 'BEGIN{ OFS="\\t" } { sum = 0; for (i = 4; i <= 3 + n; i++) sum += \$i; print \$1, \$2, \$3, sum / n }' \\
        | sort -T '.' -k1,1 -k2,2n \\
        > ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.id}"
    """
    touch ${prefix}.bedGraph

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bedtools: \$(bedtools --version | sed -e "s/bedtools v//g")
    END_VERSIONS
    """
}
