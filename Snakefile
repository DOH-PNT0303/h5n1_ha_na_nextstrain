#include: "rules/common.smk"

#SUPPORTED_LOCAL_SOURCES = ["ncbi", "andersen-lab", "joined-ncbi"]

#if LOCAL_INGEST:
#    assert INGEST_SOURCE in SUPPORTED_LOCAL_SOURCES, \
#        f"Full genome build is only set up for locat ingest from {SUPPORTED_LOCAL_SOURCES}."
#else:
#    assert S3_SRC.startswith("s3://nextstrain-data/"), \
#        "Full genome build is only set up for data from the public S3 bucket"

#import json

# -------------------- notes --------------------
# The approach here is to align each segment and join them (and their annotations) into a genome
# We don't join the metadata - we assume (!) all the available data is in the metadata for the HA segment
# For rules from tree onwards there is a lot of duplication between this snakefile and the
# per-segment snakefile. A config YAML would help abstract some of this out.
# -----------------------------------------------

# Segment order determines how the full genome annotation (entropy panel) is set up
# using the canonical ordering <https://viralzone.expasy.org/6>
SEGMENTS = ["ha", "na"]
assert len(set(SEGMENTS))==len(SEGMENTS), "Duplicate segment detected - check 'SEGMENTS' list"

BUILD_NAME = ['h5n1-oregon-washington']

# We parameterise the build by build_name, but we often refer to upstream files / sources by the subtype
def subtype(build_name):
    assert build_name=='h5n1-oregon-washington', "HA-NA genome build for 'h5n1-oregon-washington' "
    return 'h5n1'

def subtypes_by_subtype_wildcard(wildcards):
    db = {
        'h5nx': ['h5n1', 'h5n2', 'h5n3', 'h5n4', 'h5n5', 'h5n6', 'h5n7', 'h5n8', 'h5n9'],
        'h5n1': ['h5n1'],
        'h7n9': ['h7n9'],
        'h9n2': ['h9n2'],
    }
    return(db[wildcards.subtype])

rule all:
    input: expand("auspice/avian-flu_{build_name}_genome.json", build_name=BUILD_NAME)

# This must be after the `all` rule above since it depends on its inputs
#include: "rules/deploy.smk"

rule files:
    params:
        reference = lambda w: f"config/reference_{subtype(w.build_name)}_{{segment}}.gb",
        sequences = "outbreak_files_OR_WA/{segment}_sequences.fasta",
        metadata = "outbreak_files_OR_WA/metadata_segment_HA_NA.tsv",
        include = "outbreak_files_OR_WA/include_strains.txt",
        #dropped_strains = "config/dropped_strains_{build_name}.txt",
        colors = "config/colors_h5n1-franklin-county-outbreak.tsv",
        #lat_longs =  lambda w: f"config/lat_longs_{subtype(w.build_name)}.tsv",
        auspice_config = "config/auspice_config_h5n1-franklin-county-outbreak.json",
        #description = "config/description_{build_name}.md"

files = rules.files.params

rule filter:
    input:
        sequences = files.sequences,
        metadata = files.metadata,
        include = files.include
        #exclude = files.dropped_strains
    #params:
        #min_date = "2024-01-01",
        #query = 'region == "North America"'
    output:
        sequences = "results/{build_name}/genome/sequences_{segment}.fasta"
    log: "logs/{build_name}/genome/sequences_{segment}.txt"
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --output-log {log} \
            --output-sequences {output.sequences} \
            --include {input.include} \
        """

rule align:
    input:
        sequences ="results/{build_name}/genome/sequences_{segment}.fasta",
        reference = files.reference,
    output:
        alignment = "results/{build_name}/genome/aligned_{segment}.fasta",
    threads:
        8
    shell:
        """
        augur align \
            --sequences {input.sequences} \
            --reference-sequence {input.reference} \
            --output {output.alignment} \
            --remove-reference \
            --fill-gaps \
            --nthreads {threads}
        """

rule join_sequences:
    input:
        alignment = expand("results/{{build_name}}/genome/aligned_{segment}.fasta", segment=SEGMENTS),
    output:
        alignment = "results/{build_name}/genome/aligned.fasta",
    shell:
        """
        python scripts/join-segments.py \
            --segments {input.alignment} \
            --output {output.alignment}
        """

rule join_genbank:
    input:
        genbank_files = lambda w: [f"config/reference_{subtype(w.build_name)}_{segment}.gb" for segment in SEGMENTS],
    output:
        genbank = "results/{build_name}/genome/reference.gb",
    shell:
        """
        python scripts/join-genbank.py \
            --genbank {input.genbank_files} \
            --output {output.genbank}
        """

rule tree:
    message: "Building tree"
    input:
        alignment = "results/{build_name}/genome/aligned.fasta",
    output:
        tree = "results/{build_name}/genome/tree-raw.nwk",
    params:
        method = "iqtree"
    threads:
        8
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --method {params.method} \
            --nthreads {threads}
        """

def clock_rate(w):
    assert subtype(w.build_name)=='h5n1', 'Clock rates only available for H5N1'
    # These parameters taken from the main Snakefile
    clock_rates_h5n1 = {
        'pb2': 0.00287,
        'pb1': 0.00264,
        'pa': 0.00248,
        'ha': 0.00455,
        'np': 0.00252,
        'na': 0.00349,
        'mp': 0.00191,
        'ns': 0.00249
    }
    lengths = {
        'pb2': 2341,
        'pb1': 2341,
        'pa': 2233,
        'ha': 1760,
        'np': 1565,
        'na': 1458,
        'mp': 1027,
        'ns': 865
    }
    mean = sum([cr * lengths[seg] for seg,cr in clock_rates_h5n1.items()])/sum(lengths.values())
    stdev = mean/2
    return f"--clock-rate {mean} --clock-std-dev {stdev}"


rule refine:
    message:
        """
        Refining tree
          - estimate timetree
          - use {params.coalescent} coalescent timescale
          - estimate {params.date_inference} node dates
        """
    input:
        tree = "results/{build_name}/genome/tree-raw.nwk",
        alignment = "results/{build_name}/genome/aligned.fasta",
        metadata = files.metadata,
    output:
        tree = "results/{build_name}/genome/tree.nwk",
        node_data = "results/{build_name}/genome/branch-lengths.json"
    params:
        coalescent = "const",
        date_inference = "marginal",
        clock_rate = clock_rate,
        # Using the closest outgroup as the root
        # Make sure this strain is force included via augur filter --include
        root_strain = "A/junglecrow/Iwate/0304I001/2022"
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --metadata-id-columns 'strain'\
            --output-tree {output.tree} \
            --output-node-data {output.node_data} \
            --timetree \
            --keep-polytomies \
            --coalescent {params.coalescent} \
            --root {params.root_strain} \
            --date-confidence \
            --date-inference {params.date_inference} \
            {params.clock_rate} 
        """

rule ancestral:
    input:
        tree = "results/{build_name}/genome/tree.nwk",
        alignment = "results/{build_name}/genome/aligned.fasta",
        root_sequence = "results/{build_name}/genome/reference.gb",
    output:
        node_data = "results/{build_name}/genome/nt-muts.json"
    params:
        inference = "joint"
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --output-node-data {output.node_data} \
            --inference {params.inference} \
            --root-sequence {input.root_sequence} \
            --keep-ambiguous
        """

rule translate:
    message: "Translating amino acid sequences"
    input:
        tree = "results/{build_name}/genome/tree.nwk",
        node_data = "results/{build_name}/genome/nt-muts.json",
        reference = "results/{build_name}/genome/reference.gb",
    output:
        node_data = "results/{build_name}/genome/aa-muts.json"
    shell:
        """
        augur translate \
            --tree {input.tree} \
            --ancestral-sequences {input.node_data} \
            --reference-sequence {input.reference} \
            --output {output.node_data}
        """


rule traits:
    message: "Inferring ancestral traits for {params.columns!s}"
    input:
        tree = rules.refine.output.tree,
        metadata = files.metadata
    output:
        node_data = "results/{build_name}/genome/traits.json"
    params:
        columns = 'host',
        sampling_bias_correction = 5
    shell:
        """
        augur traits \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --output {output.node_data} \
            --columns {params.columns} \
            --sampling-bias-correction {params.sampling_bias_correction} \
            --confidence
        """

rule cleavage_site:
    input:
        ha_alignment = "results/{build_name}/genome/aligned_ha.fasta",
    output:
        cleavage_site_annotations = "results/{build_name}/genome/cleavage-site_ha.json",
        cleavage_site_sequences = "results/{build_name}/genome/cleavage-site-sequences_ha.json"
    shell:
        """
        python scripts/annotate-ha-cleavage-site.py \
            --alignment {input.ha_alignment} \
            --furin_site_motif {output.cleavage_site_annotations} \
            --cleavage_site_sequence {output.cleavage_site_sequences}
        """

rule export:
    input:
        tree = "results/{build_name}/genome/tree.nwk",
        metadata = files.metadata,
        node_data = [
            rules.refine.output.node_data,
            rules.ancestral.output.node_data,
            rules.translate.output.node_data,
            rules.traits.output.node_data,
            rules.cleavage_site.output.cleavage_site_annotations,
            rules.cleavage_site.output.cleavage_site_sequences
        ],
        colors = files.colors,
        #lat_longs = files.lat_longs,
        auspice_config = files.auspice_config,
        #description = files.description,
    output:
        auspice_json = "auspice/avian-flu_{build_name}_genome.json"
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --node-data {input.node_data} \
            --colors {input.colors} \
            --auspice-config {input.auspice_config} \
            --include-root-sequence-inline \
            --output {output.auspice_json}
        """
