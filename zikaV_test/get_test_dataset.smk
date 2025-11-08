conda: "workflow_env.yaml"

# private functions
def extract_gene_names(gff3_path):
    """
    Parses a GFF3 file and returns a list of gene names found in the attributes column.
    Parameters:
       gff3_path (str) - Path to a GFF3-format file. The file is read line by line.
    
    Behavior:
       - For each line containing the substring "gene_name", the line is split on tab characters.
       - If there are at least 9 fields, the ninth field (attributes) is processed.
       - The attributes string is split on '=' and the text after the last '=' is taken as the gene name.
       - Each extracted gene name is appended to the result list in the order encountered.
    
    Returns:
       list[str] - List of gene name strings.
    """
    gene_names = []
    with open(gff3_path, "r") as f:
        for line in f:
            if "gene_name" in line:
                fields = line.strip().split("\t")
                if len(fields) >= 9:
                    attributes = fields[8]
                    gene_name = attributes.split("=")[-1]
                    gene_names.append(gene_name)
    return gene_names

# reference and output
REFERENCE_CODE = "NC_035889.1"
OUTPUT_DIR = "dataset_test"

# get and sequences info
TAXON_ID = 64320
MIN_LENGTH = "9714" # 90% of genome length
MAX_SEQS = "1100"
ALLOWED_DIVERGENCE = "1500"

# phylonetic info
ROOTING = "mid_point"  # alternative root using outgroup, e.g. the reference "NC_035889.1"

rule all:
    input:
        "test_out"

rule get_reference:
    params:
        output_dir=OUTPUT_DIR,
        reference_code=REFERENCE_CODE,
    output:
        reference=f"{OUTPUT_DIR}/reference.fasta",
        reference_annotation=f"{OUTPUT_DIR}/genome_annotation.gff3",
        reference_gb=f"{OUTPUT_DIR}/reference.gb"
    shell:
        """
        python scripts/generate_from_genbank.py \
            --reference {params.reference_code} \
            --include-parent-features \
            --output-dir {params.output_dir}
        """

rule fetch_ncbi_dataset_package:
    output:
        dataset_package=temp("data/ncbi_dataset.zip"),
    retries: 5
    shell:
        """
        datasets download virus genome taxon {TAXON_ID} \
            --no-progressbar \
            --filename {output.dataset_package}
        """


rule extract_ncbi_dataset_sequences:
    input:
        dataset_package="data/ncbi_dataset.zip",
    output:
        ncbi_dataset_sequences="data/sequences.fasta",
    shell:
        """
        unzip -jp {input.dataset_package} \
            ncbi_dataset/data/genomic.fna \
        | seqkit seq -i -w0 \
        > {output.ncbi_dataset_sequences}
        """


rule format_ncbi_dataset_report:
    input:
        dataset_package="data/ncbi_dataset.zip",
    output:
        ncbi_dataset_tsv=temp("data/metadata_raw.tsv"),
    params:
        fields_to_include="accession,isolate-collection-date,geo-location,geo-region,isolate-lineage,release-date,submitter-affiliation,submitter-names",
    shell:
        """
        dataformat tsv virus-genome \
            --package {input.dataset_package} \
            --fields {params.fields_to_include:q} \
            > {output.ncbi_dataset_tsv}
        """


rule rename_columns:
    input:
        ncbi_dataset_tsv="data/metadata_raw.tsv",
    output:
        ncbi_dataset_tsv="data/metadata.tsv",
    shell:
        """
        sed -i -e 's/"//g' {input.ncbi_dataset_tsv}
        csvtk rename -t -f Accession,"Isolate Collection date" -n strain,date \
            {input.ncbi_dataset_tsv} \
            > {output.ncbi_dataset_tsv}
        """


rule add_reference_to_include:
    """
    Create an include file for augur filter
    """
    input:
        include="resources_test/include.txt",
    output:
        "results_dataset_test/include.txt",
    shell:
        """
        cp {input.include} results_dataset_test/include.txt
        echo "{REFERENCE_CODE}" >> results_dataset_test/include.txt
        """


rule filter:
    """
    Subsample to {MAX_SEQS} sequences.
    Only take sequences longer than {MIN_LENGTH}
    """
    input:
        sequences="data/sequences.fasta",
        metadata="data/metadata.tsv",
        include="results_dataset_test/include.txt",
    output:
        filtered_sequences="results_dataset_test/filtered_sequences_raw.fasta",
        filtered_metadata="results_dataset_test/filtered_metadata_raw.tsv",
    params:
        min_length="" if MIN_LENGTH == "" else "--min-length " + MIN_LENGTH,
        max_seqs=MAX_SEQS,
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            {params.min_length} \
            --include {input.include} \
            --subsample-max-sequences {params.max_seqs} \
            --output {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata}
        """


rule align:
    input:
        sequences="results_dataset_test/filtered_sequences_raw.fasta",
        reference=rules.get_reference.output.reference,
        annotation=rules.get_reference.output.reference_annotation,
    output:
        alignment="results_dataset_test/aligned.fasta",
        tsv="results_dataset_test/nextclade.tsv",
    params:
        translation_template=lambda w: "results_dataset_test/translations/cds_{cds}.translation.fasta",
    shell:
        """
        nextclade run \
            {input.sequences} \
            --input-ref {input.reference} \
            --input-annotation {input.annotation} \
            --output-translations {params.translation_template} \
            --output-tsv {output.tsv} \
            --output-fasta {output.alignment}
        """


rule get_outliers:
    """
    Automatically identify sequences with >{ALLOWED_DIVERGENCE} substitutions
    (likely to be sequencing errors or low quality/misannotated sequences) and put them in outliers.txt
    """
    input:
        nextclade="results_dataset_test/nextclade.tsv",
    output:
        outliers="results_dataset_test/outliers.txt",
        tmp="tmp/outliers.txt",
    params:
        allowed_divergence=lambda w: ALLOWED_DIVERGENCE,
    shell:
        """
        tsv-filter -H -v --is-numeric totalSubstitutions {input.nextclade} \
        > {output.tmp}
        tsv-filter -H \
            --is-numeric totalSubstitutions \
            --gt totalSubstitutions:{params.allowed_divergence} \
            {input.nextclade} \
        | tail -n +2 >> {output.tmp}
        cat {output.tmp} \
        | tsv-select -H -f seqName \
        | tail -n +2 > {output.outliers}
        """


rule exclude:
    """
    Rule to allow for manual and automatic exclusion of sequences
    without triggering a new subsampling that could
    surface new bad sequences resulting in an infinite loop
    """
    input:
        sequences="results_dataset_test/aligned.fasta",
        metadata="data/metadata.tsv",
        exclude="resources_test/exclude.txt",
        outliers="results_dataset_test/outliers.txt",
    output:
        filtered_sequences="results_dataset_test/filtered_aligned.fasta",
        filtered_metadata="results_dataset_test/filtered_metadata.tsv",
        strains="results_dataset_test/tree_strains.txt",
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} {input.outliers} \
            --output {output.filtered_sequences} \
            --output-metadata {output.filtered_metadata} \
            --output-strains {output.strains}
        """


rule tree:
    input:
        alignment="results_dataset_test/filtered_aligned.fasta",
    params:
        tree_argument = "--seed 123"
    output:
        tree="results_dataset_test/tree_raw.nwk",
    shell:
        """
        augur tree \
            --alignment {input.alignment} \
            --output {output.tree} \
            --tree-builder-args '{params.tree_argument}'
        """


rule refine:
    input:
        tree="results_dataset_test/tree_raw.nwk",
        alignment="results_dataset_test/filtered_aligned.fasta",
    output:
        tree="results_dataset_test/tree.nwk",
        node_data="results_dataset_test/branch_lengths.json",
    shell:
        """
        augur refine \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --root {ROOTING} \
            --keep-polytomies \
            --divergence-units mutations \
            --output-node-data {output.node_data} \
            --output-tree {output.tree}
        """

#GENES = extract_gene_names(rules.get_reference.output.reference_annotation)
rule ancestral:
    input:
        tree="results_dataset_test/tree.nwk",
        alignment="results_dataset_test/filtered_aligned.fasta",
        annotation=rules.get_reference.output.reference_gb,
        gff3=rules.get_reference.output.reference_annotation,
    params:
        translation_template=r"results_dataset_test/translations/cds_%GENE.translation.fasta",
        output_translation_template=r"results_dataset_test/translations/cds_%GENE.ancestral.fasta",
        genes=lambda wc, input: " ".join(extract_gene_names(input.gff3))
    output:
        node_data="results_dataset_test/muts.json",
    shell:
        """
        augur ancestral \
            --tree {input.tree} \
            --alignment {input.alignment} \
            --annotation {input.annotation} \
            --root-sequence {input.annotation} \
            --genes {params.genes} \
            --translations {params.translation_template} \
            --output-node-data {output.node_data} \
            --output-translations {params.output_translation_template}
        """


rule dummy_clades:
    """
    Nextclade requires clade membership to be specified for each node
    in the tree. This rule creates a dummy clade membership for each node
    """
    input:
        json="results_dataset_test/branch_lengths.json",
        clades="resources_test/node_clades.tsv"  # formato: node<TAB>clade
    output:
        "results_dataset_test/dummy_clades.json",
    run:
        import json

        clade_map = {}
        with open(input.clades) as f:
            for line in f:
                if line.strip():
                    node, clade = line.strip().split("\t")
                    clade_map[node] = clade

        with open(input.json) as f:
            data = json.load(f)

        nodes = data.get("nodes", {})
        new_nodes = {
            node: {"clade_membership": clade_map.get(node, "Unknown")}
            for node in nodes
        }

        data["nodes"] = new_nodes

        with open(output[0], "w") as f:
            json.dump(data, f, indent=2)


rule export:
    input:
        tree="results_dataset_test/tree.nwk",
        metadata="results_dataset_test/filtered_metadata.tsv",
        mutations="results_dataset_test/muts.json",
        branch_lengths="results_dataset_test/branch_lengths.json",
        clades="results_dataset_test/dummy_clades.json",
        auspice_config="resources_test/auspice_config.json",
    output:
        auspice="results_dataset_test/auspice.json",
    shell:
        """
        augur export v2 \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --auspice-config {input.auspice_config} \
            --node-data {input.mutations} {input.branch_lengths} {input.clades} \
            --output {output.auspice}
        """


rule subsample_example_sequences:
    input:
        all_sequences="data/sequences.fasta",
        tree_strains="results_dataset_test/tree_strains.txt",
    output:
        example_sequences="results_dataset_test/example_sequences.fasta",
    shell:
        """
        # Exclude tree sequences from all sequences
        seqkit grep -v -f {input.tree_strains} {input.all_sequences} \
        | seqkit sample -n 100 -s 42 > {output.example_sequences}
        """


rule assemble_dataset:
    input:
        tree="results_dataset_test/auspice.json",
        reference=rules.get_reference.output.reference,
        annotation=rules.get_reference.output.reference_annotation,
        sequences="results_dataset_test/example_sequences.fasta",
        pathogen="resources_test/pathogen.json",
    output:
        tree="zika_dataset_test/tree.json",
        reference="zika_dataset_test/reference.fasta",
        annotation="zika_dataset_test/genome_annotation.gff3",
        sequences="zika_dataset_test/sequences.fasta",
        pathogen="zika_dataset_test/pathogen.json",
        readme="zika_dataset_test/README.md",
        changelog="zika_dataset_test/CHANGELOG.md",
        dataset_zip="zika_dataset_test.zip",
    shell:
        """
        cp {input.tree} {output.tree}
        cp {input.reference} {output.reference}
        cp {input.annotation} {output.annotation}
        cp {input.sequences} {output.sequences}
        cp {input.pathogen} {output.pathogen}
        touch README.md {output.readme}
        touch CHANGELOG.md {output.changelog}
        zip -rj zika_dataset_test.zip  zika_dataset_test/*
        """


rule test:
    input:
        dataset="zika_dataset_test.zip",
        sequences="zika_dataset_test/sequences.fasta",
    output:
        output=directory("test_out"),
    shell:
        """
        nextclade run \
            --input-dataset {input.dataset} \
            --output-all {output.output} \
            {input.sequences}
        """