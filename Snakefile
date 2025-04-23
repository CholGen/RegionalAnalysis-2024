workdir: "analyses/"

# Inputs.
METADATA = ["../data/supplemental_data1.csv", "../data/supplemental_data2.csv"]
ML_TREE = "../data/2025.04.23_ML_GTR.tree"
BEAST_LOG = "../beast-analyses/2025-05-28_constant_relaxed.combined.log"
BEAST_MCC = "../beast-analyses/2025-05-28_constant_relaxed.mcc.tree"
BEAST_TREES = "../beast-analyses/2025-05-28_constant_relaxed.combined.down.trees"
BEAST_DISCRETE_LOG = "../beast-analyses/2025-06-24_constant_relaxed_discrete.log"
BEAST_DISCRETE_MCC = "../beast-analyses/2025-06-24_constant_relaxed_discrete.mcc.tree"
BEAST_DISCRETE_TREES = "../beast-analyses/2025-06-24_constant_relaxed_discrete.trees"
BEAST_MARKOV_JUMPS = "../beast-analyses/2025-06-24_constant_relaxed_discrete.Location.csv"
VCF = "../data/background_2025-04-23.bcf.gz"
REFERENCE = "../data/reference.fasta"
MASK = "../data/recombinant_regions.gff"
GENBANK = "../data/cholera_reference_NCBI.gb"
AUSPICE_CONFIG = "../nextstrain/auspice_config.json"
LAT_LONGS = "../nextstrain/lat_long_addendum.tsv"
COUNTRIES = ["Cameroon", "DRC", "Malawi", "Mozambique", "Nigeria", "Uganda", "Zambia"]


rule all:
    input:
        #one_a = "plots/figure1-sequences-map.pdf",
        one_b = "plots/figure1-sequences-cases-over-time.pdf",
        two_a = "plots/figure2-mcc-tree.pdf",
        two_b = "plots/figure2-lineages-per-country-over-time.pdf",
        two_c = "plots/figure2-AMR-resistance-per-year.pdf",
        three_a = "plots/figure3-median-markov-jumps-heatmap.pdf",
        three_b = "plots/figure3-markov-jumps-sequences-supported-over-time.pdf",
        three_c = "plots/figure3-sequences-vs-jumps.pdf",
        four = expand( "plots/figure4-surveillance-assessment-{country}.pdf", country=COUNTRIES ),
        s_one_pdf = "plots/figureS1-ML-tree.pdf",
        s_one_png = "plots/figureS1-ML-tree.png",
        s_two_pdf = "plots/figureS2-root-to-tip.pdf",
        s_two_png = "plots/figureS2-root-to-tip.png",
        s_three_pdf = "plots/figureS3-lineage-substitution-rates.pdf",
        s_three_png = "plots/figureS3-lineage-substitution-rates.png",
        s_four_pdf = "plots/figureS4-AMR-resistance-over-time.pdf",
        s_four_png= "plots/figureS4-AMR-resistance-over-time.png",
        s_five_pdf= "plots/figureS5-phylogenetic-tree-AMR.pdf",
        s_five_png= "plots/figureS5-phylogenetic-tree-AMR.png",
        s_one= "../figures/figureS1.pdf",
        s_two= "../figures/figureS2.pdf",
        s_three="../figures/figureS3.pdf",
        s_four= "../figures/figureS4.pdf",
        s_five= "../figures/figureS5.pdf",
        nucleotide_substitutions = "substitutions/substitutions.csv",
        nextstrain_tre = "../nextstrain/auspice.json",


rule generate_figure_one:
    input:
        #map_script = "analyses/figure1-map-sequences-CholGen.ipynb",
        metadata = METADATA,
        sequences_script = "figure1-sequences-available.ipynb"
    output:
        #one_a = "analyses/plots/figure1-sequences-map.pdf",
        one_b = "plots/figure1-sequences-cases-over-time.pdf"
    shell:
        """
        jupyter execute {input.sequences_script}
        """


rule generate_figure_two:
    input:
        metadata = METADATA,
        mcc_script = "figure2-mcc-tree.ipynb",
        mcc_tree= BEAST_MCC,
        T15_members= "../data/t15_members.txt",
        lineage_script = "figure2-lineage-counts-over-time.ipynb",
        amr_script = "figure2-amr-resistance.ipynb",
        amr_refgenes = "../data/amr_refgenes.csv"
    output:
        two_a = "plots/figure2-mcc-tree.pdf",
        two_b = "plots/figure2-lineages-per-country-over-time.pdf",
        two_c = "plots/figure2-AMR-resistance-per-year.pdf",
        s_four_pdf = "plots/figureS4-AMR-resistance-over-time.pdf",
        s_four_png= "plots/figureS4-AMR-resistance-over-time.png",
    shell:
        """
        jupyter execute {input.mcc_script} {input.lineage_script} {input.amr_script}
        """

rule combine_metadata:
    input:
        supp1 = "../data/supplemental_data1.csv",
        supp2 = "../data/supplemental_data2.csv"
    output:
        md = temp( "../data/vc_metadata_workshop.csv" )
    run:
        import pandas as pd

        md = list()
        for file in [input.supp1, input.supp2]:
            df = pd.read_csv( file, usecols=["taxa", "collection_date"] )
            df["workshop"] = (file == input.supp1)
            md.append( df )

        md = pd.concat( md )
        md.to_csv( output.md, index=False )


rule extract_jumps:
    input:
        trees = BEAST_DISCRETE_TREES,
        metadata = rules.combine_metadata.output.md
    output:
        introductions = "estimate-introductions/introductions.csv"
    shell:
        """
        python estimate-introductions/estimate-introductions.py \
            --trees {input.trees} \
            --burnin 100 \
            --states Nigeria Cameroon "Democratic Republic of the Congo" Uganda Zambia Malawi Mozambique \
            --trait Location \
            --metadata {input.metadata} \
            --output {output.introductions}
        """


rule generate_figure_three:
    input:
        markov_script = "figure3-markov-jumps.ipynb",
        markov_jumps = BEAST_MARKOV_JUMPS
    output:
        three_a = "plots/figure3-median-markov-jumps-heatmap.pdf",
        three_b = "plots/figure3-markov-jumps-sequences-supported-over-time.pdf",
        three_c = "plots/figure3-sequences-vs-jumps.pdf"
    shell:
        """
        jupyter execute {input.markov_script}
        """


rule estimate_substitutions_per_site:
    input:
        trees = BEAST_TREES,
        subs_script = "population-estimates/sub-per-site-gained.ipynb"
    output:
        subs_per_site = "population-estimates/subs-site-per-country.csv"
    shell:
        """
        jupyter execute {input.subs_script}
        """


rule rarefaction_analysis:
    input:
        introductions = rules.extract_jumps.output.introductions,
        rarefaction_script = "population-estimates/estimate-population-size.ipynb"
    output:
        rarefaction_curves = "population-estimates/rarefaction_curve.csv",
        rarefaction_model_fit = "population-estimates/rarefaction-model-fit.csv"
    shell:
        """
        jupyter execute {input.rarefaction_script}
        """


rule generate_figure_four:
    input:
        rarefaction_curves = rules.rarefaction_analysis.output.rarefaction_curves,
        rarefaction_fit = rules.rarefaction_analysis.output.rarefaction_model_fit,
        subs_per_site = rules.estimate_substitutions_per_site.output.subs_per_site,
        surveillance_script = "figure4-surveillance-assessment.ipynb"
    output:
        fours = expand( "plots/figure4-surveillance-assessment-{country}.pdf", country=COUNTRIES )
    shell:
        """
        jupyter execute {input.surveillance_script}
        """


rule generate_supp_figure_one:
    input:
        metadata = METADATA,
        tree = ML_TREE,
        ml_script = "figureS1-ML-tree.ipynb"
    output:
        s_tree_pdf = "plots/figureS1-ML-tree.pdf",
        s_tree_png = "plots/figureS1-ML-tree.png",
    shell:
        """
        jupyter execute {input.ml_script}
        """



rule generate_supp_figure_two:
    input:
        tree = ML_TREE,
        rtt_script = "figureS2-root-to-tip.ipynb"
    output:
        s_two_pdf = "plots/figureS2-root-to-tip.pdf",
        s_two_png = "plots/figureS2-root-to-tip.png"
    shell:
        """
        jupyter execute {input.rtt_script}
        """


rule extract_lineage_rates:
    input:
        tree = BEAST_MCC,
        extract_script = "lineage-rates/extract-rates.ipynb"
    output:
        rates = "lineage-rates/lineage_rates.csv"
    shell:
        """
        jupyter execute {input.extract_script}
        """


rule generate_supp_figure_three:
    input:
        rates = rules.extract_lineage_rates.output.rates,
        rates_script = "figureS3-lineage-rates.ipynb"
    output:
        s_three_pdf = "plots/figureS3-lineage-substitution-rates.pdf",
        s_three_png = "plots/figureS3-lineage-substitution-rates.png"
    shell:
        """
        jupyter execute {input.rates_script}
        """


rule generate_supp_figure_five:
    input:
        tree = BEAST_MCC,
        inc_script = "figureS5-IncAC-tree.ipynb"
    output:
        s_five_pdf = "plots/figureS5-phylogenetic-tree-AMR.pdf",
        s_five_png = "plots/figureS5-phylogenetic-tree-AMR.png"
    shell:
        """
        jupyter execute {input.inc_script}
        """

rule move_supplemental_data:
    input:
        s_one = rules.generate_supp_figure_one.output.s_tree_pdf,
        s_two = rules.generate_supp_figure_two.output.s_two_pdf,
        s_three = rules.generate_supp_figure_three.output.s_three_pdf,
        s_four = rules.generate_figure_two.output.s_four_pdf,
        s_five = rules.generate_supp_figure_five.output.s_five_pdf,
    output:
        s_one = "../figures/figureS1.pdf",
        s_two= "../figures/figureS2.pdf",
        s_three="../figures/figureS3.pdf",
        s_four= "../figures/figureS4.pdf",
        s_five= "../figures/figureS5.pdf",
    shell:
        """
        cp {input.s_one} {output.s_one}
        cp {input.s_two} {output.s_two}
        cp {input.s_three} {output.s_three}
        cp {input.s_four} {output.s_four}
        cp {input.s_five} {output.s_five}
        """

rule generate_fasta_alignment:
    input:
        vcf = VCF,
        reference = REFERENCE
    output:
        alignment = temp( "../nextstrain/intermediates/whole_alignment.fasta" )
    shell:
        """
        bcftools index -f {input.vcf} &&\
        python ../nextstrain/vcf_to_fasta.py \
            --vcf {input.vcf} \
            --reference {input.reference} \
            --output {output.alignment}
        """

rule mask_fasta_alignment:
    input:
        alignment = rules.generate_fasta_alignment.output.alignment,
        mask = MASK
    output:
        masked_alignment = temp( "../nextstrain/intermediates/whole_alignment.masked.fasta" )
    shell:
        """
        python ../nextstrain/mask_gubbins_aln.py \
            --aln {input.alignment} \
            --gff {input.mask} \
            --out {output.masked_alignment}
        """

rule add_reference_to_mask:
    input:
        alignment = rules.mask_fasta_alignment.output.masked_alignment,
        reference = REFERENCE
    output:
        ref_alignment = temp( "../nextstrain/intermediates/whole_alignment.ref.fasta" )
    shell:
        """
        cat {input.reference} {input.alignment} > {output.ref_alignment}
        """

rule generate_masked_vcf:
    input:
        alignment = rules.add_reference_to_mask.output.ref_alignment
    output:
        vcf = temp( "../nextstrain/intermediates/masked_alignment.vcf.gz" )
    shell:
        """
        snp-sites -v {input.alignment} | bcftools view -Oz -o {output.vcf}
        bcftools index {output.vcf}
        """

rule scale_phylogeny:
    input:
        tree = ML_TREE
    params:
        scale_value = 1.0 / 12192.0
    output:
        scaled_tree = temp( "../nextstrain/intermediates/scaled.tree" )
    shell:
        """
        gotree brlen scale \
            -f {params.scale_value} \
            -i {input.tree} \
            -o {output.scaled_tree}
        """


rule format_metadata:
    input:
        metadata = rules.combine_metadata.output.md,
        tree = ML_TREE
    output:
        dates=temp( "../nextstrain/intermediates/dates.csv" )
    run:
        import pandas as pd

        taxa = set( i for i in shell( f"gotree labels -i {input.tree}",iterable=True ) )

        md = pd.read_csv( input.metadata,usecols=["taxa", "collection_date"] )
        md = md.loc[md["taxa"].isin( taxa )]

        assert md.shape[0] == len( taxa )

        md = md.rename( columns={ "taxa": "name", "collection_date": "date" } )
        md.to_csv( output.dates,index=False )



rule combine_all_metadata:
    input:
        supp1 = METADATA[0],
        supp2 = METADATA[1]
    output:
        md = temp( "../nextstrain/intermediates/metadata.csv" )
    run:
        import pandas as pd

        md = list()
        for file in [input.supp1, input.supp2]:
            df = pd.read_csv( file )
            df["workshop"] = (file == input.supp1)
            md.append( df )

        md = pd.concat( md )
        md = md.rename( columns={"IncA/C" : "IncAC"} )
        md["accession"] = md["accession"].replace( {"CP001485/CP001486" : "CP001485", "AE003852/AE003853" : "AE003852" } )
        md.to_csv( output.md, index=False )


rule refine:
    input:
        tree = rules.scale_phylogeny.output.scaled_tree,
        alignment = rules.generate_masked_vcf.output.vcf,
        metadata = rules.format_metadata.output.dates,
        reference = REFERENCE,
    output:
        time_tree = temp( "../nextstrain/time.tree" ),
        branch_lengths = temp( "../nextstrain/branch_lengths.json" )
    shell:
        """
        augur refine \
            --alignment {input.alignment} \
            --tree {input.tree} \
            --metadata {input.metadata} \
            --clock-filter-iqd 4 \
            --vcf-reference {input.reference} \
            --divergence-units mutations \
            --timetree \
            --use-fft \
            --stochastic-resolve \
            --output-tree {output.time_tree} \
            --output-node-data {output.branch_lengths}
        """

rule ancestral:
    input:
        time_tree = rules.refine.output.time_tree,
        alignment = rules.generate_masked_vcf.output.vcf,
        reference = REFERENCE
    output:
        nucleotide_substitutions = "../nextstrain/nt-muts.json",
        ancestral_sequences = temp( "../nextstrain/intermediates/ancestral.vcf.gz" )
    shell:
        """
        augur ancestral \
            --tree {input.time_tree} \
            --alignment {input.alignment} \
            --vcf-reference {input.reference} \
            --output-node-data {output.nucleotide_substitutions} \
            --output-vcf {output.ancestral_sequences}
        """

rule extract_mutations:
    input:
        nucleotide_substitutions = rules.ancestral.output.nucleotide_substitutions,
        metadata = METADATA,
        tree = rules.refine.output.time_tree,
        extraction_script = "substitutions/extract_substitutions.ipynb"
    output:
        substitutions = "substitutions/substitutions.csv"
    shell:
        """
        jupyter execute {input.extraction_script}
        """


rule translate:
    input:
        time_tree = rules.refine.output.time_tree,
        ancestral_sequences = rules.ancestral.output.ancestral_sequences,
        reference = REFERENCE,
        genbank = GENBANK,
    output:
        amino_acid_substitutions = temp( "../nextstrain/aa-muts.json" )
    shell:
        """
        augur translate \
            --tree {input.time_tree} \
            --ancestral-sequences {input.ancestral_sequences} \
            --reference-sequence {input.genbank} \
            --vcf-reference {input.reference} \
            --output-node-data {output.amino_acid_substitutions}
        """

rule traits:
    input:
        time_tree = rules.refine.output.time_tree,
        metadata = rules.combine_all_metadata.output.md
    params:
        id_column = "taxa",
        columns = ["continent", "country", "IncAC"]
    output:
        traits_data = temp( "../nexstrain/traits.json" )
    shell:
        """
        augur traits \
            --tree {input.time_tree} \
            --metadata {input.metadata} \
            --metadata-id-columns {params.id_column} \
            --columns {params.columns} \
            --output-node-data {output.traits_data}
        """

rule export:
    input:
        time_tree = rules.refine.output.time_tree,
        branch_lengths = rules.refine.output.branch_lengths,
        nucleotide_subs = rules.ancestral.output.nucleotide_substitutions,
        amino_acid_subs = rules.translate.output.amino_acid_substitutions,
        traits = rules.traits.output.traits_data,
        metadata = rules.combine_all_metadata.output.md,
        auspice_config = AUSPICE_CONFIG,
        lat_longs = LAT_LONGS
    output:
        auspice_file = "../nextstrain/auspice.json"
    shell:
        """
        augur export v2 \
            --tree {input.time_tree} \
            --auspice-config {input.auspice_config} \
            --lat-longs {input.lat_longs} \
            --metadata {input.metadata} \
            --metadata-id-columns "taxa" \
            --node-data {input.branch_lengths} {input.nucleotide_subs} {input.amino_acid_subs} {input.traits} \
            --output {output.auspice_file}
        """
