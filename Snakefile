workdir: "analyses/"

# Inputs.
METADATA = ["../data/supplemental_data1.csv", "../data/supplemental_data2.csv"]
ML_TREE = "../data/2024.08.06_ML_GTR.tree"
BEAST_LOG = "../beast-analyses/2024-08-06_constant_relaxed.combined.log"
BEAST_MCC = "../beast-analyses/2024-08-06_constant_relaxed.mcc.tree"
BEAST_TREES = "../beast-analyses/2024-08-06_constant_relaxed.combined.down.trees"
BEAST_DISCRETE_LOG = "../beast-analyses/2024-08-06_constant_relaxed_discrete.log"
BEAST_DISCRETE_MCC = "../beast-analyses/2024-08-06_constant_relaxed_discrete.mcc.tree"
BEAST_DISCRETE_TREES = "../beast-analyses/2024-08-06_constant_relaxed_discrete.trees"
BEAST_MARKOV_JUMPS = "../beast-analyses/2024-08-06_constant_relaxed_discrete.Location.csv"
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
