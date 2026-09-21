# Serotypes
serotype = ["Dengue_1", "Dengue_2"]

# "even" or "proportional" for "background_method"
subsampling_counts = {
    "Dengue_1": {
        "local": 300,
        "background": 200,
        "background_method": "even"
    },
    "Dengue_2": {
        "local": "all",
        "background": 200,
        "background_method": "even"
    }
}

# extra metadata and fasta
extra_metadata = "_data/extra/test_dengue.tsv"
extra_fasta = "_data/extra/genomic_extra_dengue.fna"


rule all:
    input:
        # expand("_results/s06_add_genotype_information/Unaligned_{serotype}_infoTbl_with_genotype.csv", serotype = serotype),
        # expand("_results/s07_sequence_alignment/Aligned_{serotype}_nextalign.aligned.fasta", serotype = serotype)
        expand("_results/s11_treebuilding/subsampled_{serotype}_cleaned.treefile", serotype = serotype)
    output:
        "dag.png"
    log:
        "_logs/step00/all.log"
    shell:
        """
        snakemake --dag | dot -Tpng > dag.png 
        """

# Rule for creating necessary directories
rule create_directories:
    output:
        directory("_logs"),
        directory("_data"),
        directory("_results")
    shell:
        "mkdir -p _logs _data _results"

# Step 1: Acquisition of Genomic Data and Metadata from GenBank
rule s01_acquire_data:
    output:
        zip="_data/genbank_data.zip",
        metadata="_data/metadata.tsv",
        fasta="_data/ncbi_dataset/data/genomic.fna"
    params:
        date = "2005-01-01",
    log:
        "_logs/s01_acquire_data/s01_acquire_data.log"
    message:
        "Acquiring Genomic Data and Metadata from NCBI database"
    conda:
        "config/ncbi_datasets_env.yaml"
    shell:
        """
        datasets download virus genome taxon "Dengue Virus" --filename {output.zip} --released-after {params.date} &> {log}
        datasets summary virus genome taxon "Dengue Virus" --released-after {params.date} --as-json-lines | dataformat tsv virus-genome > {output.metadata} 2>> {log}
        
        # Unzip the downloaded data
        unzip -o -d _data {output.zip}
        """

# Step 2: Clean metadata and FASTA
rule s02_process_genbank_data:
    input:
        script="code/02_clean_metadata_and_fasta_general.R",
        metadata="_data/metadata.tsv",
        fasta="_data/ncbi_dataset/data/genomic.fna",
        extra_metadata=extra_metadata,
        extra_fasta=extra_fasta
    output:
        fasta_files ="_results/s02_process_genbank_data/Unaligned.fasta",
        info_tables_txt ="_results/s02_process_genbank_data/infoTbl.txt",
        info_tables_csv ="_results/s02_process_genbank_data/infoTbl.csv"
    params:
        start_date = "2000-01-02",
        end_date = "2025-01-01",
        host = "Homo sapiens"
    log:
        "_logs/s02_process_genbank_data/s02_process_genbank_data.log"
    message:
        "Processing and cleaning data downloaded from NCBI"
    shell:
        """
        Rscript {input.script} \
            --metadata {input.metadata} \
            --fasta {input.fasta} \
            --extra_metadata {input.extra_metadata} \
            --extra_fasta {input.extra_fasta} \
            --start_date {params.start_date} \
            --end_date {params.end_date} \
            --host "{params.host}" \
            --outfile_fasta {output.fasta_files} \
            --outfile_csv {output.info_tables_csv} \
            --outfile_tsv {output.info_tables_txt} \
            > {log} 2>&1
        """

# Step 3: Filter for sequences from SEA
rule s03_filter_dengue_data:
    input:
        script="code/03_filter_SEA.R",
        metadata="_results/s02_process_genbank_data/infoTbl.csv",
        metadata_vietnam="_data/sequences_vietnam.csv",
        metadata_china="_data/sequences_china.csv",
        metadata_extra="_data/extra/metadata_extra.tsv",
        fasta="_results/s02_process_genbank_data/Unaligned.fasta",
    output:
        fasta = "_results/s03_filter_dengue_data/Unaligned_SEA.fasta",
        csv = "_results/s03_filter_dengue_data/Unaligned_SEA_infoTbl.csv",
        tsv = "_results/s03_filter_dengue_data/Unaligned_SEA_infoTbl.tsv"
    params:
        out_prefix="_results/s03_filter_dengue_data/Unaligned_SEA"
    log:
        "_logs/s03_filter_dengue_data/filter_data.log"
    message:
        "Filter data for countries in SEA and select only sequences from china and Vietnam with known geo-coded sequences"
    shell:
        """
        Rscript {input.script} \
            --metadata {input.metadata} \
            --metadata_vietnam {input.metadata_vietnam} \
            --metadata_china {input.metadata_china} \
            --metadata_extra {input.metadata_extra} \
            --fasta {input.fasta} \
            --outfile {params.out_prefix} \
            > {log} 2>&1
        """

# Step 4: Split into serotype, add serotypes to sequence name and generate sequence specific metadata
rule s04_split_dengue_data:
    input:
        script="code/04_split_dengue.R",
        metadata="_results/s03_filter_dengue_data/Unaligned_SEA_infoTbl.csv",
        fasta="_results/s03_filter_dengue_data/Unaligned_SEA.fasta"
    output:
        fasta = "_results/s04_split_dengue_data/Unaligned_{serotype}.fasta",
        csv = "_results/s04_split_dengue_data/Unaligned_{serotype}_infoTbl.csv",
        tsv = "_results/s04_split_dengue_data/Unaligned_{serotype}_infoTbl.txt"
    params:
        out_prefix="_results/s04_split_dengue_data/Unaligned_",
        serotype = "{serotype}"
    log:
        "_logs/s04_split_dengue_data/process_data_{serotype}.log"
    message:
        "Split into {wildcards.serotype}, add serotypes to sequence name and generate sequence specific metadata"
    shell:
        """
        Rscript {input.script} \
            --metadata {input.metadata} \
            --fasta {input.fasta} \
            --serotype {params.serotype} \
            --outfile {params.out_prefix} \
            > {log} 2>&1
        """


# Step 5: Assign Genotype Using Nextclade
rule s05_assign_serotype_and_genotype_nextclade:
    input:
        dataset="_data/nextclade/denvLineages/{serotype}",
        fasta = "_results/s04_split_dengue_data/Unaligned_{serotype}.fasta"
    output:
        csv="_results/s05_assign_serotype_and_genotype_nextclade/nextclade_{serotype}.csv"
    log:
        "_logs/s05_assign_serotype_and_genotype_nextclade/assign_serotype_and_genotype_{serotype}.log"
    message:
        "Assign genotype for {wildcards.serotype} with Nextclade"
    shell:
        """
        nextclade run \
            --input-dataset {input.dataset} \
            --output-csv {output.csv} \
            {input.fasta} \
            > {log} 2>&1
        """


# Step 6: Add genotype information to metadata
rule s06_add_genotype_information:
    input:
        script="code/06_add_genotype_information_to_metadata.R",
        metadata="_results/s04_split_dengue_data/Unaligned_{serotype}_infoTbl.csv",
        genotype="_results/s05_assign_serotype_and_genotype_nextclade/nextclade_{serotype}.csv"
    output:
        csv="_results/s06_add_genotype_information/Unaligned_{serotype}_infoTbl_with_genotype.csv"
    log:
        "_logs/s06_add_genotype_information/Unaligned_{serotype}_add_genotype.log"
    message:
        "Add genotype information to metadata for {wildcards.serotype}"
    shell:
        """
        Rscript {input.script} \
            --metadata {input.metadata} \
            --genotype {input.genotype} \
            --outfile_csv {output.csv} \
            > {log} 2>&1
        """

        
# Step 7: Sequence alignment
rule s07_sequence_alignment:
    input:
        sequences="_results/s04_split_dengue_data/Unaligned_{serotype}.fasta",
        reference="_data/nextclade/denvLineages/{serotype}/reference.fasta",
        genemap="_data/nextclade/denvLineages/{serotype}/genome_annotation.gff3"
    output:
        fasta="_results/s07_sequence_alignment/Aligned_{serotype}_nextalign.aligned.fasta"
    log:
        "_logs/s07_sequence_alignment/sequence_alignment_{serotype}.log"
    conda:
        "config/alignment.yaml"
    message:
        "Running sequence alignment for {wildcards.serotype}"
    shell:
        """
        nextalign run \
            {input.sequences} \
            --reference={input.reference} \
            --genemap={input.genemap} \
            --output-fasta={output.fasta} \
            > {log} 2>&1
        """
        
# nextclade run \
#             {input.sequences} \
#             --input-ref={input.reference} \
#             --input-annotation={input.genemap}\
#             --output-fasta={output.fasta} \
#             > {log} 2>&1

# Step 7: Sequence alignment
# rule old_s07_sequence_alignment:
#     input:
#         sequences="_results/s04_split_dengue_data/Unaligned_{serotype}.fasta",
#         reference="_data/reference_genomes/reference_{serotype}.fasta",
#         genemap="_data/genemap/genemap_{serotype}.gff"
#     output:
#         fasta="_results/s07_sequence_alignment/Aligned_{serotype}_nextalign.aligned.fasta"
#     log:
#         "_logs/s07_sequence_alignment/sequence_alignment_{serotype}.log"
#     conda:
#         "config/alignment.yaml"
#     message:
#         "Running sequence alignment for {wildcards.serotype}"
#     shell:
#         """
#         nextalign run \
#             {input.sequences} \
#             --reference={input.reference} \
#             --genemap={input.genemap} \
#             --output-fasta={output.fasta} \
#             > {log} 2>&1
#         """


# Step 8: Segregating E gene and Whole Genomes and performing quaility control
rule s08_split_genome_and_QC:
    input:
        script = "code/08_seperate_EG_and_WG.R",
        fasta = "_results/s07_sequence_alignment/Aligned_{serotype}_nextalign.aligned.fasta"
    output:
        E_gene_dir = "_results/s08_split_genome_and_QC/{serotype}_EG.fasta",
        WG_gene_dir = "_results/s08_split_genome_and_QC/{serotype}_WG.fasta"
    params:
        wg_threshold = 0.31,
        eg_threshold = 0.31
    log:
        "_logs/s08_split_genome_and_QC/Segregating_{serotype}.log"
    message:
        "Segregating E gene and whole genomes from aligned {wildcards.serotype} sequences and performing quality control."
    shell:
        """
        Rscript {input.script} \
            --fasta {input.fasta} \
            --WG_threshold {params.wg_threshold} \
            --EG_threshold {params.eg_threshold} \
            --outfile_fasta_eg {output.E_gene_dir} \
            --outfile_fasta_wg {output.WG_gene_dir} \
            > {log} 2>&1
        """

# Step 9: Subsampling DENVs
rule s09_subsample_denv:
    input:
        script = "code/09_subsampler.R",
        metadata_file = "_results/s06_add_genotype_information/Unaligned_{serotype}_infoTbl_with_genotype.csv",
        fasta_file = "_results/s08_split_genome_and_QC/{serotype}_EG.fasta",
        location_local = "_data/locations_of_interest.csv"
    output:
        subsample_fasta = "_results/s09_subsample_denv/subsampled_{serotype}.fasta",
        subsample_csv = "_results/s09_subsample_denv/subsampled_{serotype}_infoTbl.csv",
        subsample_txt = "_results/s09_subsample_denv/subsampled_{serotype}_infoTbl.tsv"
    params:
        number_sequences_local = lambda wc: subsampling_counts[wc.serotype]["local"],
        number_sequences_background = lambda wc: subsampling_counts[wc.serotype]["background"],
        sampling_method = lambda wc: subsampling_counts[wc.serotype]["background_method"],
        time_interval = "Year",
        serotype = "{serotype}",
        out_prefix = "_results/s09_subsample_denv/subsampled_{serotype}",
        out_dir = "_results/s09_subsample_denv/plots/"
    log:
        "_logs/s09_subsample_denv/subsample_{serotype}.log"
    message:
        "Subsampling {wildcards.serotype} virus E gene sequences based on specified criteria."
    shell:
        """
        mkdir -p {params.out_dir}
        Rscript {input.script} \
            --metadata {input.metadata_file} \
            --fasta {input.fasta_file} \
            --location_local {input.location_local} \
            --time_interval {params.time_interval} \
            --number_sequences_local {params.number_sequences_local} \
            --number_sequences_background {params.number_sequences_background} \
            --sampling_method {params.sampling_method} \
            --serotype {params.serotype} \
            --out_prefix {params.out_prefix} \
            --output_dir {params.out_dir} \
            > {log} 2>&1
        """

# Step 9b: remove sequences with stop codons
rule s09b_hyphy_cln:
    input:
        fasta_file = "_results/s09_subsample_denv/subsampled_{serotype}.fasta"
    output:
        cleaned_fasta = "_results/s09b_hyphy_cln/subsampled_{serotype}_cleaned.fasta",
    log:
        "_logs/s09b_hyphy_cln/hyphy_cln_{serotype}.log"    
    shell:
        """
        hyphy cln \
            Universal \
            {input.fasta_file} \
            "Yes/Yes" \
            {output.cleaned_fasta} \
            > {log} 2>&1
        """

# Step 10: Correct metadata and fasta files into the correct format for iqtree and treetime  
rule s10_reformatting:
    input:
        script = "code/10_reformatting_iqtree_treetime.R",
        # fasta_file = "_results/s09_subsample_denv/subsampled_{serotype}.fasta",
        fasta_file = "_results/s09b_hyphy_cln/subsampled_{serotype}_cleaned.fasta",
        metadata_file = "_results/s09_subsample_denv/subsampled_{serotype}_infoTbl.csv",
    output:
        cleaned_fasta = "_results/s10_reformatting/subsampled_{serotype}_cleaned.fasta",
        cleaned_metadata = "_results/s10_reformatting/subsampled_{serotype}_cleaned_infoTbl.csv"    
    log:
        "_logs/s10_reformatting/reformatting_{serotype}.log"    
    message:
        "Correct {wildcards.serotype} metadata and fasta files into the correct format for iqtree and treetime"
    shell:
        """
        Rscript {input.script} \
            --metadata {input.metadata_file} \
            --fasta {input.fasta_file} \
            --output_dir_fasta {output.cleaned_fasta} \
            --output_dir_csv {output.cleaned_metadata} \
            > {log} 2>&1
        """

# Step 11: Treebuilding
rule s11_treebuilding:
    input:
        aln = "_results/s10_reformatting/subsampled_{serotype}_cleaned.fasta"
    output:
        tree = "_results/s11_treebuilding/subsampled_{serotype}_cleaned.treefile"
    params:
        prefix = "_results/s11_treebuilding/subsampled_{serotype}_cleaned",
        model = "GTR+I+G4"
    conda:
        "config/iqtree.yaml"
    log:
        "_logs/s11_treebuilding/iqtree_{serotype}.log"
    message:
        "Inferring maximum likelihood phylogenetic trees for {wildcards.serotype} using IQ-TREE."
    shell:
        """
        iqtree2 \
            -nt AUTO \
            -s {input.aln} \
            -m {params.model} \
            -pre {params.prefix} \
            -redo \
            > {log} 2>&1
        """

# Step 12: Remove duplicates with hyphy analyses
rule s12_hyphy_remove_dupes:
    input:
        fasta_file = "_results/s10_reformatting/subsampled_{serotype}_cleaned.fasta",
        tree_file = "_results/s11_treebuilding/subsampled_{serotype}_cleaned.treefile"
    output:
        outfile = "_results/s12_hyphy_remove_dupes/subsampled_{serotype}_deduped_tree.fasta"
    params:
        workdir = "_results/s12_hyphy_remove_dupes/temp/"
    conda:
        "config/hyphy.yaml"
    log:
        "_logs/s12_hyphy_remove_dupes/hyphy_remove_dupes_{serotype}.log"
    shell:
        """
        mkdir -p {params.workdir}
        cp {input.fasta_file} {params.workdir}
        cp {input.tree_file} {params.workdir}
        cp hyphy-analyses/remove-duplicates/remove-duplicates.bf {params.workdir}

        (cd {params.workdir};
        hyphy remove-duplicates.bf \
            --msa $(basename {input.fasta_file}) \
            --tree $(basename {input.tree_file}) \
            --output ../$(basename {output.outfile}) \
            ENV="DATA_FILE_PRINT_FORMAT=9" \
        ) > {log} 2>&1
        """

# Step 13: Filter strains from VN only
rule s13_filter_vn_strains:
    input:
        script = "code/13_filter_vn_strains.R",
        csv_file = "_results/s10_reformatting/subsampled_{serotype}_cleaned_infoTbl.csv"
    output:
        outfile = "_results/s13_filter_vn_strains/subsampled_{serotype}_filtered_vn_strains.txt"
    log:
        "_logs/s13_filter_vn_strains/filter_vn_strains_{serotype}.log"
    shell:
        """
        Rscript {input.script} \
            --csv {input.csv_file} \
            --outfile {output.outfile} \
            > {log} 2>&1
        """

# Step 14: Label tree with hyphy analyses
rule s14_hyphy_label_trees:
    input:
        treefile = "_results/s12_hyphy_remove_dupes/subsampled_{serotype}_deduped_tree.fasta",
        filtered_strains = "_results/s13_filter_vn_strains/subsampled_{serotype}_filtered_vn_strains.txt"
    output:
        labelled_tree = "_results/s14_hyphy_label_trees/subsampled_{serotype}_deduped_labelled_tree.fasta"
    params:
        fasta_dat = "_results/s14_hyphy_label_trees/temp/fasta_seq_{serotype}.temp",
        temp_file = "_results/s14_hyphy_label_trees/temp/tree_{serotype}.temp",
        workdir = "_results/s14_hyphy_label_trees/temp/"
    conda:
        "config/hyphy.yaml"
    log:
        "_logs/s14_hyphy_label_trees/hyphy_label_trees_{serotype}.log"
    shell:
        """
        mkdir -p {params.workdir}
        cp {input.filtered_strains} {params.workdir}
        cp hyphy-analyses/LabelTrees/label-tree.bf {params.workdir}
        head -n -1 {input.treefile} > {params.workdir}$(basename {params.fasta_dat})
        tail -1 {input.treefile} > {params.workdir}$(basename {input.treefile})
        
        (cd {params.workdir};
        hyphy label-tree.bf \
            --tree $(basename {input.treefile}) \
            --list $(basename {input.filtered_strains}) \
            --output $(basename {params.temp_file}) \
            --internal-nodes None;
        echo "Finished labelling. Adding labelled nwk to fasta.";
        cat $(basename {params.fasta_dat}) > ../$(basename {output.labelled_tree});
        cat $(basename {params.temp_file}) >> ../$(basename {output.labelled_tree});
        ) > {log} 2>&1
        """


# Step 15a: HyPhy selection pressure w/ MEME
rule s15a_hyphy_meme:
    input:
        alignment_file = "_results/s14_hyphy_label_trees/subsampled_{serotype}_deduped_labelled_tree.fasta"
    output:
        outfile = "_results/s15a_hyphy_meme/meme_{serotype}.json"
    params:
        branches = "Foreground"
    threads:
        os.cpu_count() * (1/len(serotype))
    conda:
        "config/hyphy.yaml"
    log:
        "_logs/s15a_hyphy_meme/hyphy_meme_{serotype}.log"
    shell:
        """
        mpirun -np {threads} \
            HYPHYMPI meme \
            --alignment {input.alignment_file} \
            --branches {params.branches} \
            --output {output.outfile} \
            > {log} 2>&1
        """

# Step 15b: HyPhy selection pressure w/ FUBAR
rule s15b_hyphy_fubar:
    input:
        alignment_file = "_results/s14_hyphy_label_trees/subsampled_{serotype}_deduped_labelled_tree.fasta"
    output:
        outfile = "_results/s15b_hyphy_fubar/fubar_{serotype}.json",
        cache = "_results/s15b_hyphy_fubar/fubar_{serotype}.cache"
    threads:
        os.cpu_count() * (1/len(serotype))
    conda:
        "config/hyphy.yaml"
    log:
        "_logs/s15b_hyphy_fubar/hyphy_fubar_{serotype}.log"
    shell:
        """
        mkdir -p _results/s15b_hyphy_fubar

        mpirun -np {threads} \
            HYPHYMPI fubar \
            --alignment {input.alignment_file} \
            > {log} 2>&1
        
        mv {input.alignment_file}.FUBAR.json {output.outfile}
        mv {input.alignment_file}.FUBAR.cache {output.cache}
        """

# Step 15c: HyPhy selection pressure w/ SLAC
rule s15c_hyphy_slac:
    input:
        alignment_file = "_results/s14_hyphy_label_trees/subsampled_{serotype}_deduped_labelled_tree.fasta"
    output:
        outfile = "_results/s15c_hyphy_slac/slac_{serotype}.json",
    threads:
        os.cpu_count() * (1/len(serotype))
    conda:
        "config/hyphy.yaml"
    log:
        "_logs/s15c_hyphy_slac/hyphy_slac_{serotype}.log"
    shell:
        """
        mpirun -np {threads} \
            HYPHYMPI slac \
            --alignment {input.alignment_file} \
            --output {output.outfile} \
            > {log} 2>&1
        """

# Step 15d: HyPhy selection pressure w/ FEL
rule s15d_hyphy_fel:
    input:
        alignment_file = "_results/s14_hyphy_label_trees/subsampled_{serotype}_deduped_labelled_tree.fasta"
    output:
        outfile = "_results/s15d_hyphy_fel/fel_{serotype}.json",
    threads:
        os.cpu_count() * (1/len(serotype))
    conda:
        "config/hyphy.yaml"
    log:
        "_logs/s15d_hyphy_fel/hyphy_fel_{serotype}.log"
    shell:
        """
        mpirun -np {threads} \
            HYPHYMPI fel \
            --alignment {input.alignment_file} \
            --output {output.outfile} \
            > {log} 2>&1
        """

# Step 15e: HyPhy selection pressure w/ Constrast-FEL
rule s15e_hyphy_contrast_fel:
    input:
        alignment_file = "_results/s14_hyphy_label_trees/subsampled_{serotype}_deduped_labelled_tree.fasta"
    output:
        outfile = "_results/s15e_hyphy_contrast_fel/contrast_fel_{serotype}.json"
    params:
        branch_set = "Foreground"
    threads:
        os.cpu_count() * (1/len(serotype))
    conda:
        "config/hyphy.yaml"
    log:
        "_logs/s15e_hyphy_contrast_fel/hyphy_contrast_fel_{serotype}.log"
    shell:
        """
        mpirun -np {threads} \
            HYPHYMPI contrast-fel \
            --alignment {input.alignment_file} \
            --branch-set {params.branch_set} \
            --output {output.outfile} \
            > {log} 2>&1
        """
