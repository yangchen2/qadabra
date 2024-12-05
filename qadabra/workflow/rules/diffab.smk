import os

da_args = ["table", "metadata"]
da_params = ["factor_name", "target_level", "reference_level", "confounders"]

rule biom_to_qza:
    input:
        unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        qza="results/{dataset}/input_data/qza/{dataset}.qza"
    log:
        "log/{dataset}/biom_to_qza.log"
    conda:
        "../envs/q2-birdman-dev.yaml"
    shell:
        """
        # Log the current date and time
        echo "Start time: $(date)" >> {log} 2>&1

        # Log the QIIME 2 version
        echo "Running QIIME 2 in environment:" >> {log} 2>&1
        qiime --version >> {log} 2>&1
        
        # Log all details from `qiime info`
        echo "Complete QIIME 2 environment details:" >> {log} 2>&1
        qiime info >> {log} 2>&1

        # Check if Birdman is listed in the plugins
        echo "Running qiime tools import..." >> {log} 2>&1

        qiime tools import \
          --input-path {input.table} \
          --type 'FeatureTable[Frequency]' \
          --output-path {output.qza}

        # Log the end date and time
        echo "Complete! QZA file outputed." >> {log} 2>&1  
        """


rule birdman:
    input:
        table="results/{dataset}/input_data/qza/{dataset}.qza",
        metadata=lambda wc: get_dataset_cfg(wc, da_args)["metadata"]
    output:
        raw_results = "results/{dataset}/tools/birdman/raw_results.qza",
    log:
        "log/{dataset}/birdman.log"
    params:
        lambda wc: get_dataset_cfg(wc, da_params),
        formula=get_birdman_formula,
        outdir=lambda wc, output: os.path.dirname(output[0])
    shell: """
    echo 'Initializing Conda...' >> {log} 2>&1 &&
    source ~/miniforge3/etc/profile.d/conda.sh &&

    set +u
    conda activate q2-birdman-dev
    set -u

    echo 'Activated Conda Environment: q2-birdman-dev' >> {log} 2>&1 &&

    # Log the current date and time
    echo 'Start time: $(date)' >> {log} 2>&1 &&

    # Log the QIIME 2 version
    echo 'Running QIIME 2 in environment:' >> {log} 2>&1 &&
    qiime --version >> {log} 2>&1 &&

    # Log all details from `qiime info`
    echo 'Complete QIIME 2 environment details:' >> {log} 2>&1 &&
    qiime info >> {log} 2>&1 &&

    echo 'Removing existing directory if it exists...' >> {log} 2>&1 &&
    rm -rf {output} &&

    # Run the Birdman analysis
    echo 'Running Birdman analysis...' >> {log} 2>&1 &&
    qiime birdman run \
        --i-table {input[0]} \
        --m-metadata-file {input[1]} \
        --p-formula "{params.formula}" \
        --o-output-dir {output} \
        --verbose >> {log} 2>&1 &&

    # Confirm success
    echo 'Birdman analysis completed successfully!' >> {log} 2>&1
    """


rule export_birdman_output:
    input:
        output_qza = "results/{dataset}/tools/birdman/raw_results.qza"
    output:
        output_dir = directory("results/{dataset}/tools/birdman/raw_results"),
        output_file = "results/{dataset}/tools/birdman/raw_results/metadata.tsv"
    log:
        "log/{dataset}/export_birdman_output.log"
    conda:
        "../envs/q2-birdman-dev.yaml"
    shell:
        """
        # Log the current date and time
        echo "Start time: $(date)" >> {log} 2>&1

        # Log the QIIME 2 version
        echo "Running QIIME 2 in environment:" >> {log} 2>&1
        qiime --version >> {log} 2>&1
        
        # Log all details from `qiime info`
        echo "Complete QIIME 2 environment details:" >> {log} 2>&1
        qiime info >> {log} 2>&1

        qiime tools export \
          --input-path {input.output_qza} \
          --output-path {output.output_dir}

        # Log the end date and time
        echo "QZA file uzipped!" >> {log} 2>&1  
        """


rule extract_birdman_output:
    input:
        input_file="results/{dataset}/tools/birdman/raw_results/metadata.tsv"
    output:
        output_file="results/{dataset}/tools/birdman/raw_results/differentials_raw.tsv"
    log:
        "log/{dataset}/extract_birdman_output.log"
    shell:
        """
        echo "Start time: $(date)" >> {log} 2>&1
        mv {input.input_file} {output.output_file}
        echo "Extraction complete." >> {log} 2>&1
        """


rule remove_birdman_string_row:
    input:
        "results/{dataset}/tools/birdman/raw_results/differentials_raw.tsv"
    output:
        "results/{dataset}/tools/birdman/differentials.tsv"
    log:
        "log/{dataset}/remove_birdman_string_row.log"
    shell:
        """
        echo "Start time: $(date)" >> {log} 2>&1
        sed '2d' {input} > {output}
        echo "Row removal complete." >> {log} 2>&1
        """



rule deseq2:
    input:
        unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/deseq2/differentials.tsv",
        "results/{dataset}/tools/deseq2/results.rds",
    log:
        "log/{dataset}/deseq2.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/deseq2.R"


rule ancombc:
    input:
        unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
         "results/{dataset}/tools/ancombc/differentials.tsv",
         "results/{dataset}/tools/ancombc/results.rds",
    log:
        "log/{dataset}/ancombc.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/ancombc.R"


rule aldex2:
    input:
        unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/aldex2/differentials.tsv",
        "results/{dataset}/tools/aldex2/results.rds",
    log:
        "log/{dataset}/aldex2.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/aldex2.R"


rule edger:
    input:
        unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/edger/differentials.tsv",
        "results/{dataset}/tools/edger/results.rds",
    log:
        "log/{dataset}/edger.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/edger.R"


rule songbird:
    input:
       unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/songbird/differentials.tsv",
    log:
        "log/{dataset}/songbird.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params),
        epochs=config["songbird_params"]["epochs"],
        diff_prior=config["songbird_params"]["differential_prior"],
        formula=get_songbird_formula,
        outdir=lambda wc, output: os.path.dirname(output[0]),
    conda:
        "../envs/qadabra-songbird.yaml"
    shell:
        """
        songbird multinomial \
            --input-biom {input.table} \
            --metadata-file {input.metadata} \
            --formula "{params.formula}" \
            --epochs {params.epochs} \
            --differential-prior {params.diff_prior} \
            --summary-interval 1 \
            --min-feature-count 0 \
            --min-sample-count 0 \
            --random-seed 1 \
            --summary-dir {params.outdir} > {log} 2>&1
        """


rule maaslin2:
    input:
       unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        diff_file="results/{dataset}/tools/maaslin2/differentials.tsv",
        out_dir=directory("results/{dataset}/tools/maaslin2/output"),
    log:
        "log/{dataset}/maaslin2.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/maaslin2.R"


rule metagenomeseq:
    input:
       unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/metagenomeseq/differentials.tsv",
        "results/{dataset}/tools/metagenomeseq/results.rds",
    log:
        "log/{dataset}/metagenomeseq.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/metagenomeseq.R"


rule corncob:
    input:
       unpack(lambda wc: get_dataset_cfg(wc, da_args))
    output:
        "results/{dataset}/tools/corncob/differentials.tsv",
        "results/{dataset}/tools/corncob/results.rds",
    log:
        "log/{dataset}/corncob.log",
    params:
        lambda wc: get_dataset_cfg(wc, da_params)
    conda:
        "../envs/qadabra-da-R.yaml"
    script:
        "../scripts/R/corncob.R"


rule process_differentials:
    input:
        "results/{dataset}/tools/{tool}/differentials.tsv",
    output:
        "results/{dataset}/tools/{tool}/differentials.processed.tsv",
    log:
        "log/{dataset}/process_differentials.{tool}.log",
    params:
        col=get_diffab_tool_columns
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_differentials.py"


rule concatenate_differentials:
    input:
        expand(
            "results/{{dataset}}/tools/{tool}/differentials.processed.tsv",
            tool=config["tools"]
        ),
    output:
        "results/{dataset}/concatenated_differentials.tsv",
    log:
        "log/{dataset}/concatenate_differentials.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/concatenate_differentials.py"


rule process_pvalues:
    input:
        "results/{dataset}/tools/{tool}/differentials.tsv",
    output:
        "results/{dataset}/tools/{tool}/pvalues.processed.tsv",
    log:
        "log/{dataset}/process_pvalues.{tool}.log",
    params:
        col=get_pvalue_tool_columns
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/process_pvalues.py"


rule concatenate_pvalues:
    input:
        expand(
            "results/{{dataset}}/tools/{tool}/pvalues.processed.tsv",
            tool=config["ptools"]
        ),
    output:
        "results/{dataset}/concatenated_pvalues.tsv",
    log:
        "log/{dataset}/concatenate_pvalues.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/concatenate_pvalues.py"


rule concatenate_all_results:
    input:
        concatenated_pvalues="results/{dataset}/concatenated_pvalues.tsv",
        concatenated_differentials="results/{dataset}/concatenated_differentials.tsv",
    output:
        "results/{dataset}/qadabra_all_result.tsv",
    log:
        "log/{dataset}/qadabra_all_result.log",
    conda:
        "../envs/qadabra-default.yaml"
    script:
        "../scripts/concatenate_all_results.py"