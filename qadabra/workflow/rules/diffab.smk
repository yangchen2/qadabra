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
        source activate /Users/yangchen/miniforge3/envs/q2-birdman-dev
        qiime tools import \
          --input-path {input.table} \
          --type 'FeatureTable[Frequency]' \
          --output-path {output.qza}
        """


rule birdman:
    input:
        table="results/{dataset}/input_data/qza/{dataset}.qza",
        metadata=lambda wc: get_dataset_cfg(wc, da_args)["metadata"]
    output:
        directory("results/{dataset}/tools/birdman/raw_results/"),
        beta_var="results/{dataset}/tools/birdman/raw_results/results/beta_var.tsv"
    log:
        "log/{dataset}/birdman.log"
    params:
        lambda wc: get_dataset_cfg(wc, da_params),
        formula=get_birdman_formula,
    conda:
        "../envs/q2-birdman-dev.yaml"
    shell:
        """
        source activate /Users/yangchen/miniforge3/envs/q2-birdman-dev

        # Log the current date and time
        echo "Start time: $(date)" >> {log} 2>&1

        # Log the QIIME 2 version
        echo "Running QIIME 2 in environment:" >> {log} 2>&1
        qiime --version >> {log} 2>&1
        
        # Log all details from `qiime info`
        echo "Complete QIIME 2 environment details:" >> {log} 2>&1
        qiime info >> {log} 2>&1

        # Check if Birdman is listed in the plugins
        echo "Checking for Birdman plugin:" >> {log} 2>&1
        if ! qiime info | grep birdman; then
            echo "Error: Birdman plugin is not installed or unavailable in this environment!" >> {log} 2>&1
            echo "End time: $(date)" >> {log} 2>&1
            exit 1
        fi

        # Run the Birdman analysis
        echo "Running Birdman analysis:" >> {log} 2>&1
        qiime birdman run \
            --i-table {input.table} \
            --m-metadata-file {input.metadata} \
            --p-formula "{params.formula}" \
            --o-output-dir {output[0]} >> {log} 2>&1


        # Confirm success
        if [ $? -eq 0 ]; then
            echo "Birdman analysis completed successfully!" >> {log} 2>&1
        else:
            echo "Error: Birdman analysis failed!" >> {log} 2>&1
            echo "End time: $(date)" >> {log} 2>&1
            exit 1
        fi

        # Log the end date and time
        echo "End time: $(date)" >> {log} 2>&1
        """


rule extract_birdman_results:
    input:
        directory("results/{dataset}/tools/birdman/raw_results/results/beta_var.tsv")
    output:
        beta_var="results/{dataset}/tools/birdman/differentials.tsv",
    log:
        "log/{dataset}/extract_birdman_results.log"
    shell:
        """
        cp {input}/results/beta_var.tsv {output.beta_var}
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