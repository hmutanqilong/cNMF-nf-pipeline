#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/****************************
 * paremeter definitions
 ***************************/

// input parameters
params.analysisDir = params.analysisDir ?: "./analysis"
params.count_fp = params.count_fp ?: error("The parameter 'count_fp' must be provided.")
params.name = params.name ?: "cnmf_run"

// parallelization parameters
params.total_workers = 20    // 总共20个worker

// k values for CNMF
params.k_values = params.k_values ?: "5,10,15"

// iterations for CNMF
params.n_iter = params.n_iter ?: 100

// number of high variance genes to use
params.num_genes = params.num_genes ?: 2000

// gene selection
params.cnmf_gene_selection = "top${params.num_genes}VariableGenes"

// output directory
params.outdir = params.outdir ?: "${params.analysisDir}/${params.cnmf_gene_selection}"

// computation time
params.prepare_time = params.prepare_time ?: "2h"
params.prepare_generalized_time = params.perpare_generalized_time ?: "2h"
params.factorize_worker_time = params.factorize_worker_time ?: "35h"
params.combine_time = params.combine_time ?: "2h"
params.kplot_time = params.kplot_time ?: "2h"

// computation memory
params.prepare_mem = params.prepare_mem ?: "10.G"
params.prepare_generalized_mem = params.perpare_generalized_mem ?: "2h"
params.factorize_worker_mem = params.factorize_worker_mem ?: "35h"
params.combine_mem = params.combine_mem ?: "2h"
params.kplot_mem = params.kplot_mem ?: "2h"


/****************************
 *        FUNCTIONS
 ***************************/

// function to parse k values
def parseKValues(k_values) {
    return (k_values instanceof java.util.List) ? k_values : k_values.split(',').collect { it as Integer }
}

def kList = parseKValues(params.k_values)

println ">>> Starting CNMF pipeline with the following parameters:"
println ">>> Analysis directory: ${params.analysisDir}"
println ">>> Count file path: ${params.count_fp}"
println ">>> Output directory: ${params.outdir}"
println ">>> K values: ${kList}"

/****************************
 *       PROCESSES
 ***************************/
// PROCESS A: Prepare initial data for CNMF
process prepare_initial {
    tag "prepare_initial"

    // Resource directives
    time params.prepare_time
    memory params.prepare_mem

    publishDir params.outdir, mode: 'copy', overwrite: true

    input:
    path h5ad_file

    output:
    path "prepare.done", emit: prepare_initial_done
    path "${params.name}/**", emit: prepared_files  // 输出整个目录结构

    script:
    """
    echo ">>> Running prepare_initial, output to: ${params.outdir}"

    cnmf prepare \
        --output-dir . \
        --name ${params.name} \
        -c ${params.count_fp} \
        -k 1 \
        --n-iter ${params.n_iter} \
        --numgenes ${params.num_genes}
    
    touch prepare.done
    """    
}

// PROCESS B: Run CNMF for generated k values parameter
process prepare_generalized {
    tag "prepare_generalized"
    
    time params.prepare_generalized_time
    memory params.prepare_generalized_mem

    publishDir params.outdir, mode: 'copy', overwrite: true // Ensure outputs are copied to the same output directory

    input:
    path h5ad_file
    path prepare_initial_done
    path prepared_files  // 接收prepare_initial的文件

    output:
    path "prepare_generalized.done", emit: prepare_generalized_done
    path "${params.name}/cnmf_tmp/**", emit: updated_files  // 输出更新后的目录结构


    script:
    """
    sleep 1m
    
    echo ">>> Running prepare_generalized"
    echo ">>> K values: ${kList.join(',')}"
    
        
    python ${projectDir}/scripts/cnmf.modified.py prepare_generalize \
        --output-dir . \
        --name ${params.name} \
        -c ${params.count_fp} \
        -k ${kList.join(' ')} \
        --n-iter ${params.n_iter} \
        --numgenes ${params.num_genes}
    

    echo ">>> prepare_generalized completed successfully"
    touch prepare_generalized.done

    """    
}

// PROCESS C: Run CNMF factorization for generated k values parameter
process factorize_worker {
    tag "factorize_k${k}_worker${worker_idx}"
    
    time params.factorize_worker_time
    memory params.factorize_worker_mem

    // 所有worker共享同一个输出目录
    stageInMode 'rellink'  // 使用相对链接减少复制
    scratch false  // 不使用scratch目录

    input:
    tuple val(k), val(worker_idx)  // k值和worker索引的组合
    each path(h5ad_file)           // 确保每个并行任务都能获得h5ad文件
    each path(prepare_generalized_done)  // 确保每个并行任务都能获得完成标记

    output:
    tuple val(k), val(worker_idx), path("factorize_k${k}_worker${worker_idx}.done"), emit: factorize_worker_done

    script:
    """

    echo ">>> Running factorization for k=${k}, worker=${worker_idx}"
    echo ">>> Working from: \$(pwd)"
    echo ">>> params.outdir: ${params.outdir}"
    echo ">>> Resolved output directory: ${file(params.outdir).toAbsolutePath()}"
    
    # 显示任务编号（帮助追踪并行度）
    total_combinations=\$((${params.k_values.split(',').size()} * ${params.total_workers}))
    current_task=\$((\$((${k} - ${kList.min()})) * ${params.total_workers} + ${worker_idx} + 1))
    echo ">>> Task \$current_task of \$total_combinations"

    # 验证必要文件存在
    param_file="${file(params.outdir).toAbsolutePath()}/${params.name}/cnmf_tmp/${params.name}.nmf_params.df.npz"
    if [ ! -f "\$param_file" ]; then
        echo ">>> ERROR: Parameter file not found: \$param_file"
        echo ">>> Listing contents of expected directory:"
        ls -la ${file(params.outdir).toAbsolutePath()}/${params.name}/cnmf_tmp/ 2>/dev/null || echo "Directory does not exist"
        exit 1
    fi
    
    echo ">>> Parameter file verified: \$param_file"
    echo ">>> Starting factorization..."


    # 直接使用publishDir目录，所有worker共享同一个目录结构
    cnmf factorize \
        --output-dir ${file(params.outdir).toAbsolutePath()} \
        --name ${params.name} \
        --worker-index ${worker_idx} \
        --total-workers ${params.total_workers} \
        --skip-completed-runs

    # 检查生成的文件
    echo ">>> Checking generated iteration files:"
    iter_files=\$(ls ${file(params.outdir).toAbsolutePath()}/${params.name}/cnmf_tmp/*spectra.k_${k}.iter_*.df.npz 2>/dev/null | wc -l || echo 0)
    echo ">>> Found \$iter_files iteration files for k=${k}"
    
    echo ">>> Worker ${worker_idx} completed for k=${k}"
    touch factorize_k${k}_worker${worker_idx}.done

    """
}

// PROCESS D: Collect results from all workers for a specific k value
process factorize_complete {
    tag "factorize_complete_k${k}"
    
    input:
    tuple val(k), path(worker_done_files)  // 收集所有worker的完成文件

    output:
    tuple val(k), path("factorize_complete_k${k}.done"), emit: factorize_complete

    script:
    """
    echo ">>> All workers completed for k=${k}"
    echo ">>> Worker files: ${worker_done_files}"
    
    # 检查所有worker是否都完成了
    worker_count=\$(ls factorize_k${k}_worker*.done | wc -l)
    echo ">>> Completed workers for k=${k}: \$worker_count"
    
    touch factorize_complete_k${k}.done
    """
}

// PROCESS E: cNMF combine results
process combine_results {
    tag "combine_results"
    
    time params.combine_time
    memory params.combine_mem

    input:
    path(factorize_complete_files)

    output:
    path("cnmf_combine.done"), emit: combine_done
    
    script:
    """
    echo ">>> Running cNMF combine step"
    
    cnmf combine \
        --output-dir ${file(params.outdir).toAbsolutePath()} \
        --name ${params.name}
    
    touch cnmf_combine.done
    """
}

// PROCESS F: k_selection plot
process kplot {
    tag "k_selection_plot"
    
    time params.kplot_time
    memory params.kplot_mem

    input:
    path(combine_done)  // 依赖combine步骤完成
    
    output:
    path("k_selection_plot.done"), emit: k_selection_plot_done
    
    script:
    """
    echo ">>> Starting CNMF k selection plot generation"
    
    echo ">>> Running CNMF k selection plot..."
    
    # 运行k selection plot命令
    cnmf k_selection_plot \
        --output-dir ${file(params.outdir).toAbsolutePath()} \
        --name ${params.name}
    
    touch cnmf_k_selection_plot.done
    """
}


/****************************
 *       WORKFLOWS
 ***************************/
workflow {
    h5ad_ch = Channel.fromPath(params.count_fp)
    // step 1
    prepare_initial_out = prepare_initial(h5ad_ch)

    // step 2
    prepare_generalized_out = prepare_generalized(h5ad_ch, prepare_initial_out.prepare_initial_done, prepare_initial_out.prepared_files)

    // step 3 create k and worker index combinations
    k_worker_combinations = Channel.from(kList)
        .combine(Channel.from(0..<params.total_workers))
        .view { k, worker -> ">>> Scheduling: k=${k}, worker=${worker}" }
    
    // 显示总任务数
    k_worker_combinations
        .count()
        .view { count -> ">>> Total parallel tasks: ${count}" }    
    
    // step 4 并行运行所有worker
    factorize_worker_out = factorize_worker(
        k_worker_combinations,
        h5ad_ch, 
        prepare_generalized_out.prepare_generalized_done
    )

    // step 5: 按k值分组，等待所有worker完成
    factorize_grouped = factorize_worker_out.factorize_worker_done
        .map { k, worker_idx, done_file -> [k, done_file] }  // 移除worker_idx
        .groupTuple()  // 按k值分组
    
    factorize_complete_out = factorize_complete(factorize_grouped)

    // 打印完成信息
    factorize_complete_out.factorize_complete.subscribe { k, done_file ->
        println ">>> All factorization completed for k=${k}"
    }

    // step 6: combine
    cnmf_combine_out = combine_results(
        factorize_complete_out.factorize_complete
            .map { k, done_file -> done_file }  // 只保留完成文件，移除k值
            .collect()  // 收集所有k值的完成标记文件
    )
    // 打印最终完成信息
    cnmf_combine_out.combine_done.subscribe { done_file ->
        println ">>> All CNMF analysis completed successfully!"
        println ">>> Results available in: ${params.outdir}/${params.name}/"
    }

    // step 7: kplot
    cnmf_kplot_out = kplot(cnmf_combine_out.combine_done)

}
