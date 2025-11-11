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
params.prepare_generalized_time = params.prepare_generalized_time ?: "2h"
params.factorize_worker_time = params.factorize_worker_time ?: "35h"
params.combine_time = params.combine_time ?: "2h"
params.kplot_time = params.kplot_time ?: "2h"

// computation memory
params.prepare_mem = params.prepare_mem ?: "10.G"
params.prepare_generalized_mem = params.prepare_generalized_mem ?: "2G"
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

// Get output directory as absolute path
def getOutDir() {
    return file(params.outdir).toAbsolutePath().toString()
}

// Log pipeline parameters
def logParameters(kList, outDir, params) {
    println """
    ╔════════════════════════════════════════════════════════════╗
    ║           cNMF Pipeline Configuration                      ║
    ╠════════════════════════════════════════════════════════════╣
    ║ Analysis Dir:    ${params.analysisDir}
    ║ Count File:      ${params.count_fp}
    ║ Output Dir:      ${outDir}
    ║ K Values:        ${kList.join(', ')}
    ║ N Iterations:    ${params.n_iter}
    ║ Num Genes:       ${params.num_genes}
    ║ Total Workers:   ${params.total_workers}
    ╚════════════════════════════════════════════════════════════╝
    """.stripIndent()
}


// Parse and log configuration
def kList = parseKValues(params.k_values)
def outDir = getOutDir()


/****************************
 *       PROCESSES
 ***************************/
// PROCESS A: Prepare initial data for CNMF
process prepare_initial {
    tag "prepare_initial"

    // Resource directives
    time params.prepare_time
    memory params.prepare_mem

    publishDir outDir, mode: 'copy', overwrite: true

    input:
    path h5ad_file

    output:
    path "prepare.done", emit: prepare_initial_done
    path "${params.name}/**", emit: prepared_files  // Output the whole directory structure

    script:
    """
    echo ">>> Running prepare_initial, output to: ${outDir}"

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

    publishDir outDir, mode: 'copy', overwrite: true // Ensure outputs are copied to the same output directory

    input:
    path h5ad_file
    path prepare_initial_done
    path prepared_files  // 接收prepare_initial的文件

    output:
    path "prepare_generalized.done", emit: prepare_generalized_done
    path "${params.name}/cnmf_tmp/**", emit: updated_files  // 输出更新后的目录结构


    script:
    """
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
    def outDirAbs = outDir
    def paramFile = "${outDirAbs}/${params.name}/cnmf_tmp/${params.name}.nmf_params.df.npz"
    def totalTasks = kList.size() * params.total_workers
    def currentTask = (kList.indexOf(k) * params.total_workers) + worker_idx + 1
    
    """
    echo ">>> Task ${currentTask} of ${totalTasks}: factorization k=${k}, worker=${worker_idx}"
    
    # Verify parameter file exists
    if [ ! -f "${paramFile}" ]; then
        echo ">>> ERROR: Parameter file not found: ${paramFile}"
        ls -la ${outDirAbs}/${params.name}/cnmf_tmp/ 2>/dev/null || echo ">>> Directory does not exist"
        exit 1
    fi

    # 直接使用publishDir目录，所有worker共享同一个目录结构
    cnmf factorize \
        --output-dir ${outDirAbs} \
        --name ${params.name} \
        --worker-index ${worker_idx} \
        --total-workers ${params.total_workers} \
        --skip-completed-runs

    # Verfiy output
    iter_count=\$(ls ${outDirAbs}/${params.name}/cnmf_tmp/*spectra.k_${k}.iter_*.df.npz 2>/dev/null | wc -l || echo 0)
    echo ">>> Worker ${worker_idx} completed for k=${k} (\${iter_count} iterations)"

    touch factorize_k${k}_worker${worker_idx}.done
    """
}

// PROCESS D: Collect results from all workers for a specific k value
process factorize_complete {
    tag "factorize_complete_k${k}"
    

    input:
    tuple val(k), path(worker_done_files)  // 收集所有worker的完成文件

    output:
    val(k), emit: k_completed

    script:
    def outDirAbs = outDir
    def output_dir = "${outDirAbs}/${params.name}/cnmf_tmp"
    def worker_count = worker_done_files.size()
    def expected_workers = params.total_workers

    """
    # Count the number of completed worker files
    echo ">>> Verifying factorization completion for k=\${k}"
    echo ">>> Expected workers: \${expected_workers}"
    echo ">>> Completed workers: \${worker_count}"

    # Verify that all workers have completed
    if [ "\${worker_count}" -ne "\${expected_workers}" ]; then
        echo ">>> ERROR: Not all workers completed for k=\${k}"
        echo ">>> Expected: \${expected_workers}, Got: \${worker_count}"
        exit 1
    fi
    
    spectra_count=\$(ls "${output_dir}"/*spectra.k_\${k}.iter_*.df.npz 2>/dev/null | wc -l || echo 0)
    if [ "\${spectra_count}" -eq 0 ]; then
        echo ">>> WARNING: No spectra files found for k=\${k}"
    else
        echo ">>> SUCCESS: Found \${spectra_count} spectra files for k=\${k}"
    fi
    
    echo ">>> All workers completed successfully for k=\${k}"
    """
}

// PROCESS E: cNMF combine results
process combine_results {
    tag "combine_results"
    
    time params.combine_time
    memory params.combine_mem

    publishDir outDir, mode: 'copy', overwrite: true

    input:
    val(k_list)  // 接收所有k值列表

    output:
    path("cnmf_combine.done"), emit: combine_done
    
    script:
    def outDirAbs = outDir
    """
    echo ">>> Running cNMF combine step"
    
    cnmf combine \
        --output-dir ${outDirAbs} \
        --name ${params.name}
    
    touch cnmf_combine.done
    """
}

// PROCESS F: k_selection plot
process kplot {
    tag "k_selection_plot"
    
    time params.kplot_time
    memory params.kplot_mem

    publishDir outDir, mode: 'copy', overwrite: true

    input:
    path(combine_done)  // 依赖combine步骤完成
    
    output:
    path("k_selection_plot.done"), emit: k_selection_plot_done
    
    script:
    def outDirAbs = outDir
    """
    echo ">>> Starting CNMF k selection plot generation"
    
    # 运行k selection plot命令
    cnmf k_selection_plot \
        --output-dir ${outDirAbs} \
        --name ${params.name}
    
    touch k_selection_plot.done
    """
}


/****************************
 *       WORKFLOWS
 ***************************/
workflow {

    logParameters(kList, outDir, params)

    h5ad_ch = Channel.fromPath(params.count_fp)
    
    // step 1
    prepare_initial_out = prepare_initial(h5ad_ch)

    // step 2
    prepare_generalized_out = prepare_generalized(
        h5ad_ch, 
        prepare_initial_out.prepare_initial_done, 
        prepare_initial_out.prepared_files
    )

    // step 3 create k and worker index combinations
   k_worker_combinations = Channel.from(kList)
        .combine(Channel.from(0..<params.total_workers))
    
    // 显示总任务数
    println """
    ╔═══════════════════════════════════════╗
    ║   Parallel Execution Plan             ║
    ╠═══════════════════════════════════════╣
    ║ K Values:        ${kList.join(', ')}
    ║ Total Workers:   ${params.total_workers}
    ║ Total Tasks:     ${kList.size() * params.total_workers}
    ╚═══════════════════════════════════════╝
    """.stripIndent() 
   
    // step 4 并行运行所有worker
    factorize_worker_out = factorize_worker(
        k_worker_combinations,
        h5ad_ch, 
        prepare_generalized_out.prepare_generalized_done
    )

    // step 5: 按k值分组，等待所有worker完成
    factorize_grouped = factorize_worker_out.factorize_worker_done
        .map { k, worker_idx, done_file -> [k, done_file] } 
        .groupTuple(by: 0)  // 按k值分组
    
    factorize_complete_out = factorize_complete(factorize_grouped)

    // 打印完成信息
    factorize_complete_out.k_completed.subscribe { k ->
        println ">>> All factorization completed for k=${k}"
    }

    // step 6: combine
    cnmf_combine_out = combine_results(
        factorize_complete_out.k_completed.collect()   // 收集所有k值成列表
    )

    // 打印最终完成信息
    cnmf_combine_out.combine_done.subscribe { done_file ->
        println ">>> All CNMF analyses completed successfully!"
        println ">>> Results available in: ${params.outdir}/${params.name}/"
    }

    // step 7: kplot
    cnmf_kplot_out = kplot(cnmf_combine_out.combine_done)

}
