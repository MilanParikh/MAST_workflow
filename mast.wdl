version 1.0

workflow mast {
    input {
    	String output_directory
        File anndata_file
        String celltype_col = 'broad_clusters'
        String sample_col = 'sample'
        #general parameters
        Int cpu = 16
        String memory = "64G"
        String docker = "mparikhbroad/mast_workflow:latest"
        Int preemptible = 2
    }

    String output_directory_stripped = sub(output_directory, "/+$", "")

    call generate_model {
        input:
            output_dir = output_directory_stripped,
            anndata_file = anndata_file,
            celltype_col = celltype_col,
            sample_col = sample_col,
            cpu=cpu,
            memory=memory,
            docker=docker,
            preemptible=preemptible
    }

    scatter(celltype in generate_model.celltypes_array) {
        call summarize_condition {
            input:
                output_dir = output_directory_stripped,
                zlmCond_file = generate_model.zlmCond_file,
                celltype = celltype,
                cpu=cpu,
                memory=memory,
                docker=docker,
                preemptible=preemptible
        }
    }
    
    output {
        Array[File] result_files = summarize_condition.result_file
    }
}

task generate_model {

    input {
        String output_dir
        File anndata_file
        String celltype_col
        String sample_col
        String memory
        Int cpu
        String docker
        Int preemptible
    }

    command <<<
        set -e

        R << CODE
        library(zellkonverter)
        library(MAST)

        sce <- readH5AD("~{anndata_file}")
        sca <- SceToSingleCellAssay(sce, class = "SingleCellAssay")

        gc()

        sca <- sca[freq(sca)>0.1,]

        cdr2 <- colSums(assay(sca)>0)
        colData(sca)\$ngeneson <- scale(cdr2)

        sampledata <- factor(colData(sca)\$~{sample_col})
        colData(sca)\$sampledata <- sampledata

        celltype <- factor(colData(sca)\$~{celltype_col})
        colData(sca)\$celltype <- celltype

        zlmCond <- zlm(formula = ~ngeneson + celltype + (1 | sampledata), 
               sca=sca, 
               method='glmer', 
               ebayes=F, 
               strictConvergence=F,
               fitArgsD=list(nAGQ = 0))

        saveRDS(zlmCond, 'zlmCond.Rds')

        fileConn<-file("celltypes.txt")
        writeLines(levels(celltype), fileConn)
        close(fileConn)
        CODE

        gsutil -m cp zlmCond.Rds ~{output_dir}/
        
    >>>

    output {
        File zlmCond_file = "zlmCond.Rds"
        Array[String] celltypes_array = read_lines("celltypes.txt")
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(anndata_file, "GB")*4) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}

task summarize_condition {

    input {
        String output_dir
        File zlmCond_file
        String celltype
        String memory
        Int cpu
        String docker
        Int preemptible
    }

    command {
        set -e

        R --no-save << CODE
        library(MAST)

        zlmCond <- readRDS('~{zlmCond_file}')

        summaryCond <- summary(zlmCond, doLRT='celltype~{celltype}')
        summaryDt <- summaryCond\$datatable
        result <- merge(summaryDt[contrast=='celltype~{celltype}' & component=='H',.(primerid, `Pr(\>Chisq)`)],
                        summaryDt[contrast=='celltype~{celltype}' & component=='logFC', .(primerid, coef)],
                        by='primerid') # logFC coefficients
        result[,coef:=result[,coef]/log(2)]
        result[,FDR:=p.adjust(`Pr(\>Chisq)`, 'fdr')]
        result = result[result\$FDR<0.01,, drop=F]

        result <- stats::na.omit(as.data.frame(result))

        saveRDS(summaryCond, '~{celltype}_summaryCond.Rds')
        write.csv(result, "~{celltype}.csv")

        gsutil -m cp ~{celltype}_summaryCond.Rds ~{output_dir}/
        gsutil -m cp ~{celltype}.csv ~{output_dir}/

        CODE
        
    }

    output {
        File summaryCond_file = "~{celltype}_summaryCond.Rds"
        File result_file = "~{celltype}.csv"
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(zlmCond_file, "GB")*4) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}