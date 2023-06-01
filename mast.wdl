version 1.0

workflow mast {
    input {
    	String output_directory
        File anndata_file
        String celltype_col = 'broad_clusters'
        Boolean parallel = false
        #general parameters
        Int cpu = 16
        String memory = "64G"
        String docker = "mparikhbroad/mast_workflow:latest"
        Int preemptible = 2
    }

    String output_directory_stripped = sub(output_directory, "/+$", "")

    if(parallel) {
        call convert_anndata {
            input:
                anndata_file = anndata_file,
                celltype_col = celltype_col,
                cpu=cpu,
                memory=memory,
                docker=docker,
                preemptible=preemptible
        }
        scatter(celltype in convert_anndata.celltypes_array) {
            call run_celltype_MAST {
                input:
                    output_dir = output_directory_stripped,
                    seurat_rds = convert_anndata.seurat_rds,
                    celltype = celltype,
                    cpu=cpu,
                    memory=memory,
                    docker=docker,
                    preemptible=preemptible
            }
        }
    }

    if(!parallel) {
        call run_MAST {
            input:
                output_dir = output_directory_stripped,
                anndata_file = anndata_file,
                celltype_col = celltype_col,
                cpu=cpu,
                memory=memory,
                docker=docker,
                preemptible=preemptible
        }
    }
    
    output {
        File? markers_file = run_MAST.markers_file
        Array[File]? celltype_markers_file = run_celltype_MAST.celltype_markers_file
    }
}

task convert_anndata {

    input {
        File anndata_file
        String celltype_col
        String memory
        Int cpu
        String docker
        Int preemptible
    }

    command <<<
        set -e

        R << CODE
        library(zellkonverter)
        library(Seurat)

        sce <- readH5AD("~{anndata_file}")
        seuratobj <- as.Seurat(sce, counts=NULL, data='X')
        Idents(object = seuratobj) <- seuratobj@meta.data\$~{celltype_col}
        saveRDS(seuratobj, 'seurat_obj.Rds')
        fileConn<-file("celltypes.txt")
        writeLines(levels(seuratobj@meta.data\$~{celltype_col}), fileConn)
        close(fileConn)
        CODE
        
    >>>

    output {
        File seurat_rds = "seurat_obj.Rds"
        Array[String] celltypes_array = read_lines("celltypes.txt")
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(anndata_file, "GB")*2) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}

task run_celltype_MAST {

    input {
        String output_dir
        File seurat_rds
        String celltype
        String memory
        Int cpu
        String docker
        Int preemptible
    }

    command <<<
        set -e

        R << CODE
        library(MAST)
        library(Seurat)

        seuratobj <- readRDS('~{seurat_rds}')
        markers <- FindAllMarkers(seuratobj, ident.1 = '~{celltype}', test.use="MAST")
        write.csv(markers, '~{celltype}_markers.csv')
        CODE

        gsutil -m cp ~{celltype}_markers.csv ~{output_dir}/
        
    >>>

    output {
        File celltype_markers_file = "~{celltype}_markers.csv"
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(seurat_rds, "GB")*2) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}

task run_MAST {

    input {
        String output_dir
        File anndata_file
        String celltype_col
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
        library(Seurat)

        sce <- readH5AD("~{anndata_file}")
        seuratobj <- as.Seurat(sce, counts=NULL, data='X')
        Idents(object = seuratobj) <- seuratobj@meta.data\$~{celltype_col}
        markers <- FindAllMarkers(seuratobj, test.use="MAST")
        write.csv(markers, 'markers.csv')
        CODE

        gsutil -m cp markers.csv ~{output_dir}/
        
    >>>

    output {
        File markers_file = "markers.csv"
    }

    runtime {
        docker: docker
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk " + ceil(size(anndata_file, "GB")*2) + " HDD"
        cpu: cpu
        preemptible: preemptible
    }

}