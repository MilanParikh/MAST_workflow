version 1.0

workflow mast {
    input {
    	String output_directory
        File anndata_file
        String celltype_col = 'broad_clusters'
        #general parameters
        Int cpu = 16
        String memory = "64G"
        String docker = "mparikhbroad/mast_workflow:latest"
        Int preemptible = 2
    }

    String output_directory_stripped = sub(output_directory, "/+$", "")

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
    
    output {
        File markers_file = run_MAST.markers_file
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