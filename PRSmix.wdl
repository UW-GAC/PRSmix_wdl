version 1.0

workflow PRSmix {
    input {
        Boolean a_runStep1 = false
        Boolean a_runStep2 = false
        File? s1_weight_file
        File? score_inp
    }

    Int disk = 100
    Int ncores = 2
    Int memory = 8

    if (a_runStep1) {
        call s1_harmonize_SNPeffects {
            input: 
                disk = disk,
                ncores = ncores,
                memory = memory
        }
    }

    if (a_runStep2) {
        call s2_computePRS {
            input: 
                weight_file = select_first([s1_harmonize_SNPeffects.weight_out_file, s1_weight_file, ""]),
                disk = disk,
                ncores = ncores,
                memory = memory
        }
    }

    call s3_combine_PRS {
        input: 
            #original_beta_files_list = select_first([s1_harmonize_SNPeffects.weight_out_file, s1_weight_file, ""]),
            score_files_list = select_first([s2_computePRS.score_out, score_inp]),
            disk = disk,
            ncores = ncores,
            memory = memory
    }

    meta {
        author : "Buu Truong"
        email : "btruong@broadinstitute.org"
        description : "This workflow runs PRSmix - see the README on the github for more information - https://github.com/buutrg/PRSmix"
    }

    output {
        File? weight_out_file = s1_harmonize_SNPeffects.weight_out_file
        File? score_file = s2_computePRS.score_out
        File prsmix_output = s3_combine_PRS.prsmix_output
    }

}


task s1_harmonize_SNPeffects {
    input {
        File? s1_weight_file

        File? ref_file
        File? pgs_folder_tar 
        File? pgs_list
        Int snp_col = 1
        Int a1_col = 2
        Int beta_col = 3
        Boolean isheader = true
        Int chunk_size = 20
        String out = "harmonized_snpeff.txt"

        Int ncores
        Int disk
        Int memory
    }


    command <<<
        R --no-save << RSCRIPT

            library(PRSmix)
            library(data.table)
            options(datatable.fread.datatable=FALSE)

            print("s1_harmonize_SNPeffects")

            ref_file = "~{ref_file}"
            pgs_folder_tar = "~{pgs_folder_tar}"
            pgs_list = "~{pgs_list}"
            ncores = ~{ncores}
            snp_col = as.numeric(~{snp_col})
            a1_col = as.numeric(~{a1_col})
            beta_col = as.numeric(~{beta_col})
            isheader = as.logical("~{isheader}")
            chunk_size = as.numeric(~{chunk_size})
            out = "~{out}"
            s1_weight_file = "~{s1_weight_file}"
            ls()

            print(ref_file)
            print(pgs_folder_tar)

            system(paste0("tar -xvf ", pgs_folder_tar))
            system("ls")
            system("ls score_folder")
            pgs_folder = "score_folder"

            system(paste0("ls ", pgs_folder))

            harmonize_snpeffect_toALT(
                ref_file = ref_file, 
                pgs_folder = pgs_folder,
                pgs_list = pgs_list,
                snp_col = snp_col,
                a1_col = a1_col,
                beta_col = beta_col,
                isheader = isheader,
                chunk_size = chunk_size,
                ncores = ncores,
                out = out 
            )

        RSCRIPT

    >>>

    output {
        File weight_out_file = out
    }


    runtime {
        docker: "buutrg/prsmix:1.0.0"
        disks: "local-disk ${disk} HDD"
        memory: "${memory} GB"
        cpu : "${ncores}"
        bootDiskSizeGb: 100
    }

}


task s2_computePRS {
    input {
        File? score_inp

        File? bed 
        File? bim 
        File? fam 
        File weight_file
        String out = "out"

        Int disk = 30
        Int ncores = 1
        Int memory = 30
    }

    command <<<

        R --no-save << RSCRIPT

            library(PRSmix)
            library(data.table)
            options(datatable.fread.datatable=FALSE)

            print("s2_computePRS")

            bed = "~{bed}"
            geno = substring(bed, 1, nchar(bed) - 4)
            weight_file = "~{weight_file}"
            out = "~{out}"
            score_inp = "~{score_inp}"
            ls()

            compute_PRS(geno = geno, weight_file = weight_file, out = out, plink2_path="/usr/bin/plink2")

        RSCRIPT

    >>>

    output {
        File score_out = "${out}.sscore"
    }

    runtime {
        docker: "buutrg/prsmix:1.0.0"
        disks: "local-disk ${disk} HDD"
        memory: "${memory} GB"
        cpu : "${ncores}"
        bootDiskSizeGb: 100
    }

}


task s3_combine_PRS {
    input {
        File pheno_file
        File covariate_file
        File score_files_list
        File trait_specific_score_file
        String pheno_name
        Boolean isbinary
        String out = "out"
        Boolean liabilityR2 = false
        String IID_pheno = "IID"
        String covar_list
        String cat_covar_list = "none"
        Boolean is_extract_adjSNPeff = false
        File? original_beta_files_list
        String train_size_list = "300"
        String power_thres_list = "0.95"
        String pval_thres_list = "0.05"
        Boolean read_pred_training = false
        Boolean read_pred_testing = false

        Int ncores
        Int disk = 30
        Int memory = 30
    }

    command <<<

        R --no-save << RSCRIPT

            library(PRSmix)
            library(data.table)
            options(datatable.fread.datatable=FALSE)

            print("combine PRS")

            pheno_file = "~{pheno_file}"
            covariate_file = "~{covariate_file}"
            score_files_list = "~{score_files_list}"
            trait_specific_score_file = "~{trait_specific_score_file}"
            pheno_name = "~{pheno_name}"
            out = "~{out}"
            isbinary = as.logical("~{isbinary}")
            liabilityR2 = as.logical("~{liabilityR2}")
            IID_pheno = "~{IID_pheno}"
            covar_list = "~{covar_list}"
            cat_covar_list = "~{cat_covar_list}"
            ncores = as.numeric("~{ncores}")
            is_extract_adjSNPeff = as.logical("~{is_extract_adjSNPeff}")
            original_beta_files_list = ~{default="NULL" original_beta_files_list}
            train_size_list = "~{train_size_list}"
            power_thres_list = "~{power_thres_list}"
            pval_thres_list = "~{pval_thres_list}"
            read_pred_training = as.logical("~{read_pred_training}")
            read_pred_testing = as.logical("~{read_pred_testing}")
            ls()

            if (cat_covar_list == "none") { 
                cat_covar_list = NULL
            } else {
                cat_covar_list = unlist(strsplit(cat_covar_list, split=","))
            }

            covar_list = unlist(strsplit(covar_list, split=","))
            pval_thres_list = as.numeric(unlist(strsplit(pval_thres_list, split=",")))
            power_thres_list = as.numeric(unlist(strsplit(power_thres_list, split=",")))
            train_size_list = as.numeric(unlist(strsplit(train_size_list, split=",")))

            combine_PRS(
                pheno_file = pheno_file,
                covariate_file = covariate_file,
                score_files_list = score_files_list,
                trait_specific_score_file = trait_specific_score_file,
                pheno_name = pheno_name,
                out = out,
                isbinary = isbinary,
                liabilityR2 = liabilityR2,
                IID_pheno = IID_pheno,
                covar_list = covar_list,
                cat_covar_list = cat_covar_list,
                ncores = ncores,
                is_extract_adjSNPeff = is_extract_adjSNPeff,
                original_beta_files_list = original_beta_files_list,
                train_size_list = train_size_list,
                power_thres_list = power_thres_list,
                pval_thres_list = pval_thres_list,
                read_pred_training = read_pred_training, 
                read_pred_testing = read_pred_testing
                )

            system("ls")
            system("mkdir prsmix_output")

            system("mv *_train_allPRS.txt prsmix_output")
            system("mv *_test_allPRS.txt prsmix_output")

            system("mv *_train_df.txt prsmix_output")
            system("mv *_test_df.txt prsmix_output")

            system("mv *_test_summary_traitPRS_withPRSmix.txt prsmix_output")
            system("mv *_test_summary_traitPRS_withPRSmixPlus.txt prsmix_output")

            system("mv *_time_PRSmix.txt prsmix_output")
            system("mv *_time_PRSmixPlus.txt prsmix_output")

            system("tar -cvzf prsmix_output.tar prsmix_output")

        RSCRIPT
    >>>

    output {
        File prsmix_output = "prsmix_output.tar.gz"
    }

    runtime {
        docker: "buutrg/prsmix:1.0.0"
        disks: "local-disk ${disk} HDD"
        memory: "${memory} GB"
        cpu : "${ncores}"
        bootDiskSizeGb: 50
    }

}
