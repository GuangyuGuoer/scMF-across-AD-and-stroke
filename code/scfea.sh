python /home/linux/GY/Fan/scFEA/src/scFEA.py --data_dir data --input_dir input \
                    --test_file ad_1212.csv \
                    --moduleGene_file module_gene_complete_mouse_m168.csv \
                    --stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv \
                    --cName_file cName_complete_mouse_c70_m168.csv \
                    --sc_imputation True \
                    --output_flux_file output/ad1212_flux.csv \
                    --output_balance_file output/ad1212_balance.csv
              

python /home/linux/GY/Fan/scFEA/src/scFEA.py  --data_dir data --input_dir input \
                    --test_file stroke_1212.csv \
                    --moduleGene_file module_gene_complete_mouse_m168.csv \
                    --stoichiometry_matrix cmMat_complete_mouse_c70_m168.csv \
		    --cName_file cName_complete_mouse_c70_m168.csv \
                    --sc_imputation True \
                    --output_flux_file output/stroke1212_flux.csv \
                    --output_balance_file output/stroke1212_balance.csv
