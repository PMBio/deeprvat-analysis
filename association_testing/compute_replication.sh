
python ${DEEPRVAT_ANALYSIS_DIR}/association_testing/compute_replication.py --phenotypes all --out-file replication/replication_all_phenotypes.parquet ./

python ${DEEPRVAT_ANALYSIS_DIR}/association_testing/compute_replication.py --phenotypes train --out-file replication/replication_training_phenotypes.parquet ./

python ${DEEPRVAT_ANALYSIS_DIR}/association_testing/compute_replication.py --phenotypes new --out-file replication/replication_new_phenotypes.parquet ./

