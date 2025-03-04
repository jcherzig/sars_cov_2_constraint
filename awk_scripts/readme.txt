These awk scripts are used to sample and process data prior to model training.

match_ids_date is used to generate subsamples for train/test sets based on the deposition dates of sequences from GISAID.

leave_one_out_sample_v3 is used to generate the LOO samples

leave_one_out_split_v1 splits the resulting file into individual samples

leave_one_out_postprocess_presplit_v1 cleans and processes the resulting samples ready for downstream applications