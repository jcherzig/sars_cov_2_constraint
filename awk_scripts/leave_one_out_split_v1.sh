awk '

# set field separator to null to split at blank lines
BEGIN{
RS=""
}

# print each sample to its own file
{print > ("loo_sample_" NR ".fasta")}

'  $*