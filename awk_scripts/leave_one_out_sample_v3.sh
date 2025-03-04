awk '
#NOTE resulting files need processing with 'leave_one_out_postprocess' script to ensure correct fasta formatting. If you change the value of n you MUST also change it in the post processing script
#usage leave_one_out_sample_v3 FASTAFILE.fasta > OUTPUT.fasta (note above, output is not correctly formatted for fasta and needs post processing)

#initialize variables
# alter your value of n here; this is number of sequences left out of each sample and therefore also the resulting number of samples
BEGIN { n=10 
        line_num=1
}

# read in fasta file to array indexed by rownumber
{  sequence[NR] = $0
}

# first concatenate pairs of headers and sequences together with a separator we can later use to split
END{  for(i=1;i<=NR;i++){
          c_seq[i]= sequence[i] "-SPLIT-" sequence[++i]
         }

    # now iterate over all values up to n
      for(i=1;i<=n;i++){
          # reset line_num
          line_num=1
             for (l in c_seq){

                 # for values where i does not equal n we print all values not divisible by i
                 if (i != n){
                    if ((line_num % n) != i) {
                         print c_seq[l]
                    }
                    # add a line break after each full sample set
                    if (line_num == length(c_seq)) {
		     print ""
		    }
                    line_num++

                 # for values where i is equal to n we print out value divisible by n
                 } else {
                    if ((line_num % n) != 0) {
                         print c_seq[l]
                    }
                    # add a line break after each full sample set
	            if (line_num == length(c_seq)) {
		         print ""
		    }
                    line_num++
                 }
           }
     }
}

'  $*