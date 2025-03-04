awk '
#process fasta samples created by 'leave_one_out_sample_v3' script
#NOTE the value of n MUST be the same as in the sampling script (at the moment this only works for n=10 without adding extra conditionals to the write anyway)

#modify filter length depending on structure
#modify percent ambiguous positions
BEGIN{
filter_length=1271
filter_ambig=0.01
}

{   if ((sep=index($0,"-SPLIT-")) > 0) {
        header = substr($0,1,sep-1)
        sequence = substr($0, sep+7, length($0))
	if(length(sequence)<filter_length || gsub(/X/,"")>(length(sequence)*filter_ambig)){
          header=""
	  sequence=""
	}
	out_header[NR] = header
        out_sequence[NR] = sequence
      }
}

END{  for(i = 1; i<=NR; i++){
	if(out_header[i]!=""){
	   print out_header[i] > FILENAME
	   print out_sequence[i] > FILENAME
	}
    }
}
'  $*