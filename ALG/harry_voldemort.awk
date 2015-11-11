#!/usr/bin/awk -f
BEGIN { har = 0;
	    vol = 0;
	    W = 0;
	  };
{W = W + NF} # For each line count the number of words
/Harry|Harry Potter|Potter/ {har++}  # add the number of times Harry appears
/Voldemort/ {vol++} # Add the number of times Voldemort appears
END {print har/W "\t" vol/W} # Return the result
