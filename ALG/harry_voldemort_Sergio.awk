#!/usr/bin/gawk -f

{ 
    BNAME = FILENAME;
    sub("-.+", "", BNAME);
    for (i = 1; i<=NF; i++) { 
        if ($i ~ /Harry/) { 
            harry[BNAME]++;
        }; 
        if ($i ~ /Voldemort/) { 
            voldemort[BNAME]++;
        }; 
        if ($i ~ /\w/) {
            words[BNAME]++
        }
    } 
} 
END{
    print "BOOK", "WORDS", "HARRY", "HARRY/WORDS", "VOLDEMORT", "VOLDEMORT/WORDS"; 
    for (book in voldemort) {
        print book, words[book], harry[book], harry[book]/words[book], voldemort[book], voldemort[book]/words[book];
    }
}