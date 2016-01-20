#!/usr/bin/python3
# encoding: utf-8

import operator


def get_sequence_from_FASTA_file(filename, residue, threshold=0.3):
    """returns an integer corresponding to the total number of protein
    sequences having a relative frequency higher
    or equal than a given threshold for a given residue.

    >>> print(get_sequence_from_FASTA_file("example_fasta_file.fasta.txt", "G", 0.1))
    120
    """
    result = 0
    residues = 0
    residue_c = 0
    with open(filename, "r") as file_h:
        for line in file_h:
            if line.startswith(">"):
                if residues != 0:
                    if residue_c / residues >= threshold:
                        result += 1
                residues = 0
                residue_c = 0
            else:
                residues += len(line.rstrip())
                residue_c += line.count(residue)
        else:
            if residues != 0:
                if residue_c / residues >= threshold:
                    result += 1
    return(result)


def print_sequence_tails(filename, first_n=10, last_m=10):
    """Print on the screen the protein
    identifier, the first N-aminoacids and the last M-aminoacids. The three
    fields must be separated by a tabulator, and one protein by line"""
    assert last_m > 0
    assert first_n > 0
    seq = ""  # To prevent error on the first iteration
    with open(filename, "r") as file_h:
        seq_id = ""
        for line in file_h:
            if line.startswith(">"):
                if seq != "":
                    print("{}\t{}\t{}".format(
                        seq_id, seq[0:first_n], seq[-last_m:]))
                seq_id = line[1:].rstrip()
                seq = ""
            else:
                seq += line.rstrip()


def calculate_aminoacid_frequencies(fasta_f, subsequences_f, output_f):
    """For all the subsequences on subsequence_f check how many sequences in
    the fasta_f file has each subsequence and which percentage it is.
    Return a summary with the number of proteins, the number of subsequences
    and for each subsequence the number of sequences the subsequences appears
    and the proportions."""

    def read_fasta(fasta_f):
        """Read a fasta file and return a dictionary with id as a key and
        sequence as a value"""
        file_fasta = {}
        with open(fasta_f) as fasta:
            for seq in fasta:
                if seq.startswith(">"):
                    id_seq = seq.rstrip()
                    file_fasta[id_seq] = ""
                else:
                    file_fasta[id_seq] += seq.rstrip()
        return file_fasta

    def read_subseq(subsequences_f):
        """Read a file with subsequences and yield them"""
        with open(subsequences_f) as subsequences:
            readed_subseq = [seq for seq in subsequences]
        return readed_subseq

    fasta_seq = read_fasta(fasta_f)
    subseqs = read_subseq(subsequences_f)
    output_fh = open(output_f, "w+")
    output_fh.write("#Number of proteins: {}\n".format(len(fasta_seq.items())))
    output_fh.write("#Number of subsequences: {}\n".format(len(subseqs)))
    output_fh.write("#subsequence proportions:\n")
    count_subseq = {}
    count_subseq.fromkeys(subseqs)
    for subseq in subseqs:
        count = 0
        for seq in fasta_seq.values():
            subseq = subseq.rstrip()
            if subseq in seq:
                count += 1
        else:
            count_subseq[subseq] = count
    else:
        pass
    sorted_x = sorted(count_subseq.items(),
                      key=operator.itemgetter(1),
                      reverse=True)
    for (val, count) in sorted_x:
        output_fh.write("{0}\t{1:>10}\t{2:.4f}\n".format(
                        val, count, count / len(subseqs)))
    output_fh.close()
