import os
from optparse import OptionParser
from x_clip_disc_filter import *


def parse_option():
    parser = OptionParser()
    parser.add_option("-i", "--input", dest="input",
                      help="input file ", metavar="FILE")
    parser.add_option("-r", "--reference", dest="reference",
                      help="The reference file ", metavar="FILE")
    parser.add_option("-c", "--cns", dest="cns",
                      help="The consensus file ", metavar="FILE")
    parser.add_option("-b", "--bam", dest="bam",
                      help="Input bam file", metavar="FILE")
    parser.add_option("-o", "--output", dest="output",
                      help="The output file", metavar="FILE")
    parser.add_option("-p", "--path", dest="wfolder", type="string",
                      help="Working folder")
    parser.add_option("-n", "--cores", dest="cores", type="int",
                      help="number of cores")
    parser.add_option("-e", "--extend", dest="extend", type="int",
                      help="extend length")

    (options, args) = parser.parse_args()
    return (options, args)

####
def get_candidate_list_from_TEA_mem_output(sf_input, sf_output):
    with open(sf_input) as fin_list, open(sf_output, "w") as fout_list:
        for line in fin_list:
            fields=line.split()
            if fields[0]=="sample":
                continue
            chrm=fields[1]
            clip_pos=fields[7]#by default, use the right clip position
            nlclip=int(fields[15])
            nrclip=int(fields[16])
            if nlclip>=nrclip:
                clip_pos=fields[6]
            fout_list.write(chrm+"\t"+clip_pos+"\n")


if __name__ == '__main__':
    (options, args) = parse_option()

    sf_bam_list = options.bam  ###read in a bam list file
    s_working_folder = options.wfolder
    n_jobs = options.cores
    sf_ref = options.reference  ###reference genome, some cram file require this file to open

    sf_candidate_list = options.input
    sf_rep_cns = options.cns  ####repeat copies here, with the flank regions
    sf_output = options.output

    iextnd = 300  ###for each site, re-collect reads in range [-iextnd, iextnd], this around ins +- 3*derivation
    bin_size = 50000000  # block size for parallelization
    i_flank_lenth = 500
    bmapped_cutoff = 0.5
    i_concord_dist = 550  # this should be the 3*std_derivation
    f_concord_ratio = 0.55
    n_clip_cutoff = 3
    n_disc_cutoff = 6


    sf_candidate_list_brkpnt=sf_candidate_list+".brkpnt"
    get_candidate_list_from_TEA_mem_output(sf_candidate_list, sf_candidate_list_brkpnt)

    x_cd_filter = XClipDiscFilter(sf_bam_list, s_working_folder, n_jobs, sf_ref)
    x_cd_filter.call_MEIs_consensus(sf_candidate_list_brkpnt, iextnd, bin_size, sf_rep_cns, "sf_flank", i_flank_lenth,
                                    bmapped_cutoff, i_concord_dist, f_concord_ratio, n_clip_cutoff, n_disc_cutoff,
                                    sf_output)

####zh