import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool

FLANK_FOLDER = "flank_regions"

def unwrap_gnrt_flank_regions(arg, **kwarg):
    return XReference.run_gnrt_flank_region_for_chrm(*arg, **kwarg)


class XReference():
    # def __init__(self, sf_ref, s_working_folder):
    #     self.sf_ref = sf_ref
    #     self.s_working_folder = s_working_folder
    #     if self.s_working_folder[-1] != "/":
    #         self.s_working_folder += "/"

    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "b_with_chr"
    def process_chrm_name(self, chrm, b_with_chr):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True

        # print chrm, self.b_with_chr, b_chrm_with_chr ###############################################################
        if b_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

    def load_in_sites(self, sf_sites):
        m_sites = {}
        with open(sf_sites) as fin_sites:
            for line in fin_sites:
                fields = line.split()
                chrm = fields[0]
                pos = int(fields[1])
                if chrm not in m_sites:
                    m_sites[chrm] = {}
                m_sites[chrm][pos] = 1
        return m_sites

    # gnrt the left and right flank regions for each candidate site
    def run_gnrt_flank_region_for_chrm(self, record):
        chrm = record[0]
        sf_sites = record[1]
        sf_ref = record[2]
        i_extend = int(record[3])
        working_folder = record[4]
        if working_folder[-1] != "/":
            working_folder += "/"

        m_sites = self.load_in_sites(sf_sites)
        f_fa = pysam.FastaFile(sf_ref)
        m_ref_chrms = {}
        for tmp_chrm in f_fa.references:
            m_ref_chrms[tmp_chrm] = 1
        b_with_chr = False
        if "chr1" in m_ref_chrms:
            b_with_chr = True

        for pos in m_sites[chrm]:
            if pos - i_extend < 0:
                continue
            istart = pos - i_extend
            iend = pos + i_extend

            ref_chrm = self.process_chrm_name(chrm, b_with_chr)
            s_left_region = f_fa.fetch(ref_chrm, istart, pos)
            s_right_region = f_fa.fetch(ref_chrm, pos + 1, iend)

            sf_flank_fa = working_folder + "{0}_{1}_flanks.fa".format(chrm, pos)
            with open(sf_flank_fa, "w") as fout_flank:
                fout_flank.write(">left\n")
                fout_flank.write(s_left_region + "\n")
                fout_flank.write(">right\n")
                fout_flank.write(s_right_region + "\n")
        f_fa.close()

    def gnrt_flank_region_for_sites(self, sf_sites, i_extend, n_jobs, sf_ref, s_working_folder):
        m_sites = self.load_in_sites(sf_sites)

        if s_working_folder[-1] != "/":
            s_working_folder += "/"
        flank_folder = s_working_folder + FLANK_FOLDER
        if os.path.exists(flank_folder) == False:
            cmd = "mkdir {0}".format(flank_folder)
            Popen(cmd, shell=True, stdout=PIPE).communicate()

        l_records = []
        for chrm in m_sites:
            l_records.append((chrm, sf_sites, sf_ref, i_extend, flank_folder))

        pool = Pool(n_jobs)
        pool.map(unwrap_gnrt_flank_regions, zip([self] * len(l_records), l_records), 1)
        pool.close()
        pool.join()

    #given chromosome lengths, break to bins
    def break_ref_to_bins(self, m_chrm_length, bin_size):
        m_bins={}
        for chrm in m_chrm_length:
            if chrm not in m_bins:
                m_bins[chrm]=[]
            chrm_lenth=m_chrm_length[chrm]
            tmp_pos=0
            while tmp_pos<chrm_lenth:
                block_end=tmp_pos+bin_size
                if block_end>chrm_lenth:
                    block_end=chrm_lenth
                m_bins[chrm].append((tmp_pos, block_end))
                tmp_pos=block_end
        return m_bins

    ##get index by chrm position
    def get_bin_by_pos(self, chrm, pos, bin_size, m_bins):
        bin_index=-1
        if chrm not in m_bins:
            return bin_index
        bin_index = pos/bin_size
        return m_bins[chrm][bin_index] #return [start, end) of the bin

    #