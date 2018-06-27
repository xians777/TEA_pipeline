##03/11/2018
##@@author: Simon (Chong) Chu, DBMI, Harvard Medical School
##@@contact: chong_chu@hms.harvard.edu

##The input sites are guranteed by the first two steps that:
# 1) Have enough clipped reads support (nearby region)
# 2) Several discordant reads support
# 3) Clipped reads and discordant read aligned to repeat copies

##This module mainly for filtering out false positive events by:
##1) Check the alignment position on consensus, whether form a cluster
##2) And the distance between the discordant reads aligned to consensus is within 3*std_dev of the clip_pos

##Expected output:
##1) Smaller number of candidate by filter out most of the false positives
##2) For each candidate site, output the TSD if exist
##3) Transduction analysis

##Steps needed:
# 1) For each input site, retrieve all the discordant and clipped reads in the focal region  (Need a seperate module?)
##This step will output:
# I. Clipped part in fastq format,
# II. Discordant positions on repeat consensus
# III. Peak left_clip_pos, and right_clip_pos, and candidate TSD
# 2) Re-align the clipped part to repeat copies
# This step will output:
# I. "Clip" positions on consensus repeats
# II. Check the consistency between peak "clip position" and discordant positions
# III. Call out candidate transduction? (optional)


###Hard code at:
#  def check_disc_consistency(self, m_clip_checked, m_disc_pos, m_polyA, i_dist, ratio):
# def call_MEIs(self, sf_candidate_list, extnd, bin_size, sf_rep_copies, i_flank_lenth, bmapped_cutoff,
#                  sf_annotation, i_concord_dist, f_concord_ratio, sf_final_list):
#

import os
import sys
import pysam
from subprocess import *
from multiprocessing import Pool
from x_alignments import *
from x_annotation import *
##from union_find_set import *

BWA_T_CUTOFF = 32  ##the minimum clipped length
BWA_REALIGN_CUTOFF = 9
NEARBY_REGION = 5
CLIP_FREQ = 10
TRIM_CLIP_FREQ = 2
MAX_CLIP_CLIP_LEN = 6
MINIMUM_DISC_MAPQ = 20
DISC_THRESHOLD = 2000
TSD_CUTOFF = 50
MINIMAL_TRANSDUCT_MAPQ = 30
N_MIN_A_T = 5  # minimu number of consecutive "A" or "T"
NEARBY_CLIP=50

BWA_PATH = "bwa"
SAMTOOLS_PATH = "samtools"
CLIP_FQ_SUFFIX = ".clipped.fq"
CLIP_BAM_SUFFIX = ".clipped.sam"
CLIP_POS_SUFFIX = ".clip_pos"
DISC_NAME_SUFFIX = ".disc_names"
DISC_POS_SUFFIX = ".disc_pos"
CLIP_RE_ALIGN_POS_SUFFIX = ".clip_realign_pos"
ALLELE_FREQUENCY_SUFFIX='.af'
OUTPUT_BAM_SUFFIX = ".out_bam"
OUTPUT_BAM_HEADER = ".bam_header.sam"
FLAG_LEFT_CLIP = "L"
FLAG_RIGHT_CLIP = "R"
SEPERATOR = '~'
####

def unwrap_self_collect_clip_disc_reads(arg, **kwarg):
    return XClipDisc.collect_clipped_disc_reads_by_region(*arg, **kwarg)

def unwrap_self_calc_AF_by_clip_reads(arg, **kwarg):
    return XClipDisc.calc_AF_of_site(*arg, **kwarg)

class XClipDisc():
    def __init__(self, sf_bam, working_folder, n_jobs, sf_ref):
        self.sf_bam = sf_bam
        self.working_folder = working_folder
        self.n_jobs = n_jobs
        self.sf_reference=sf_ref

    def get_realignment_for_sites(self):
        # Get the realignment position of the clipped parts from the output of the first step
        sf_bam_name = os.path.basename(self.sf_bam)
        sf_algnmt = self.working_folder + sf_bam_name + CLIP_BAM_SUFFIX

        ####given specific chrm, get all the related clipped reads

    ####
    ## "self.b_with_chr" is the format gotten from the alignment file
    ## all other format should be changed to consistent with the "self.b_with_chr"
    def _process_chrm_name(self, b_tmplt_with_chr, chrm):
        b_chrm_with_chr = False
        if len(chrm) > 3 and chrm[:3] == "chr":  ##Here remove the "chr"
            b_chrm_with_chr = True
        # print chrm, self.b_with_chr, b_chrm_with_chr #################################################################

        if b_tmplt_with_chr == True and b_chrm_with_chr == True:
            return chrm
        elif b_tmplt_with_chr == True and b_chrm_with_chr == False:
            return "chr" + chrm
        elif b_tmplt_with_chr == False and b_chrm_with_chr == True:
            return chrm[3:]
        else:
            return chrm

    ###Problem here: 1. chrm in "candidate_list" may not consistent with chrm in bam file
    ###2. all should follow the style in candidate list

    def collect_clipped_disc_reads_by_region(self, record):
        chrm = record[0][0]  ##this is the chrm style in candidate list
        insertion_pos = record[0][1]
        extnd = record[0][2]
        start_pos = insertion_pos - extnd
        if start_pos <= 0:
            start_pos = 1
        end_pos = insertion_pos + extnd
        sf_bam = record[1]
        working_folder = record[2]

        bam_info = BamInfo(sf_bam, self.sf_reference)
        b_with_chr = bam_info.is_chrm_contain_chr()
        chrm_in_bam = self._process_chrm_name(b_with_chr, chrm)

        # load the reads, and write the related clipped part into file
        s_pos_info = "{0}_{1}".format(chrm, insertion_pos)
        sf_clip_fq = working_folder + s_pos_info + CLIP_FQ_SUFFIX  # this is to save the clipped part for re-alignment
        f_clip_fq = open(sf_clip_fq, "w")
        # sf_clip_pos = working_folder + s_pos_info + CLIP_POS_SUFFIX  # this is to save clip positions
        # f_clip_pos = open(sf_clip_pos, "w")

        #sf_disc_names = working_folder + s_pos_info + DISC_NAME_SUFFIX
        #f_disc_names = open(sf_disc_names, "w")
        sf_disc_pos = working_folder + s_pos_info + DISC_POS_SUFFIX  # this is to save the discordant positions
        f_disc_pos = open(sf_disc_pos, "w")
        samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)
        n_cnt_clip = 0  # as a index to set a different read id for collected reads clip at the same position
        for algnmt in samfile.fetch(chrm_in_bam, start_pos, end_pos):  ##fetch reads mapped to "chrm:start_pos-end_pos"
            ##here need to skip the secondary and supplementary alignments?
            # if algnmt.is_secondary or algnmt.is_supplementary:
            #     continue
            if algnmt.is_duplicate == True:  ##duplciate
                continue
            b_first = True
            if algnmt.is_read2 == True:
                b_first = False
            if algnmt.is_unmapped == True:  #unmapped
                continue
            l_cigar = algnmt.cigar
            if len(l_cigar) < 1:  #wrong alignment
                continue
            if algnmt.mapping_quality < MINIMUM_DISC_MAPQ:
                continue

            # b_fully_mapped=False
            # if len(l_cigar) == 1 and l_cigar[0][0] == 0:  ##fully mapped
            #     b_fully_mapped=True

            query_name = algnmt.query_name
            query_seq = algnmt.query_sequence
            ##this is different from the one saved in the fastq/sam, no offset 33 to subtract
            query_quality = algnmt.query_qualities
            map_pos = algnmt.reference_start
            mate_chrm = '*'
            mate_pos = 0
            is_rc = 0
            if algnmt.is_reverse == True:  # is reverse complementary
                is_rc = 1
            if algnmt.mate_is_unmapped == False:
                mate_chrm = algnmt.next_reference_name
                mate_pos = algnmt.next_reference_start

            if l_cigar[0][0] == 4:  # left clipped
                if algnmt.is_supplementary or algnmt.is_secondary:  ###secondary and supplementary are not considered
                    continue

                # if map_pos in m_pos:
                clipped_seq = query_seq[:l_cigar[0][1]]
                if len(clipped_seq) < BWA_REALIGN_CUTOFF:
                    continue
                clipped_qulity = self._cvt_to_Ascii_quality(query_quality[:l_cigar[0][1]])
                clipped_rname = "{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}".format(chrm, SEPERATOR, map_pos, SEPERATOR,
                                                                            FLAG_LEFT_CLIP, SEPERATOR, is_rc, SEPERATOR,
                                                                            insertion_pos, SEPERATOR, n_cnt_clip)
                n_cnt_clip += 1

                f_clip_fq.write("@" + clipped_rname + "\n")
                f_clip_fq.write(clipped_seq + "\n")
                f_clip_fq.write("+\n")
                f_clip_fq.write(clipped_qulity + "\n")

            if l_cigar[-1][0] == 4:  # right clipped
                ##calculate the exact clip position
                for (type, lenth) in l_cigar[:-1]:
                    if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                        continue
                    else:
                        map_pos += lenth

                if algnmt.is_supplementary or algnmt.is_secondary:  ###secondary and supplementary are not considered
                    continue

                # if map_pos in m_pos:  #soft-clip
                clipped_rname = "{0}{1}{2}{3}{4}{5}{6}{7}{8}{9}{10}".format(chrm, SEPERATOR, map_pos, SEPERATOR,
                                                                            FLAG_RIGHT_CLIP, SEPERATOR, is_rc,
                                                                            SEPERATOR,
                                                                            insertion_pos, SEPERATOR, n_cnt_clip)
                n_cnt_clip += 1

                start_pos = -1 * l_cigar[-1][1]
                clipped_seq = query_seq[start_pos:]
                if len(clipped_seq) < BWA_REALIGN_CUTOFF:
                    continue
                clipped_qulity = self._cvt_to_Ascii_quality(query_quality[start_pos:])

                f_clip_fq.write("@" + clipped_rname + "\n")
                f_clip_fq.write(clipped_seq + "\n")
                f_clip_fq.write("+\n")
                f_clip_fq.write(clipped_qulity + "\n")

            if mate_chrm == "*":  ##unmapped reads are not interested!
                continue
            xfilter = XFilter()  ###decoy seuqence and contigs are not interested
            if xfilter.is_decoy_contig_chrms(mate_chrm) == True:
                continue

            ## here only collect the read names for discordant reads, later will re-align the discordant reads
            if self.is_discordant(chrm_in_bam, map_pos, mate_chrm, mate_pos, DISC_THRESHOLD) == True:
                # check where the mate is mapped, if within a repeat copy, then get the position on consensus
                #f_disc_names.write(query_name + "\n")
                s_mate_first = 1 #whether the mate read is the "first read" in a pair
                if b_first == True:
                    s_mate_first = 0

                # here mate_chrm must be the style in the bam file
                # And chrm must be the style in the candidate file
                s_mate_pos_info = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(chrm, map_pos, mate_chrm, mate_pos,
                                                                               query_name, s_mate_first, insertion_pos)
                f_disc_pos.write(s_mate_pos_info)

        samfile.close()
        f_clip_fq.close()
        f_disc_pos.close()


    ##whether a pair of read is discordant (for TEI only) or not
    def is_discordant(self, chrm, map_pos, mate_chrm, mate_pos, is_threshold):
        # b_disc=False
        if chrm != mate_chrm:  ###of different chroms
            return True
        else:
            if abs(mate_pos - map_pos) > is_threshold:  # Of same chrom, but insert size are quite large
                return True
                # if first_rc==second_rc: ##direction are abnormal, by default consider (F,R) as right mode
                #     return True
        return False

    ####This function: Given a list of insertion sites, get all the clipped parts (in fastq) and disc reads
    ##Parameters: "sf_candidate_list" saves all the candidate sites
    #             for each site, will collect reads within [-extnd, +extnd] region
    #             "sf_rep_copies" is the repeat copy files in fasta format (including flank regions)
    def collect_clipped_disc_reads_of_given_list(self, sf_candidate_list, extnd, bin_size, sf_all_clip_fq, sf_disc_fa):
        l_chrm_records = []
        with open(sf_candidate_list) as fin_list:
            for line in fin_list:
                fields = line.split()
                chrm = fields[0]
                pos = int(fields[1])  # candidate insertion site
                l_chrm_records.append(((chrm, pos, extnd), self.sf_bam, self.working_folder))
        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_collect_clip_disc_reads, zip([self] * len(l_chrm_records), l_chrm_records), 1)
        pool.close()
        pool.join()

        ###I. For clipped reads
        ## get all the clipped parts for each candidate site
        ###II. For disc reads
        ### here merge all the read_id files and pos files
        ####in format:chrm, map_pos, mate_chrm, mate_pos, query_name, s_mate_first, insertion_pos
        ####it is possible there are duplicate records are merged
        sf_all_disc_pos = self.working_folder + "all_disc_pos" + DISC_POS_SUFFIX #this is to save all the disc pos
        with open(sf_all_clip_fq, "w") as fout_merged_clip_fq, open(sf_all_disc_pos, "w") as fout_merged_disc_pos:
            #m_read_names={}####in case duplicate reads are saved
            for record in l_chrm_records:
                chrm = record[0][0]
                insertion_pos = record[0][1]
                s_pos_info = "{0}_{1}".format(chrm, insertion_pos)
                sf_clip_fq = self.working_folder + s_pos_info + CLIP_FQ_SUFFIX
                with open(sf_clip_fq) as fin_clip_fq:
                    for line in fin_clip_fq:
                        fout_merged_clip_fq.write(line)

                sf_disc_pos = self.working_folder + s_pos_info + DISC_POS_SUFFIX  # this is to save the disc positions
                with open(sf_disc_pos) as fin_disc_pos:
                    for line in fin_disc_pos:
                        fout_merged_disc_pos.write(line)

            # for read_name in m_read_names:
            #     fout_merged_disc_names.write(read_name)

        #now, need to retrieve the reads according to the read names and disc positions
        bam_info = BamInfo(self.sf_bam, self.sf_reference)
        bam_info.extract_mate_reads_by_name(sf_all_disc_pos, bin_size, self.working_folder, self.n_jobs, sf_disc_fa)

####
    '''
        for each alignment in the alignment_list:
    '''
    def _cvt_to_Ascii_quality(self, l_score):
        new_score = [x + 33 for x in l_score]
        return ''.join(map(chr, new_score))

    #calculate AF for one site
    #by default: cnt # of clipped reads within [-30, 30] of the insertion site
    ##calculate the average # of reads cover a site
    #then calculate the AF
    def calc_AF_of_site(self, record):
        chrm = record[0]  ##this is the chrm style in candidate list
        sf_bam_list=record[1]
        sf_candidate_list=record[2]
        extend=int(record[3])
        clip_extnd=int(record[4])
        af_cutoff=float(record[5])
        working_folder=record[6]
        if working_folder[-1]!="/":
            working_folder+="/"

        l_bams=[]
        with open(sf_bam_list) as fin_bam_list:
            for line in fin_bam_list:
                l_bams.append(line.rstrip())

        m_candidate_sites={}
        with open(sf_candidate_list) as fin_candidates:
            for line in fin_candidates:
                fields=line.split()
                tmp_chrm=fields[0]
                pos=int(fields[1])
                if tmp_chrm not in m_candidate_sites:
                    m_candidate_sites[tmp_chrm]=[]
                m_candidate_sites[tmp_chrm].append(pos)

        m_rslts={}
        for sf_bam in l_bams:
            bam_info = BamInfo(sf_bam, self.sf_reference)
            b_with_chr = bam_info.is_chrm_contain_chr()
            chrm_in_bam = self._process_chrm_name(b_with_chr, chrm)
            samfile = pysam.AlignmentFile(sf_bam, "rb", reference_filename=self.sf_reference)

            for insertion_pos in m_candidate_sites[chrm]:
                start_pos = insertion_pos - extend
                if start_pos <= 0:
                    start_pos = 1
                end_pos = insertion_pos + extend
                cnt_all_reads = 0
                cnt_clip = 0
                cnt_disc=0
                cnt_all_2=0
                lmost_region = -1
                rmost_region = -1
                for algnmt in samfile.fetch(chrm_in_bam, start_pos, end_pos):  ##fetch reads mapped to "chrm:start_pos-end_pos"
                    if algnmt.is_duplicate == True:  ##duplciate
                        continue
                    if algnmt.is_unmapped == True:  # unmapped
                        continue
                    l_cigar = algnmt.cigar
                    if len(l_cigar) < 1:  # wrong alignment
                        continue
                    if algnmt.mapping_quality < MINIMUM_DISC_MAPQ:
                        continue
                    if algnmt.is_supplementary or algnmt.is_secondary:
                        continue

                    map_pos = algnmt.reference_start
                    if lmost_region==-1:
                        lmost_region=map_pos
                    if rmost_region==-1:
                        rmost_region=map_pos

                    if map_pos<lmost_region:
                        lmost_region=map_pos
                    if map_pos>rmost_region:
                        rmost_region=map_pos

                    cnt_all_reads += 1
                    if abs(insertion_pos - map_pos) < extend:
                        cnt_all_2+=1

                    mate_chrm = '*'
                    mate_pos = 0
                    if algnmt.mate_is_unmapped == False:
                        mate_chrm = algnmt.next_reference_name
                        mate_pos = algnmt.next_reference_start

                    ##check disc information
                    if mate_chrm!="*":
                        if self.is_discordant(chrm_in_bam, map_pos, mate_chrm, mate_pos, DISC_THRESHOLD) == True:
                            cnt_disc+=1
                    ##check clip information
                    clip_pos=map_pos
                    if l_cigar[0][0] == 4:
                        if abs(insertion_pos-clip_pos)<=clip_extnd:
                            cnt_clip+=1
                    elif l_cigar[-1][0] == 4:
                        for (type, lenth) in l_cigar[:-1]:
                            if type == 4 or type == 5 or type == 1:  # (1 for insertion)
                                continue
                            else:
                                clip_pos += lenth
                        if abs(insertion_pos-clip_pos)<=clip_extnd:
                            cnt_clip+=1

                if insertion_pos not in m_rslts:
                    m_rslts[insertion_pos]=[]
                irange=abs(rmost_region-lmost_region)
                m_rslts[insertion_pos].append((cnt_clip, 2*clip_extnd, cnt_all_reads, irange, cnt_disc, cnt_all_2))
            samfile.close()
        sf_tmp_out=working_folder+chrm+ALLELE_FREQUENCY_SUFFIX
        with open(sf_tmp_out,"w") as fout_tmp:
            for ins_pos in m_rslts:
                cnt_clip=0
                cnt_all=0
                l_af=[]
                l_af_disc=[]
                for tmp_rcd in m_rslts[ins_pos]:
                    cnt_clip=tmp_rcd[0]
                    clip_range=tmp_rcd[1]
                    cnt_all=tmp_rcd[2]
                    all_range=tmp_rcd[3]
                    cnt_disc=tmp_rcd[4]
                    cnt_disc_all=tmp_rcd[5]

                    if cnt_all<=0:
                        l_af.append(str(0.0))
                        l_af_disc.append(str(0.0))
                    else:
                        af = float(cnt_clip * all_range) / float(cnt_all * clip_range)
                        if af>1.0:
                            af=1.0
                        l_af.append(str(af))

                        af_disc = float(cnt_disc) / float(cnt_disc_all)
                        if af_disc > 1.0:
                            af_disc = 1.0
                        l_af_disc.append(str(af_disc))

                s_af="\t".join(l_af)
                s_af_disc="\t".join(l_af_disc)

                b_clip_satisfield=False
                for clip_af in l_af:
                    if float(clip_af) <= af_cutoff:
                        b_clip_satisfield=True
                        break
                b_disc_satisfield=False
                for disc_af in l_af_disc:
                    if float(disc_af) <= af_cutoff:
                        b_disc_satisfield=True
                        break
                if b_clip_satisfield==True and b_disc_satisfield==True:
                    s_info="{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chrm, ins_pos, cnt_clip, cnt_all, s_af, s_af_disc)
                    fout_tmp.write(s_info)
####

    #Given:
    #   sf_bam_list: bam files
    #   sf_candidate_list: candidate sites
    #   #
    def calc_AF_by_clip_reads_of_given_list(self, sf_bam_list, sf_candidate_list, extnd, clip_slack, af, sf_out):
        m_chrm = {}
        with open(sf_candidate_list) as fin_list:
            for line in fin_list:
                fields = line.split()
                chrm = fields[0]
                m_chrm[chrm]=1

        l_chrm_records=[]
        for chrm in m_chrm:
            l_chrm_records.append((chrm, sf_bam_list, sf_candidate_list, extnd, clip_slack, af, self.working_folder))

        pool = Pool(self.n_jobs)
        pool.map(unwrap_self_calc_AF_by_clip_reads, zip([self] * len(l_chrm_records), l_chrm_records), 1)
        pool.close()
        pool.join()

        ####
        with open(sf_out, "w") as fout_rslt:
            for chrm in m_chrm:
                sf_tmp_out = self.working_folder + chrm + ALLELE_FREQUENCY_SUFFIX
                if os.path.isfile(sf_tmp_out)==False:
                    continue
                with open(sf_tmp_out) as fin_tmp:
                    for line in fin_tmp:
                        fout_rslt.write(line)


class XClipDiscFilter():
    def __init__(self, sf_bam_list, working_folder, n_jobs, sf_reference):
        self.sf_bam_list = sf_bam_list
        self.working_folder = working_folder
        self.n_jobs = n_jobs
        self.sf_reference=sf_reference

    def _get_n_bams(self):
        n_cnt=0
        with open(self.sf_bam_list) as fin_list:
            for line in fin_list:
                if len(line.rstrip())>0:
                    n_cnt+=1
        return n_cnt

####
    # Note: sf_rep_copies are copies with two end flank regions
    #this version align the reads to the consensus
    #then do the analysis
    def call_MEIs_consensus(self, sf_candidate_list, extnd, bin_size, sf_rep_cns, sf_flank, i_flank_lenth,
                            bmapped_cutoff, i_concord_dist, f_concord_ratio, nclip_cutoff, ndisc_cutoff, sf_final_list):
        # collect the clipped and discordant reads
        # each record in format like: @20~41951715~L~1~41952104~0~0,
        # (chrm, map_pos, FLAG_LEFT_CLIP, is_rc, insertion_pos, n_cnt_clip, sample_id)
        sf_clip_fq = self.working_folder + "candidate_sites_all_clip.fq"
        sf_disc_fa = self.working_folder + "candidate_sites_all_disc.fa"
        self.collect_clipped_disc_reads(sf_candidate_list, extnd, bin_size, sf_clip_fq, sf_disc_fa)

        # # ##re-align the clipped reads
        sf_clip_algnmt = self.working_folder + "temp_clip.sam"
        self.realign_clipped_reads(sf_rep_cns, sf_clip_fq, sf_clip_algnmt)
        # ##re-align the disc reads
        sf_disc_algnmt = self.working_folder + "temp_disc.sam"
        self.realign_disc_reads(sf_rep_cns, sf_disc_fa, sf_disc_algnmt)

        ####analysis the re-aligned clipped reads, called out:
        # 1) left-max-clip-position, right-max-clip-position, TSD
        # 2) peak clip position on the repeat consensus
        m_ins_lpos, m_ins_rpos, m_cns_lpos, m_cns_rpos, m_clip_sample, m_polyA = \
            self.parse_clip_realignment_consensus(sf_clip_algnmt, bmapped_cutoff)


        idist = 15  # count # of clipped reads within this range
        ratio = 0.55  # at lease "ratio" percent of the clipped reads are clipped within the region
        # Each returned record in format: (l_peak_pos, r_peak_pos, cns_peak_start, cns_peak_end)
        EXD_BIN_SIZE = idist
        n_bam = self.cnt_n_bams(self.sf_bam_list)

        #check clip consistency
        m_clip_checked_list = self.check_clip_alone_consistency2(m_ins_lpos, m_ins_rpos, m_polyA, idist,
                                                                ratio, nclip_cutoff, EXD_BIN_SIZE)

        m_qlfd_clip={}
        m_qlfd_clip_no_polyA={}
        m_potential_complex={}

        for ins_chrm in m_clip_checked_list:
            for ins_pos in m_clip_checked_list[ins_chrm]:
                b_polyA=m_clip_checked_list[ins_chrm][ins_pos][-1]
                b_consist=m_clip_checked_list[ins_chrm][ins_pos][-2]
                if b_consist==True and b_polyA==True:
                    if ins_chrm not in m_qlfd_clip:
                        m_qlfd_clip[ins_chrm]={}
                    m_qlfd_clip[ins_chrm][ins_pos]=m_clip_checked_list[ins_chrm][ins_pos]

                if b_consist==True:
                    if ins_chrm not in m_qlfd_clip_no_polyA:
                        m_qlfd_clip_no_polyA[ins_chrm]={}
                    m_qlfd_clip_no_polyA[ins_chrm][ins_pos]=m_clip_checked_list[ins_chrm][ins_pos]

                b_no_lclip=m_clip_checked_list[ins_chrm][ins_pos][-4]
                b_no_rclip=m_clip_checked_list[ins_chrm][ins_pos][-3]
                if (b_no_rclip==False and b_no_lclip==True) or (b_no_rclip==True and b_no_lclip==False):
                    if ins_chrm not in m_potential_complex:
                        m_potential_complex[ins_chrm]={}
                    m_potential_complex[ins_chrm][ins_pos]=m_clip_checked_list[ins_chrm][ins_pos]

        # 3) disc map position on the repeat consensus
        m_rep_pos_disc, m_disc_sample, m_disc_polyA = self.parse_disc_algnmt_consensus(sf_disc_algnmt, bmapped_cutoff)
        m_disc_checked_list = self.check_disc_consistency(m_clip_checked_list, m_rep_pos_disc, m_disc_sample,
                                                          i_concord_dist, f_concord_ratio, ndisc_cutoff)


        # ratio_clip_disc=0.25
        # #m_clip_disc_checked_list=self.check_clip_disc_consistency(m_disc_checked_list, m_rep_pos_disc, idist, ratio)
        # m_disc_checked_list2 = self.check_clip_disc_consistency(m_disc_checked_list, m_rep_pos_disc, idist, ratio_clip_disc)

        # m_disc_checked_list2={}
        # for ins_chrm in m_disc_checked_list:
        #     for ins_pos in m_disc_checked_list[ins_chrm]:
        #         if (ins_chrm in m_disc_polyA) and (ins_pos in m_disc_polyA[ins_chrm]) \
        #                 and (m_disc_polyA[ins_chrm][ins_pos]>1):######################################################
        #             if ins_chrm not in m_disc_checked_list2:
        #                 m_disc_checked_list2[ins_chrm]={}
        #             m_disc_checked_list2[ins_chrm][ins_pos]=m_disc_checked_list[ins_chrm][ins_pos]

        ###here let's combine the information and see
        # for ins_chrm in m_disc_checked_list:
        #     for ins_pos in m_disc_checked_list[ins_chrm]:
        #         if (ins_chrm in m_qlfd_clip_no_polyA) and (ins_pos in m_qlfd_clip_no_polyA[ins_chrm]):
        #             print "Mark_01:", ins_chrm, ins_pos
        #

        # # # 4) check consistency
        # # # first, check whether transduction, if so then report this
        # # return value in format [chr][pos][(pos-in-copy-with-flank, b_3mer_transduction)]
        n_cutoff = 3  ###################################number of clip and discordant reads fall in transduction region

        # # then, check consistence for disc and clip reads
        # return value in format [chr][(ins_pos, max_left_cnt, max_right_cnt, cnt_within_dist, cnt_all)]
        # m_candidates = self.check_clip_disc_consistency(m_clip_checked_list, m_rep_pos_disc, i_concord_dist,
        #                                                 f_concord_ratio)

        # #save the final candidates, and each in format:
        # #chrm, pos, TSD, TRN (position) , potential_length, genotype
        # #
        #m_pos_tsd = self.call_TSD_pos_gntp(m_qlfd_clip, TSD_CUTOFF)

        with open(sf_final_list, "w") as fout_final_list:
            for ins_chrm in m_disc_checked_list:
                for ins_pos in m_disc_checked_list[ins_chrm]:
                    record = m_disc_checked_list[ins_chrm][ins_pos]
                    lpeak_pos = record[0]
                    rpeak_pos = record[1]
                    TSD=0
                    if lpeak_pos>0 and rpeak_pos>0:
                        TSD=abs(rpeak_pos-lpeak_pos)

                    refined_pos=lpeak_pos
                    if lpeak_pos==-1 or lpeak_pos==None:
                        refined_pos=rpeak_pos
                    sinfo = "{0}\t{1}\t{2}\n".format(ins_chrm, refined_pos, TSD)
                    fout_final_list.write(sinfo)
########


    # Note: sf_rep_copies are copies with two end flank regions
    def call_MEIs(self, sf_candidate_list, extnd, bin_size, sf_rep_copies, i_flank_lenth, bmapped_cutoff,
                  sf_annotation, i_concord_dist, f_concord_ratio, sf_final_list):
        # collect the clipped and discordant reads
        #each record in format like: @20~41951715~L~1~41952104~0~0,
        # (chrm, map_pos, FLAG_LEFT_CLIP, is_rc, insertion_pos, n_cnt_clip, sample_id)
        sf_clip_fq = self.working_folder + "candidate_sites_all_clip.fq"
        sf_disc_fa = self.working_folder + "candidate_sites_all_disc.fa"
        self.collect_clipped_disc_reads(sf_candidate_list, extnd, bin_size, sf_clip_fq, sf_disc_fa)

        # # ##re-align the clipped reads
        sf_clip_algnmt = self.working_folder + "temp_clip.sam"
        self.realign_clipped_reads(sf_rep_copies, sf_clip_fq, sf_clip_algnmt)
        # ##re-align the disc reads
        sf_disc_algnmt = self.working_folder + "temp_disc.sam"
        self.realign_disc_reads(sf_rep_copies, sf_disc_fa, sf_disc_algnmt)


        ##check whether the chrm list contain "chr"
        xannotation = XAnnotation(sf_annotation)
        b_with_chr = self.is_sam_header_contain_chr(sf_clip_algnmt)
        xannotation.set_with_chr(b_with_chr)
        xannotation.load_rmsk_annotation_L1()  # load in the annotation
        m_rmsk_annotation = xannotation.get_rmsk_annotation()
        ####analysis the re-aligned clipped reads, called out:
        # 1) left-max-clip-position, right-max-clip-position, TSD
        # 2) peak clip position on the repeat consensus
        m_ins_lpos, m_ins_rpos, m_clip_sample, m_transduction_clip, m_polyA = self.parse_clip_realignment(
            sf_clip_algnmt, bmapped_cutoff, m_rmsk_annotation, i_flank_lenth)

####1. all the bams should have poly-A reads, total number should be larger than a threshold
####2. add those sites fall in the transduction (both clip and disc)

####for test only
        for i_chrm in m_clip_sample:
            for i_pos in m_clip_sample[i_chrm]:
                print "Mark0", i_chrm, i_pos
####
####
        idist = 15  # count # of clipped reads within this range
        ratio = 0.24  # at lease "ratio" percent of the clipped reads are clipped within the region
        # Each returned record in format: (l_peak_pos, r_peak_pos, cns_peak_start, cns_peak_end)
        EXD_BIN_SIZE = idist
        n_bam=self.cnt_n_bams(self.sf_bam_list)
        m_clip_checked_list = self.check_clip_consistency(m_ins_lpos, m_ins_rpos, m_clip_sample, m_polyA, n_bam, idist,
                                                          ratio, EXD_BIN_SIZE)
####for test only
        for i_chrm in m_clip_checked_list:
            for i_pos in m_clip_checked_list[i_chrm]:
                print "Mark1", i_chrm, m_clip_checked_list[i_chrm][i_pos][0]
            ####
            #
            ###save those cases that exist in m_transduction_clip, but missed by m_clip_checked_list
####
        # 3) disc map position on the repeat consensus
        m_rep_pos_disc, m_disc_sample, m_transduction_disc, m_disc_polyA = self.parse_disc_algnmt(sf_disc_algnmt,
                                                                                                  bmapped_cutoff,
                                                                                                  m_rmsk_annotation,
                                                                                                  i_flank_lenth)
        m_candidates = self.check_disc_consistency(m_clip_checked_list, m_rep_pos_disc, m_disc_sample,
                                                   i_concord_dist, f_concord_ratio)
####for test only
        for i_chrm in m_candidates:
            for i_pos in m_candidates[i_chrm]:
                print "Mark2", i_chrm, m_candidates[i_chrm][i_pos][0]
                # m_candidates = self.check_clip_disc_consistency(m_clip_checked_list, m_rep_pos_disc, i_concord_dist,
                #                                                 f_concord_ratio)
            ####
####
        # # # 4) check consistency
        # # # first, check whether transduction, if so then report this
        # # return value in format [chr][pos][(pos-in-copy-with-flank, b_3mer_transduction)]
        n_cutoff = 3  ###################################number of clip and discordant reads fall in transduction region
        m_transduct = self.check_transduction(m_transduction_clip, m_transduction_disc, i_flank_lenth, n_cutoff)
        #
        # # then, check consistence for disc and clip reads
        # return value in format [chr][(ins_pos, max_left_cnt, max_right_cnt, cnt_within_dist, cnt_all)]
        # m_candidates = self.check_clip_disc_consistency(m_clip_checked_list, m_rep_pos_disc, i_concord_dist,
        #                                                 f_concord_ratio)

        # #save the final candidates, and each in format:
        # #chrm, pos, TSD, TRN (position) , potential_length, genotype
        # #
        m_pos_tsd_transduct = self.call_TSD_pos_transduction_gntp(m_clip_checked_list, m_transduct, TSD_CUTOFF)

        #get all the transductions from both clip and disc reads
        n_clip_transduct_cutoff=3
        m_qlfd_transduct_clip=self.call_transduction_from_clip(m_transduction_clip, n_clip_transduct_cutoff)
        n_disc_transduct_cutoff=5
        m_qlfd_transduct_disc=self.call_transduction_from_disc(m_transduction_disc, n_disc_transduct_cutoff)

        with open(sf_final_list, "w") as fout_final_list:
            for ins_chrm in m_candidates:
                for ins_pos in m_candidates[ins_chrm]:
                    # print "Mark5: ", ins_chrm, ins_pos
                    if (ins_chrm in m_pos_tsd_transduct) and (ins_pos in m_pos_tsd_transduct[ins_chrm]):
                        refined_pos = m_pos_tsd_transduct[ins_chrm][ins_pos][0]  ###for now, this is the left pos
                        TSD = m_pos_tsd_transduct[ins_chrm][ins_pos][1]
                        # if TSD=="":####filter out those without target site duplicates
                        #     continue
                        s_transduct = "null"
                        if (ins_chrm in m_qlfd_transduct_disc) and (ins_pos in m_qlfd_transduct_disc[ins_chrm]):
                            s_transduct=m_qlfd_transduct_disc[ins_chrm][ins_pos]
                        # cnt_within_dist = m_candidates[ins_chrm][ins_pos][0]
                        # cnt_all = m_candidates[ins_chrm][ins_pos][1]
                        # sinfo = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(ins_chrm, refined_pos, TSD, s_transduct,
                        #                                                 cnt_within_dist, cnt_all)
                        sinfo = "{0}\t{1}\t{2}\t{3}\n".format(ins_chrm, refined_pos, TSD, s_transduct)
                        fout_final_list.write(sinfo)
            # #here also chech those doesn't fall into m_candidates events
            # for ins_chrm in m_qlfd_transduct_disc:
            #     for ins_pos in m_qlfd_transduct_disc[ins_chrm]:
            #         if (ins_chrm in m_candidates) and (ins_pos in m_candidates[ins_chrm]):
            #             continue
            #         s_transduct=m_qlfd_transduct_disc[ins_chrm][ins_pos]
            #         sinfo = "{0}\t{1}\t{2}\t{3}\n".format(ins_chrm, ins_pos, "null", s_transduct)
            #         fout_final_list.write(sinfo)
####

    ####
    def cnt_n_bams(self, sf_bam):
        n_cnt=0
        with open(sf_bam) as fin_bam:
            for line in fin_bam:
                sinfo=line.rstrip()
                if len(sinfo)>1:
                    n_cnt+=1
        return n_cnt

    ####
    def is_sam_header_contain_chr(self, sf_sam):
        samfile = pysam.AlignmentFile(sf_sam, "r", reference_filename=self.sf_reference)
        header = samfile.header
        samfile.close()
        l_chrms = header['SQ']
        m_chrms = {}
        for record in l_chrms:
            sinfo = record['SN']
            fields = sinfo.split(SEPERATOR)
            chrm_name = fields[0]
            m_chrms[chrm_name] = 1
        if ("1" in m_chrms) or ("2" in m_chrms):
            return False
        else:
            return True

    ####find out candidate transduction events from "clip transduct" events
    #m_clip_transduct in format:
    #m_transduction[ori_chrm][ori_insertion_pos][transduct_pos]=[list of positions
    #!!!!should also consider whether clip within a region
    #also one insertion may have multiple candidate transduction positions, so find the maximum one
    def call_transduction_from_clip(self, m_clip_transduct, n_cutoff):
        m_transduct_candidates={}
        for ins_chrm in m_clip_transduct:
            for ins_pos in m_clip_transduct[ins_chrm]:
                s_picked_pos=""
                n_max=0
                for transduct_pos in m_clip_transduct[ins_chrm][ins_pos]:
                    n_support_read=len(m_clip_transduct[ins_chrm][ins_pos][transduct_pos])
                    if n_support_read>n_max:
                        n_max=n_support_read
                        s_picked_pos=transduct_pos
                if n_max>n_cutoff:
                    if ins_chrm not in m_transduct_candidates:
                        m_transduct_candidates[ins_chrm]={}
                    m_transduct_candidates[ins_chrm][ins_pos]=s_picked_pos
        return m_transduct_candidates


    ####Given the transduction candidate list from clipped reads, check the disc candidates
    ####Note that: some candidate may only exist in disc candidates, so need to add them to the finaly list
    ####This is because the "disc transduct" are gotten from the "clip list", which may not cover "clip transduct"
    #!!!!!!should also consider whether form cluster
    # also one insertion may have multiple candidate transduction positions, so find the maximum one
    def call_transduction_from_disc(self, m_disc_transduct, n_cutoff):
        m_transduct_candidates = {}
        for ins_chrm in m_disc_transduct:
            for ins_pos in m_disc_transduct[ins_chrm]:
                s_picked_pos = ""
                n_max = 0
                for trsdct_pos in m_disc_transduct[ins_chrm][ins_pos]:
                    n_support_read = len(m_disc_transduct[ins_chrm][ins_pos][trsdct_pos])
                    if n_support_read>n_max:
                        n_max=n_support_read
                        s_picked_pos=trsdct_pos
                if n_max>n_cutoff:
                    if ins_chrm not in m_transduct_candidates:
                        m_transduct_candidates[ins_chrm]={}
                    m_transduct_candidates[ins_chrm][ins_pos]=s_picked_pos
        return m_transduct_candidates

    # m_ins_lpos in format: m_ins_lpos[ori_chrm][ori_insertion_pos][ori_mpos]=cnt
    def call_TSD_pos_transduction_gntp(self, m_clip_checked_list, m_transduct, TSD_cutoff):
        # first get the peak left clip position
        # then get the peak right clip position
        # call out the TSD, if the len<=25, also need to get the "seq" of the TSD
        # check whether it is transduction, if yes, report the linked copy position
        # call out the genotype information

        # The candidate sites will have both left and right clipped reads in a region
        # so we can just traverse the m_ins_lpos
        m_pos_TSD_transduct = {}
        for ins_chrm in m_clip_checked_list:
            for ins_pos in m_clip_checked_list[ins_chrm]:
                record = m_clip_checked_list[ins_chrm][ins_pos]
                lpeak_pos = record[0]
                rpeak_pos = record[1]

                TSD = ""
                if abs(rpeak_pos - lpeak_pos) <= TSD_cutoff:
                    TSD = "{0}~{1}~{2}".format(ins_chrm, lpeak_pos, rpeak_pos)

                s_transduct = ""
                if (ins_chrm in m_transduct) and (ins_pos in m_transduct[ins_chrm]):
                    s_transduct = "{0}~{1}".format(m_transduct[ins_chrm][ins_pos][0], m_transduct[ins_chrm][ins_pos][1])

                if ins_chrm not in m_pos_TSD_transduct:
                    m_pos_TSD_transduct[ins_chrm] = {}
                if ins_pos not in m_pos_TSD_transduct[ins_chrm]:
                    m_pos_TSD_transduct[ins_chrm][ins_pos] = (lpeak_pos, TSD, s_transduct)
        return m_pos_TSD_transduct


    # m_ins_lpos in format: m_ins_lpos[ori_chrm][ori_insertion_pos][ori_mpos]=cnt
    def call_TSD_pos_gntp(self, m_clip_checked_list, TSD_cutoff):
        # first get the peak left clip position
        # then get the peak right clip position
        # call out the TSD, if the len<=25, also need to get the "seq" of the TSD
        # check whether it is transduction, if yes, report the linked copy position
        # call out the genotype information

        # The candidate sites will have both left and right clipped reads in a region
        # so we can just traverse the m_ins_lpos
        m_pos_TSD_transduct = {}
        for ins_chrm in m_clip_checked_list:
            for ins_pos in m_clip_checked_list[ins_chrm]:
                record = m_clip_checked_list[ins_chrm][ins_pos]
                lpeak_pos = record[0]
                rpeak_pos = record[1]

                TSD = ""
                if abs(rpeak_pos - lpeak_pos) <= TSD_cutoff:
                    TSD = "{0}~{1}~{2}".format(ins_chrm, lpeak_pos, rpeak_pos)

                if ins_chrm not in m_pos_TSD_transduct:
                    m_pos_TSD_transduct[ins_chrm] = {}
                if ins_pos not in m_pos_TSD_transduct[ins_chrm]:
                    m_pos_TSD_transduct[ins_chrm][ins_pos] = (lpeak_pos, TSD)
        return m_pos_TSD_transduct

    ####
    ####left and right, both should be considered!!!


    ####
    def merge_sites(self, m_1, m_2):
        m_check_sites = {}
        for ins_chrm in m_1:
            if ins_chrm not in m_check_sites:
                m_check_sites[ins_chrm] = {}
            for ins_pos in m_1[ins_chrm]:
                if ins_pos not in m_check_sites[ins_chrm]:
                    m_check_sites[ins_chrm][ins_pos] = 1
        for ins_chrm in m_2:
            if ins_chrm not in m_check_sites:
                m_check_sites[ins_chrm] = {}
            for ins_pos in m_2[ins_chrm]:
                if ins_pos not in m_check_sites[ins_chrm]:
                    m_check_sites[ins_chrm][ins_pos] = 1
        return m_check_sites

    # check whether tranduction exists
    def check_transduction(self, m_clip_transduct, m_disc_transduct, flank_length, n_cutoff):
        m_check_sites = self.merge_sites(m_clip_transduct, m_disc_transduct)
        m_candidate_transduct = {}
        for ins_chrm in m_check_sites:
            for ins_pos in m_check_sites[ins_chrm]:
                cnt_clip_support = 0
                cnt_disc_support = 0
                cnt_max_support = 0
                max_copy_pos = ""

                if (ins_chrm in m_clip_transduct) and (ins_pos in m_clip_transduct[ins_chrm]):
                    for rep_copy_pos in m_clip_transduct[ins_chrm][ins_pos]:
                        cnt_clip_support = len(m_clip_transduct[ins_chrm][ins_pos][rep_copy_pos])
                        if (ins_chrm in m_disc_transduct) and (ins_pos in m_disc_transduct[ins_chrm]) and (
                                    rep_copy_pos in m_disc_transduct[ins_chrm][ins_pos]):
                            cnt_disc_support = len(m_disc_transduct[ins_chrm][ins_pos][rep_copy_pos])
                        tmp_all = cnt_clip_support + cnt_disc_support
                        if cnt_max_support < tmp_all:
                            cnt_max_support = tmp_all
                            max_copy_pos = rep_copy_pos
                elif (ins_chrm in m_disc_transduct) and (ins_pos in m_disc_transduct[ins_chrm]):
                    for rep_copy_pos in m_disc_transduct[ins_chrm][ins_pos]:
                        cnt_disc_support = len(m_disc_transduct[ins_chrm][ins_pos][rep_copy_pos])
                        if cnt_max_support < cnt_disc_support:
                            cnt_max_support = cnt_disc_support
                            max_copy_pos = rep_copy_pos

                # find out it's 3' or 5' transduction
                # first find the peak is at first or end
                # fields=max_copy_pos.split(SEPERATOR)
                # copy_length=int(fields[-1]) - int(fields[-2])
                cnt_3 = 0
                cnt_5 = 0
                if (ins_chrm in m_clip_transduct) and (ins_pos in m_clip_transduct[ins_chrm]) and (
                            max_copy_pos in m_clip_transduct[ins_chrm][ins_pos]):
                    for map_pos_flank in m_clip_transduct[ins_chrm][ins_pos][max_copy_pos]:
                        if map_pos_flank < flank_length:
                            cnt_5 += 1
                        else:
                            cnt_3 += 1
                if (ins_chrm in m_disc_transduct) and (ins_pos in m_disc_transduct[ins_chrm]) and (
                            max_copy_pos in m_disc_transduct[ins_chrm][ins_pos]):
                    for map_pos_flank in m_disc_transduct[ins_chrm][ins_pos][max_copy_pos]:
                        if map_pos_flank < flank_length:
                            cnt_5 += 1
                        else:
                            cnt_3 += 1

                b_3mer = 1
                if cnt_5 >= cnt_3:
                    b_3mer = 0

                b_transduct = False
                if (b_3mer == True and cnt_3 > n_cutoff) or (b_3mer == False and cnt_5 > n_cutoff):
                    b_transduct = True
                # here also need to get the positions on the consensus, to get the estimated length of the insertion
                if ins_chrm not in m_candidate_transduct:
                    m_candidate_transduct[ins_chrm] = {}

                if b_transduct == True:
                    m_candidate_transduct[ins_chrm][ins_pos] = (max_copy_pos, b_3mer)
        return m_candidate_transduct

    ##find the first and second peak clip position
    # ins_chrm is the ori reported insertion chrm on the reference
    # ins_pos is the ori reported insertion pos on the reference
    def find_first_second_peak(self, m_clip_pos, ins_chrm, ins_pos, BIN_SIZE):
        ##first check whether m_clip_pos[ins_chrm][ins_pos] exist
        peak_pos = -1
        peak_acm = 0
        scnd_peak_pos = -1
        scnd_peak_acm = 0
        mpos = m_clip_pos[ins_chrm][ins_pos]  # is a dictionary {ref_pos:[(cns_pos)]}
        istart = -1 * (BIN_SIZE / 2)
        iend = BIN_SIZE / 2
        for pos in mpos:  # cns_pos may be -1
            if pos <= 0:
                continue
            cur_acm = 0
            for offset in range(istart, iend):
                tmp_pos = pos + offset
                if tmp_pos in mpos:
                    cur_acm += len(mpos[tmp_pos]) ###number of reads clipped at the ref position
            if cur_acm > peak_acm:
                scnd_peak_pos = peak_pos
                scnd_peak_acm = peak_acm
                peak_acm = cur_acm
                peak_pos = pos
            elif cur_acm > scnd_peak_acm:
                scnd_peak_acm = cur_acm
                scnd_peak_pos = pos

        l_peak_nbr = []
        l_scnd_peak_nbr = []
        for offset in range(istart, iend):
            tmp_pos = peak_pos + offset
            if tmp_pos in mpos:
                l_peak_nbr.append(tmp_pos)

            tmp_scnd_pos = scnd_peak_pos + offset
            if tmp_scnd_pos in mpos:
                l_scnd_peak_nbr.append(tmp_scnd_pos)

        return peak_pos, l_peak_nbr, peak_acm, scnd_peak_pos, l_scnd_peak_nbr, scnd_peak_acm

    # find the best combination of left and right peak
    def find_best_combination(self, l1, l2, l1a, l2a, r1, r2, r1a, r2a, TSD_cutoff):
        b_TSD = True
        if abs(l1 - r1) <= TSD_cutoff:
            return l1, r1, b_TSD, l1a, r1a
        else:
            # for combination of l1 and r2
            if abs(l1 - r2) <= TSD_cutoff and abs(l2 - r1) > TSD_cutoff:
                return l1, r2, b_TSD, l1a, r2a
            elif abs(l2 - r1) <= TSD_cutoff and abs(l1 - r2) > TSD_cutoff:
                return l2, r1, b_TSD, l2a, r1a
            elif abs(l1 - r2) <= TSD_cutoff and abs(
                            l2 - r1) <= TSD_cutoff:  # both satisfied, then check accumulated values
                if (l1a + r2a) > (l2a + r1a):
                    return l1, r2, b_TSD, l1a, r2a
                else:
                    return l2, r1, b_TSD, l2a, r1a
            else:  # none of the combination is satisfied, then check poly-A
                b_TSD = False
                return l1, r1, b_TSD, l1a, r1a
                #####

    ##find the window with the maximum of clip position
    # ins_chrm is the ori reported insertion chrm on the reference
    # ins_pos is the ori reported insertion pos on the reference
    def find_peak_window(self, m_clip_pos, ins_chrm, ins_pos, BIN_SIZE):
        ##first check whether m_clip_pos[ins_chrm][ins_pos] exist
        peak_pos = -1
        peak_acm = 0
        mpos = m_clip_pos[ins_chrm][ins_pos]  # is a dictionary {ref_pos:[cns_pos]}
        istart = -1 * BIN_SIZE / 2
        iend = BIN_SIZE / 2
        for pos in mpos:  # cns_pos may be -1
            if pos <= 0:
                continue
            cur_acm = 0
            for offset in range(istart, iend):
                tmp_pos = pos + offset
                if tmp_pos in mpos:
                    cur_acm += len(mpos[tmp_pos])
            if cur_acm > peak_acm:
                peak_acm = cur_acm
                peak_pos = pos
        l_peak_nbr = []
        for offset in range(istart, iend):
            tmp_pos = peak_pos + offset
            if tmp_pos in mpos:
                l_peak_nbr.append(tmp_pos)
        return peak_pos, l_peak_nbr, peak_acm

    ####
    # check whether the clipped parts aligned on the repeat copies are clustered in a region
    # return True/False, and the peak window
    def check_clip_rep_consistency(self, m_ref_clip_pos, ins_chrm, ins_pos, l_ngbr_pos, imax_dist, ratio_cutoff):
        m_cns_pos = {}
        # 1. collect all the positions
        for map_pos_ref in l_ngbr_pos:  # this is the map position on the reference
            for cns_pos in m_ref_clip_pos[ins_chrm][ins_pos][map_pos_ref]:
                if cns_pos not in m_cns_pos:
                    m_cns_pos[cns_pos] = 1
                else:
                    m_cns_pos[cns_pos] += 1

        peak_win_start = 0
        peak_win_end = 0
        if len(m_cns_pos) == 0:
            return False, peak_win_start, peak_win_end

        # 2. find the peak window of m_cns_pos
        l_pos = list(m_cns_pos.keys())
        l_pos.sort()  ###sort the candidate sites
        peak_win_cnt = 0
        cur_win_cnt = 0
        cur_win_start = l_pos[0]
        cur_win_end = 0
        pre_pos = 0
        for cur_pos in l_pos:
            if (cur_pos - pre_pos) > imax_dist and pre_pos != 0:
                if cur_win_cnt > peak_win_cnt:
                    peak_win_cnt = cur_win_cnt
                    peak_win_start = cur_win_start
                    peak_win_end = cur_win_end

                cur_win_start = cur_pos
                cur_win_cnt = 0

            cur_win_end = cur_pos
            cur_win_cnt += m_cns_pos[cur_pos]
            pre_pos = cur_pos

        ##for the last comparison
        if cur_win_cnt > peak_win_cnt:
            peak_win_cnt = cur_win_cnt
            peak_win_start = cur_win_start
            peak_win_end = cur_win_end

        # 3. check the consistency, the ratio of reads fall in the peak window
        n_hit = 0
        for pos in l_pos:
            if pos >= peak_win_start and pos <= peak_win_end:
                n_hit += 1

        if (float(n_hit) / float(len(l_pos))) >= ratio_cutoff:
            return True, peak_win_start, peak_win_end
        else:
            return False, peak_win_start, peak_win_end
####
####
    def check_clip_sample_consistency(self, m_sample, ins_chrm, ins_pos, n_bam, ncutoff):
        if len(m_sample[ins_chrm][ins_pos])<n_bam:
            #print "AAAA", len(m_sample[ins_chrm][ins_pos]), n_bam
            return False
        n_qualified_sample = 0
        # require all the samples must have certain number of clipped reads
        for sample_id in m_sample[ins_chrm][ins_pos]:
            n_left=m_sample[ins_chrm][ins_pos][sample_id][0]
            n_right=m_sample[ins_chrm][ins_pos][sample_id][1]
            if (n_left+n_right) >= ncutoff:
                n_qualified_sample+=1
        if n_qualified_sample<n_bam:
            #print "BBB", n_qualified_sample, n_bam
            return False
        return True
####
####
    ####All the clipped position in the repeat consensus should focus on specific positions
    # this is to check the clipped parts of those reads clipped at left/right peak postions,
    # whether they also clustered at a focal region on the repeat consensus
    # also, we need to seperate the poly-A cluster and the other parts
    # Note: m_left_clip_pos in format {chrm:{pos:{map_pos:[cns_pos]}}}
    # nlpolyA and nrpolyA are the number of polyA clipped reads, by comparison, we can know which part is the end
    # Here, we also require one end is poly-A
    def check_clip_consistency(self, m_left_clip_pos, m_right_clip_pos, m_sample, m_polyA, n_bam, idist, ratio, BSIZE):
        m_checked_list = {}

        # 1. check the polyA and position
        for ins_chrm in m_left_clip_pos:
            for ins_pos in m_left_clip_pos[ins_chrm]:
                ##Here if more than one bam file provided, then require all of them have clip reads support
                # 1.1. decide which side is poly-A
                b_left_polyA = True
                nlpolyA = 0
                nrpolyA = 0
                if (ins_chrm in m_polyA) and (ins_pos in m_polyA[ins_chrm]):
                    nlpolyA = m_polyA[ins_chrm][ins_pos][0]
                    nrpolyA = m_polyA[ins_chrm][ins_pos][1]
                if nrpolyA > nlpolyA:
                    b_left_polyA = False

                b_polyA = False  # also check whether one end is poly-A
                if nlpolyA > 0 or nrpolyA > 0:
                    b_polyA = True

                l1_pos, l1_nbr, l1_acm, l2_pos, l2_nbr, l2_acm = self.find_first_second_peak(m_left_clip_pos, ins_chrm,
                                                                                             ins_pos, BSIZE)
                if (ins_chrm not in m_right_clip_pos) or (ins_pos not in m_right_clip_pos[ins_chrm]):
                    print "Position: {0}:{1} do not have right clipped positions".format(ins_chrm, ins_pos)
                    continue
                r1_pos, r1_nbr, r1_acm, r2_pos, r2_nbr, r2_acm = self.find_first_second_peak(m_right_clip_pos, ins_chrm,
                                                                                             ins_pos, BSIZE)
                l_peak_pos, r_peak_pos, b_TSD, lacm, racm = self.find_best_combination(l1_pos, l2_pos, l1_acm, l2_acm,
                                                                                       r1_pos, r2_pos, r1_acm, r2_acm,
                                                                                       TSD_CUTOFF)

                l_peak_nbr = None
                r_peak_nbr = None
                if l_peak_pos == l2_pos:
                    l_peak_nbr = l2_nbr
                else:
                    l_peak_nbr = l1_nbr
                if r_peak_pos == r2_pos:
                    r_peak_nbr = r2_nbr
                else:
                    r_peak_nbr = r1_nbr

                # 1.2 check whether the clip reads aligned on consensus are clustered
                b_consist = False
                cns_peak_start = 0
                cns_peak_end = 0

                # make sure the peak window is at the end of consensus
                # the right part form a peak cluster
                b_lconsist, cns_lpeak_start, cns_lpeak_end = self.check_clip_rep_consistency(m_left_clip_pos, ins_chrm,
                                                                                             ins_pos, l_peak_nbr,
                                                                                             idist, ratio)
                b_rconsist, cns_rpeak_start, cns_rpeak_end = self.check_clip_rep_consistency(m_right_clip_pos,
                                                                                             ins_chrm,
                                                                                             ins_pos, r_peak_nbr,
                                                                                             idist, ratio)
                if b_lconsist == True or b_rconsist == True:
                    b_consist = True
####for test!!!!
                if b_consist==True:
                    if lacm>racm:
                        print "Mark10", ins_chrm, l_peak_pos
                    else:
                        print "Mark10", ins_chrm, r_peak_pos

                # in some TP cases, there is no polyA happen, like 6:162181268 in NA12878 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                if b_polyA == False:  ####Here also require one end is poly-A
                    b_consist = False

                # if b_polyA==True:
                #     if b_left_polyA==True:
                #         cns_peak_start=cns_rpeak_start
                #         cns_peak_end=cns_rpeak_end
                #     else:
                #         cns_peak_start = cns_lpeak_start
                #         cns_peak_end = cns_lpeak_end

                if ins_chrm not in m_checked_list:
                    m_checked_list[ins_chrm] = {}
                if b_consist == True:
                    m_checked_list[ins_chrm][ins_pos] = (
                        l_peak_pos, r_peak_pos, cns_lpeak_start, cns_lpeak_end, cns_rpeak_start, cns_rpeak_end)

        return m_checked_list
#

    #this is to check one side of the clipped reads, whether they form a cluster on the consensus repeat
    def check_one_side_clip_consistency(self, m_clip_pos, idist, ratio, BSIZE):
        m_checked_list = {}
        for ins_chrm in m_clip_pos:
            for ins_pos in m_clip_pos[ins_chrm]:
                l1_pos, l1_nbr, l1_acm, l2_pos, l2_nbr, l2_acm = self.find_first_second_peak(m_clip_pos,
                                                                                             ins_chrm, ins_pos, BSIZE)

                # 1.2 check whether the clip reads aligned on consensus are clustered
                # make sure the peak window is at the end of consensus
                # the right part form a peak cluster
                b_lconsist, cns_lpeak_start, cns_lpeak_end = self.check_clip_rep_consistency(m_clip_pos,
                                                                                             ins_chrm, ins_pos,
                                                                                             l1_nbr, idist, ratio)

                if ins_chrm not in m_checked_list:
                    m_checked_list[ins_chrm] = {}
                m_checked_list[ins_chrm][ins_pos] = (l1_pos, cns_lpeak_start, cns_lpeak_end, b_lconsist, l1_acm)
        return m_checked_list

####
    ###this version consider the left and right clip positions seperatelly, thus we can include some complex events:
    #like TE insertion with deletion
    ##Also check TSD consistency:
        #if candidate site has both left and right clipped reads
        #then the left and right peak clip position should satisfy the TSD constrain
        #otherwise, if only have left or right clip, then no need to check this TSD constrain
        #for ins_chrm in m_clip_checked_list
    def check_clip_alone_consistency2(self, m_left_clip_pos, m_right_clip_pos, m_polyA, idist, ratio, n_cutoff, BSIZE):
        m_checked_list={}
        m_left_checked=self.check_one_side_clip_consistency(m_left_clip_pos, idist, ratio, BSIZE)
        m_right_checked=self.check_one_side_clip_consistency(m_right_clip_pos, idist, ratio, BSIZE)


        for ins_chrm in m_left_checked:
            for ins_pos in m_left_checked[ins_chrm]:#this is for events have both left and right clipped reads
                # if (ins_chrm not in m_polyA) or (ins_pos not in m_polyA[ins_chrm]): ##################################
                #     continue
                lrecord=m_left_checked[ins_chrm][ins_pos]#(l1_pos, cns_lpeak_start, cns_lpeak_end, b_lconsist, n_clip)
                cnt_clip=lrecord[-1]
                if (ins_chrm in m_right_checked) and (ins_pos in m_right_checked[ins_chrm]):
                    rrecord=m_right_checked[ins_chrm][ins_pos]#(l1_pos, cns_peak_start, cns_peak_end, b_lconsist, nclip)
                    ##here make sure the TSD is consistency
                    # if abs(lrecord[0] - rrecord[0]) > TSD_CUTOFF:####
                    #     #print "MMMM", ins_chrm, ins_pos, lrecord[0], rrecord[0] #####################################
                    #     continue
                    cnt_clip += rrecord[-1]
                    if cnt_clip<n_cutoff:
                        continue
                    if ins_chrm not in m_checked_list:
                        m_checked_list[ins_chrm]={}

                    m_checked_list[ins_chrm][ins_pos]=(lrecord[0], rrecord[0], lrecord[1], lrecord[2], rrecord[1],
                                                       rrecord[2], False, False, True, False)
                else:#this is the only have left-clipped
                    if cnt_clip<n_cutoff:
                        continue
                    if ins_chrm not in m_checked_list:
                        m_checked_list[ins_chrm]={}
                    m_checked_list[ins_chrm][ins_pos] = (lrecord[0], None, lrecord[1], lrecord[2], -1,
                                                         -1, False, True, True, False)

        #this is to focus on the events only have right-clipped
        for ins_chrm in m_right_checked:
            for ins_pos in m_right_checked[ins_chrm]:
                # if (ins_chrm not in m_polyA) or (ins_pos not in m_polyA[ins_chrm]): ##################################
                #     continue
                rrecord=m_right_checked[ins_chrm][ins_pos]
                cnt_clip=rrecord[-1]
                if cnt_clip<n_cutoff:
                    continue
                if (ins_chrm in m_left_checked) and (ins_pos in m_left_checked[ins_chrm]):
                    continue
                if ins_chrm not in m_checked_list:
                    m_checked_list[ins_chrm]={}
                m_checked_list[ins_chrm][ins_pos]=(None, rrecord[0], -1, -1, rrecord[1], rrecord[2],
                                                   True, False, True, False)

        ####
        return m_checked_list
#                 m_checked_list[ins_chrm][ins_pos] = (l_peak_pos, r_peak_pos, cns_lpeak_start, cns_lpeak_end,
#                                                      cns_rpeak_start, cns_rpeak_end, b_no_lclip, b_no_rclip,
#                                                      b_consist, b_polyA)
####
#
#
    ####All the clipped position in the repeat consensus should focus on specific positions
    # this is to check the clipped parts of those reads clipped at left/right peak postions,
    # whether they also clustered at a focal region on the repeat consensus
    # also, we need to seperate the poly-A cluster and the other parts
    # Note: m_left_clip_pos in format {chrm:{pos:{map_pos:[cns_pos]}}}
    # nlpolyA and nrpolyA are the number of polyA clipped reads, by comparison, we can know which part is the end
    # Here, we also require one end is poly-A
    def check_clip_alone_consistency(self, m_left_clip_pos, m_right_clip_pos, m_polyA, idist, ratio, BSIZE):
        m_checked_list = {}
        # 1. check the polyA and position
        for ins_chrm in m_left_clip_pos:
            for ins_pos in m_left_clip_pos[ins_chrm]:
                ##Here if more than one bam file provided, then require all of them have clip reads support
                # 1.1. decide which side is poly-A
                b_left_polyA = True
                nlpolyA = 0
                nrpolyA = 0
                if (ins_chrm in m_polyA) and (ins_pos in m_polyA[ins_chrm]):
                    nlpolyA = m_polyA[ins_chrm][ins_pos][0]
                    nrpolyA = m_polyA[ins_chrm][ins_pos][1]
                if nrpolyA > nlpolyA:
                    b_left_polyA = False

                b_polyA = False  # also check whether one end is poly-A
                if nlpolyA > 0 or nrpolyA > 0:
                    b_polyA = True

                ####Note that: for complex events, like deletion with a insertion, it is possible only one side of clip
                b_no_lclip=False
                b_no_rclip=False
                if (ins_chrm not in m_left_clip_pos) or (ins_pos not in m_left_clip_pos[ins_chrm]):
                    print "Position: {0}:{1} do not have left clipped positions".format(ins_chrm, ins_pos)
                    b_no_lclip=True
                l1_pos, l1_nbr, l1_acm, l2_pos, l2_nbr, l2_acm = self.find_first_second_peak(m_left_clip_pos,
                                                                                             ins_chrm, ins_pos, BSIZE)
                if (ins_chrm not in m_right_clip_pos) or (ins_pos not in m_right_clip_pos[ins_chrm]):
                    print "Position: {0}:{1} do not have right clipped positions".format(ins_chrm, ins_pos)
                    b_no_rclip=True
 #######problem here!!!!!!!!!!!!!!!!!!!!!
                    continue

                r1_pos, r1_nbr, r1_acm, r2_pos, r2_nbr, r2_acm = self.find_first_second_peak(m_right_clip_pos,
                                                                                             ins_chrm, ins_pos, BSIZE)
                l_peak_pos, r_peak_pos, b_TSD, lacm, racm = self.find_best_combination(l1_pos, l2_pos, l1_acm,
                                                                                       l2_acm, r1_pos, r2_pos, r1_acm,
                                                                                       r2_acm, TSD_CUTOFF)

                l_peak_nbr = None
                r_peak_nbr = None
                if l_peak_pos == l2_pos:
                    l_peak_nbr = l2_nbr
                else:
                    l_peak_nbr = l1_nbr
                if r_peak_pos == r2_pos:
                    r_peak_nbr = r2_nbr
                else:
                    r_peak_nbr = r1_nbr

                # 1.2 check whether the clip reads aligned on consensus are clustered
                # make sure the peak window is at the end of consensus
                # the right part form a peak cluster
                b_lconsist, cns_lpeak_start, cns_lpeak_end = self.check_clip_rep_consistency(m_left_clip_pos,
                                                                                             ins_chrm, ins_pos,
                                                                                             l_peak_nbr, idist, ratio)
                b_rconsist, cns_rpeak_start, cns_rpeak_end = self.check_clip_rep_consistency(m_right_clip_pos,
                                                                                             ins_chrm, ins_pos,
                                                                                             r_peak_nbr, idist, ratio)
                b_consist = False
                if b_lconsist == True or b_rconsist == True:
                    b_consist = True

                # cns_peak_start = 0
                # cns_peak_end = 0
                # in some TP cases, there is no polyA happen, like 6:162181268 in NA12878 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                # if b_polyA == False:  ####Here also require one end is poly-A
                #     b_consist = False

                # if b_polyA==True:
                #     if b_left_polyA==True:
                #         cns_peak_start=cns_rpeak_start
                #         cns_peak_end=cns_rpeak_end
                #     else:
                #         cns_peak_start = cns_lpeak_start
                #         cns_peak_end = cns_lpeak_end

                if ins_chrm not in m_checked_list:
                    m_checked_list[ins_chrm] = {}

                m_checked_list[ins_chrm][ins_pos] = (l_peak_pos, r_peak_pos, cns_lpeak_start, cns_lpeak_end,
                                                     cns_rpeak_start, cns_rpeak_end, b_no_lclip, b_no_rclip,
                                                     b_consist, b_polyA)
        return m_checked_list
####

    def find_max_cover(self, l_pos, win_size):
        i_right = 0
        i_total = len(l_pos)
        i_cnt = 0
        i_max = 0
        i_start = 0
        i_end = 0
        for i_left in range(i_total):
            cur_val = l_pos[i_left]
            if i_right < i_left:
                i_right = i_left
            for i_temp in range(i_right, i_total):
                right_val = l_pos[i_temp]
                i_right = i_temp
                if (right_val - cur_val) > win_size:
                    if i_max < i_cnt:
                        i_max = i_cnt
                        i_start = i_left
                        i_end = i_right - 1
                    break
                i_cnt += 1
            if i_max < i_cnt:  # for the last one
                i_max = i_cnt
                i_start = i_left
                i_end = i_right

            i_cnt -= 1
            if i_cnt < 0:
                i_cnt = 0
            if i_right == (i_total - 1):
                break
        return i_start, i_end

    ####
    ####check whether most of the disc reads (> ratio) will form cluster and located within a focal region
    def check_disc_consistency(self, m_clip_checked, m_disc_pos, m_sample_disc, i_dist, ratio, n_disc_cutoff):
        m_selected = {}
        n_bams_all=self._get_n_bams() #all the bams
        for ins_chrm in m_clip_checked:
            for ins_pos in m_clip_checked[ins_chrm]:
                # require all the bams should at least have some discordant reads
                #n_cutoff_disc_read = 4  ####At least 2 discordant reads#########################Hard code !!!!!!!!!!!!!!!!!
                b_sample_consist = True
                tmp_n_bams=0
                if (ins_chrm in m_sample_disc) and (ins_pos in m_sample_disc[ins_chrm]):
                    for sample_id in m_sample_disc[ins_chrm][ins_pos]:
                        tmp_n_bams+=1
                        if m_sample_disc[ins_chrm][ins_pos][sample_id] < n_disc_cutoff:#all the bams have support
                            b_sample_consist = False
                            break

                if b_sample_consist == False:
                    continue

                if tmp_n_bams<n_bams_all:#require discordant reads exist at all the bams
                    continue

                if (ins_chrm in m_disc_pos) and (ins_pos in m_disc_pos[ins_chrm]):
                    # 2. find the peak window of m_cns_pos
                    l_pos = m_disc_pos[ins_chrm][ins_pos]
                    l_pos.sort()  ###sort the candidate sites
                    n_all_cnt = len(l_pos)
                    i_max_left, i_max_right = self.find_max_cover(l_pos, i_dist)
                    peak_win_cnt = i_max_right - i_max_left

                    l_sup_pos = l_pos[:i_max_left] + l_pos[i_max_right + 1:]
                    l_sup_max_left, l_sub_max_right = self.find_max_cover(l_sup_pos, i_dist)
                    second_peak_win_cnt = l_sub_max_right - l_sup_max_left

                    n_hit = peak_win_cnt + second_peak_win_cnt
                    if (float(n_hit) / float(n_all_cnt)) > ratio:
                        if ins_chrm not in m_selected:
                            m_selected[ins_chrm] = {}
                        if ins_pos not in m_selected[ins_chrm]:
                            m_selected[ins_chrm][ins_pos] = m_clip_checked[ins_chrm][ins_pos]
####
        return m_selected
####
####
####
    ####check whether most of the disc reads (> ratio) will form cluster and located within a focal region
    def check_disc_alone_consistency(self, m_disc_pos, m_sample_disc, i_dist, ratio):
        m_selected = {}
        n_bams_all=self._get_n_bams() #all the bams
        for ins_chrm in m_disc_pos:
            for ins_pos in m_disc_pos[ins_chrm]:
                # require all the bams should at least have some discordant reads
                n_cutoff_disc_read = 2  #####################################################Hard code !!!!!!!!!!!!!!!!!
                b_sample_consist = True
                tmp_n_bams=0
                if (ins_chrm in m_sample_disc) and (ins_pos in m_sample_disc[ins_chrm]):
                    for sample_id in m_sample_disc[ins_chrm][ins_pos]:
                        tmp_n_bams+=1
                        if m_sample_disc[ins_chrm][ins_pos][sample_id] < n_cutoff_disc_read:#all the bams have support
                            b_sample_consist = False
                            break
                if b_sample_consist == False:
                    continue
                if tmp_n_bams<n_bams_all:#require discordant reads exist at all the bams
                    continue

                # 2. find the peak window of m_cns_pos
                l_pos = m_disc_pos[ins_chrm][ins_pos]
                l_pos.sort()  ###sort the candidate sites
                n_all_cnt = len(l_pos)
                i_max_left, i_max_right = self.find_max_cover(l_pos, i_dist)
                peak_win_cnt = i_max_right - i_max_left

                l_sup_pos = l_pos[:i_max_left] + l_pos[i_max_right + 1:]
                l_sup_max_left, l_sub_max_right = self.find_max_cover(l_sup_pos, i_dist)
                second_peak_win_cnt = l_sub_max_right - l_sup_max_left

                n_hit = peak_win_cnt + second_peak_win_cnt
                if (float(n_hit) / float(n_all_cnt)) > ratio:
                    if ins_chrm not in m_selected:
                        m_selected[ins_chrm] = {}
                    if ins_pos not in m_selected[ins_chrm]:
                        m_selected[ins_chrm][ins_pos] = 1
        return m_selected
####

    # check whehter the the discordant position are within the some distance from the clip position
    # Input: 1) clip positions on the consensus repeat, left and right peak are the boundary of the insertion
    #        2) disc positions on the consensus repeat
    #        3) maximum distance between clip position and disc position
    #        4) minimum ratio of qualified disc
    def check_clip_disc_consistency(self, m_clip_checked, m_disc_pos, idist, ratio):
        m_candidates = {}
        ###for each candidate, there must be both left and right clipped reads in the nearby region,
        ###so, here we can only traverse the left clipped pos
        for ins_chrm in m_clip_checked:
            for ins_pos in m_clip_checked[ins_chrm]:
                # check discord reads mapped position whether consistent with clip positions
                record = m_clip_checked[ins_chrm][ins_pos]
                max_lpos = record[0]
                max_rpos = record[1]
                cns_lpeak_start = record[2]
                cns_lpeak_end = record[3]
                cns_rpeak_start = record[4]
                cns_rpeak_end = record[5]
                #
                b_consistent = False
                cnt_all = 0
                cnt_within_dist = 0
                if (ins_chrm in m_disc_pos) and (ins_pos in m_disc_pos[ins_chrm]):
                    cnt_all = len(m_disc_pos[ins_chrm][ins_pos])
                    for disc_pos in m_disc_pos[ins_chrm][ins_pos]:
                        # check whether disc_pos is within [cns_peak_start, cns_peak_end] or not
                        lwin_left_boundary = cns_lpeak_start - idist
                        lwin_right_boundary = cns_lpeak_end + idist
                        rwin_left_boundary = cns_rpeak_start - idist
                        rwin_right_boundary = cns_rpeak_end + idist
                        if lwin_left_boundary>0 and disc_pos >= lwin_left_boundary and disc_pos <= lwin_right_boundary:
                            cnt_within_dist += 1
                        elif rwin_left_boundary>0 and disc_pos >= rwin_left_boundary and disc_pos <= rwin_right_boundary:
                            cnt_within_dist += 1
                    if cnt_all > 0 and (float(cnt_within_dist) / float(cnt_all)) >= ratio:
                        b_consistent = True

                if b_consistent == True:
                    if ins_chrm not in m_candidates:
                        m_candidates[ins_chrm] = {}
                    m_candidates[ins_chrm][ins_pos]=m_clip_checked[ins_chrm][ins_pos]
                    #m_candidates[ins_chrm][ins_pos] = (cnt_within_dist, cnt_all)
        return m_candidates


    # check the clipped part is qualified aligned or not
    def is_clipped_part_qualified_algnmt(self, l_cigar, ratio_cutoff):
        if len(l_cigar) < 1:  # wrong alignment
            return False, 0
        if len(l_cigar) > 2:
            ####check the cigar
            ###if both clipped, and the clipped part is large, then skip
            b_left_clip = False
            i_left_clip_len = 0
            if l_cigar[0][0] == 4 or l_cigar[0][0] == 5:  # left clipped
                b_left_clip = True
                i_left_clip_len = l_cigar[0][1]
            b_right_clip = False
            i_right_clip_len = 0
            if l_cigar[-1][0] == 4 or l_cigar[-1][0] == 5:  # right clipped
                b_right_clip = True
                i_right_clip_len = l_cigar[-1][1]

            if b_left_clip == True and b_right_clip == True:
                if (i_left_clip_len > MAX_CLIP_CLIP_LEN) and (i_right_clip_len > MAX_CLIP_CLIP_LEN):
                    return False, 0

        ####for the alignment (of the clipped read), if the mapped part is smaller than the clipped part,
        ####then skip
        n_total = 0
        n_map = 0
        for (type, lenth) in l_cigar:
            if type == 0:
                n_map += lenth
            if type != 2:  # deletion is not added to the total length
                n_total += lenth

        if n_map < (n_total * ratio_cutoff):  ########################require at least 3/4 of the seq is mapped !!!!!!!!
            return False, 0
        return True, n_map

        ####
        # Parse the re-aligned clipped parts, and find the
        # 1) max-left(right)-clip on reference, and max-clip pos on consensus
        # 2) get the potential TSD
        # 3) transduction
        ##parameter: 1) sf_clip_alignmt: the re-aligned clipped part, 2) minimum mapped ratio for realigned clip part

    ####
    # Here only consider the clipped reads within [-extend, extend] around the original insertion position
    ####
    def parse_clip_realignment(self, sf_clip_alignmt, bmapped_cutoff, mrmsk, i_flank_lenth):
        samfile = pysam.AlignmentFile(sf_clip_alignmt, "r", reference_filename=self.sf_reference)
        m_ins_pos_left = {}  ##this is to statistic the insertion positions for left-clip on the reference
        m_ins_pos_right = {}  # for the right clipped positions
        m_rep_pos_left = {}  # this is to save the posiion on the repeat consensus of the left clipped reads
        m_rep_pos_right = {}  # this is to save the posiion on the repeat consensus of the right clipped reads
        m_transduction = {}  # save the candidate transduction positions
        # nl_polyA = 0
        # nr_polyA = 0
        m_ins_pos_sample = {}  # sample information of the clip position
        m_polyA = {}
        for algnmt in samfile.fetch():
            if algnmt.is_unmapped == True:  ####skip the unmapped reads
                continue

            # In format:chrm~map_pos~FLAG_RIGHT(LEFT)_CLIP~is_reverse_complementary~insertion_position~rid~sample_id
            read_info = algnmt.query_name
            read_info_fields = read_info.split(SEPERATOR)

            ori_chrm = read_info_fields[0]
            ref_mpos = int(read_info_fields[1])
            ori_clip_flag = read_info_fields[2]
            ori_is_rc = read_info_fields[3]  ##############this is not used here now
            ori_insertion_pos = int(read_info_fields[4])
            sample_id = read_info_fields[6]  ####sample id
            #####Hard code here#########################################################################################
            if abs(ref_mpos - ori_insertion_pos) > 50:
                continue

            b_left = True  # left clip
            if ori_clip_flag == FLAG_RIGHT_CLIP:
                b_left = False  # right clip

            # first check whether the clipped part is qualified aligned
            l_cigar = algnmt.cigar
            b_clip_qualified_algned, n_map_bases = self.is_clipped_part_qualified_algnmt(l_cigar, bmapped_cutoff)

            if b_clip_qualified_algned == False:  # skip the unqualified re-aligned parts
                continue

            # map position on repeat consensus
            map_pos_copy = algnmt.reference_start  #map position on which repeat copy
            b_clip_part_rc = algnmt.is_reverse
            # id in format: chrm~start_pos~end_pos~sub-family
            mapped_copy_id = algnmt.reference_name  # map to which repeat copy

            mapped_copy_fields = mapped_copy_id.split(SEPERATOR)
            copy_chrm = mapped_copy_fields[0]
            copy_start = int(mapped_copy_fields[1])
            copy_end = int(mapped_copy_fields[2])
            # copy_sub_family=mapped_copy_fields[3]
            clipped_seq = algnmt.query_sequence
            ##if it is mapped to the left flank regions
            # here make sure it is not poly-A
            b_polya = self.contain_poly_A_T(clipped_seq, N_MIN_A_T)  # by default, at least 5A or 5T
            mapq = algnmt.mapping_quality

            if b_polya == True:
                if ori_chrm not in m_polyA:
                    m_polyA[ori_chrm] = {}
                if ori_insertion_pos not in m_polyA[ori_chrm]:
                    m_polyA[ori_chrm][ori_insertion_pos] = []
                    m_polyA[ori_chrm][ori_insertion_pos].append(0)
                    m_polyA[ori_chrm][ori_insertion_pos].append(0)
                if b_left == True:
                    m_polyA[ori_chrm][ori_insertion_pos][0] += 1
                else:
                    m_polyA[ori_chrm][ori_insertion_pos][1] += 1

            # for those mapped to the flanks
            # each insertion have two sides, so if both sides are transduction (very rare), then will be missed
            #
            if (map_pos_copy < (i_flank_lenth - 5)) or (
                        map_pos_copy > (copy_end - copy_start + i_flank_lenth - 5)):

                if (mapq >= MINIMAL_TRANSDUCT_MAPQ) and b_polya == False:  #possible transduction
                    if ori_chrm not in m_transduction:
                        m_transduction[ori_chrm] = {}
                    if ori_insertion_pos not in m_transduction[ori_chrm]:
                        m_transduction[ori_chrm][ori_insertion_pos] = {}
                    transduct_pos = "{0}~{1}~{2}".format(copy_chrm, copy_start, copy_end)
                    if transduct_pos not in m_transduction[ori_chrm][ori_insertion_pos]:
                        m_transduction[ori_chrm][ori_insertion_pos][transduct_pos] = []

                    # TDList:convert this to ref pos
                    m_transduction[ori_chrm][ori_insertion_pos][transduct_pos].append(map_pos_copy)
                continue
            #
            # calculate the position in consensus
            pos_in_consensus = -1
            if copy_chrm in mrmsk:
                if copy_start in mrmsk[copy_chrm]:
                    copy_record = mrmsk[copy_chrm][copy_start][0]
                    cns_start = copy_record[-2]
                    cns_end = copy_record[-1]
                    b_copy_rc = copy_record[1]

                    # first, get the clip position on the copy
                    if (b_clip_part_rc and b_left) or (b_clip_part_rc == False and b_left == False):
                        map_pos_copy += n_map_bases

                    # here need to pay attention whether the copy is masked "reverse complementary" or not.
                    # If the copy is masked "reverse complementary", then the position in consensus should be reversed
                    i_offset = map_pos_copy - i_flank_lenth
                    if b_copy_rc == True:
                        i_offset = (cns_end - cns_start) - map_pos_copy + i_flank_lenth

                    # convert the map_pos on repeat copy to position on consensus
                    pos_in_consensus = cns_start + i_offset

                    # if pos_in_consensus< -2:##########################################################################
                    #     print map_pos_copy, mapped_copy_id, copy_record, pos_in_consensus

            ###save the clip position on the refernece
            if b_left == True:  # left clipped, save to "m_ins_pos_left"
                if ori_chrm not in m_ins_pos_left:
                    m_ins_pos_left[ori_chrm] = {}
                if ori_insertion_pos not in m_ins_pos_left[ori_chrm]:
                    m_ins_pos_left[ori_chrm][ori_insertion_pos] = {}
                if ref_mpos not in m_ins_pos_left[ori_chrm][ori_insertion_pos]:
                    m_ins_pos_left[ori_chrm][ori_insertion_pos][ref_mpos] = []
                m_ins_pos_left[ori_chrm][ori_insertion_pos][ref_mpos].append(pos_in_consensus)  # note  "-1"
            else:
                if ori_chrm not in m_ins_pos_right:
                    m_ins_pos_right[ori_chrm] = {}
                if ori_insertion_pos not in m_ins_pos_right[ori_chrm]:
                    m_ins_pos_right[ori_chrm][ori_insertion_pos] = {}
                if ref_mpos not in m_ins_pos_right[ori_chrm][ori_insertion_pos]:
                    m_ins_pos_right[ori_chrm][ori_insertion_pos][ref_mpos] = []
                m_ins_pos_right[ori_chrm][ori_insertion_pos][ref_mpos].append(pos_in_consensus)  # note  "-1"

            ##save the clip position on the repeat consensus
            if b_left == True:
                if ori_chrm not in m_rep_pos_left:
                    m_rep_pos_left[ori_chrm] = {}
                if ori_insertion_pos not in m_rep_pos_left[ori_chrm]:
                    m_rep_pos_left[ori_chrm][ori_insertion_pos] = {}
                if pos_in_consensus not in m_rep_pos_left[ori_chrm][ori_insertion_pos]:
                    m_rep_pos_left[ori_chrm][ori_insertion_pos][pos_in_consensus] = 1
                else:
                    m_rep_pos_left[ori_chrm][ori_insertion_pos][pos_in_consensus] += 1
            else:
                if ori_chrm not in m_rep_pos_right:
                    m_rep_pos_right[ori_chrm] = {}
                if ori_insertion_pos not in m_rep_pos_right[ori_chrm]:
                    m_rep_pos_right[ori_chrm][ori_insertion_pos] = {}
                if pos_in_consensus not in m_rep_pos_right[ori_chrm][ori_insertion_pos]:
                    m_rep_pos_right[ori_chrm][ori_insertion_pos][pos_in_consensus] = 1
                else:
                    m_rep_pos_right[ori_chrm][ori_insertion_pos][pos_in_consensus] += 1

            # save the sample information
            if ori_chrm not in m_ins_pos_sample:
                m_ins_pos_sample[ori_chrm] = {}
            if ori_insertion_pos not in m_ins_pos_sample[ori_chrm]:
                m_ins_pos_sample[ori_chrm][ori_insertion_pos] = {}
            if sample_id not in m_ins_pos_sample[ori_chrm][ori_insertion_pos]:
                m_ins_pos_sample[ori_chrm][ori_insertion_pos][sample_id] = []
                m_ins_pos_sample[ori_chrm][ori_insertion_pos][sample_id].append(0)
                m_ins_pos_sample[ori_chrm][ori_insertion_pos][sample_id].append(0)
            if b_left == True:
                m_ins_pos_sample[ori_chrm][ori_insertion_pos][sample_id][0] += 1
            else:
                m_ins_pos_sample[ori_chrm][ori_insertion_pos][sample_id][1] += 1

        samfile.close()
        return m_ins_pos_left, m_ins_pos_right, m_ins_pos_sample, m_transduction, m_polyA

    #for this version, the clipped parte are realigned to the consensus, not all the repeat copies
    def parse_clip_realignment_consensus(self, sf_clip_alignmt, bmapped_cutoff):
        samfile = pysam.AlignmentFile(sf_clip_alignmt, "r", reference_filename=self.sf_reference)
        m_ins_pos_left = {}  ##this is to statistic the insertion positions for left-clip on the reference
        m_ins_pos_right = {}  # for the right clipped positions
        m_rep_pos_left = {}  # this is to save the posiion on the repeat consensus of the left clipped reads
        m_rep_pos_right = {}  # this is to save the posiion on the repeat consensus of the right clipped reads
        #m_transduction = {}  # save the candidate transduction positions

        m_ins_pos_sample = {}  # sample information of the clip position
        m_polyA = {}
        for algnmt in samfile.fetch():
            #also check the mapping quality
            # mapq = algnmt.mapping_quality
            # if mapq<MINIMUM_DISC_MAPQ:##############Here should be very careful for SVA and Alu!!!!!!!!!!!!!!!!!!!!!!!
            #     continue

            # In format:chrm~map_pos~FLAG_RIGHT(LEFT)_CLIP~is_reverse_complementary~insertion_position~rid~sample_id
            read_info = algnmt.query_name
            read_info_fields = read_info.split(SEPERATOR)

            ori_chrm = read_info_fields[0]
            ref_mpos = int(read_info_fields[1])
            ori_clip_flag = read_info_fields[2]
            ori_is_rc = read_info_fields[3]  ##############this is not used here now
            ori_insertion_pos = int(read_info_fields[4])
            sample_id = read_info_fields[6]  ####sample id

            #####Hard code here########################################################################################
            if abs(ref_mpos - ori_insertion_pos) > NEARBY_CLIP: #only focus on the clipped reads nearby
                continue

            b_left = True  # left clip
            if ori_clip_flag == FLAG_RIGHT_CLIP:
                b_left = False  # right clip

            clipped_seq = algnmt.query_sequence
            ##if it is mapped to the left flank regions
            # here make sure it is not poly-A
            #b_polya = self.contain_poly_A_T(clipped_seq, N_MIN_A_T)  # by default, at least 5A or 5T
            b_polya=self.is_consecutive_polyA_T(clipped_seq)
            if b_polya == True:
                if ori_chrm not in m_polyA:
                    m_polyA[ori_chrm] = {}
                if ori_insertion_pos not in m_polyA[ori_chrm]:
                    m_polyA[ori_chrm][ori_insertion_pos] = []
                    m_polyA[ori_chrm][ori_insertion_pos].append(0)
                    m_polyA[ori_chrm][ori_insertion_pos].append(0)
                if b_left == True:
                    m_polyA[ori_chrm][ori_insertion_pos][0] += 1
                else:
                    m_polyA[ori_chrm][ori_insertion_pos][1] += 1

            if algnmt.is_unmapped == True:  ####skip the unmapped reads
                continue
            # first check whether the clipped part is qualified aligned
            l_cigar = algnmt.cigar
            b_clip_qualified_algned, n_map_bases = self.is_clipped_part_qualified_algnmt(l_cigar, bmapped_cutoff)

            if b_clip_qualified_algned == False:  # skip the unqualified re-aligned parts
                continue

            # calculate the position in consensus
            pos_in_consensus = algnmt.reference_start  #map position on which repeat copy
            b_clip_part_rc = algnmt.is_reverse
            # first, get the clip position on the copy
            if (b_clip_part_rc and b_left) or (b_clip_part_rc == False and b_left == False):
                pos_in_consensus += n_map_bases

            ###save the clip position on the refernece
            if b_left == True:  # left clipped, save to "m_ins_pos_left"
                if ori_chrm not in m_ins_pos_left:
                    m_ins_pos_left[ori_chrm] = {}
                if ori_insertion_pos not in m_ins_pos_left[ori_chrm]:
                    m_ins_pos_left[ori_chrm][ori_insertion_pos] = {}
                if ref_mpos not in m_ins_pos_left[ori_chrm][ori_insertion_pos]:
                    m_ins_pos_left[ori_chrm][ori_insertion_pos][ref_mpos] = []
                m_ins_pos_left[ori_chrm][ori_insertion_pos][ref_mpos].append(pos_in_consensus)  # note  "-1"
            else:
                if ori_chrm not in m_ins_pos_right:
                    m_ins_pos_right[ori_chrm] = {}
                if ori_insertion_pos not in m_ins_pos_right[ori_chrm]:
                    m_ins_pos_right[ori_chrm][ori_insertion_pos] = {}
                if ref_mpos not in m_ins_pos_right[ori_chrm][ori_insertion_pos]:
                    m_ins_pos_right[ori_chrm][ori_insertion_pos][ref_mpos] = []
                m_ins_pos_right[ori_chrm][ori_insertion_pos][ref_mpos].append(pos_in_consensus)  # note  "-1"

            ##save the clip position on the repeat consensus
            if b_left == True:
                if ori_chrm not in m_rep_pos_left:
                    m_rep_pos_left[ori_chrm] = {}
                if ori_insertion_pos not in m_rep_pos_left[ori_chrm]:
                    m_rep_pos_left[ori_chrm][ori_insertion_pos] = {}
                if pos_in_consensus not in m_rep_pos_left[ori_chrm][ori_insertion_pos]:
                    m_rep_pos_left[ori_chrm][ori_insertion_pos][pos_in_consensus] = 1
                else:
                    m_rep_pos_left[ori_chrm][ori_insertion_pos][pos_in_consensus] += 1
            else:
                if ori_chrm not in m_rep_pos_right:
                    m_rep_pos_right[ori_chrm] = {}
                if ori_insertion_pos not in m_rep_pos_right[ori_chrm]:
                    m_rep_pos_right[ori_chrm][ori_insertion_pos] = {}
                if pos_in_consensus not in m_rep_pos_right[ori_chrm][ori_insertion_pos]:
                    m_rep_pos_right[ori_chrm][ori_insertion_pos][pos_in_consensus] = 1
                else:
                    m_rep_pos_right[ori_chrm][ori_insertion_pos][pos_in_consensus] += 1

            # save the sample information
            if ori_chrm not in m_ins_pos_sample:
                m_ins_pos_sample[ori_chrm] = {}
            if ori_insertion_pos not in m_ins_pos_sample[ori_chrm]:
                m_ins_pos_sample[ori_chrm][ori_insertion_pos] = {}
            if sample_id not in m_ins_pos_sample[ori_chrm][ori_insertion_pos]:
                m_ins_pos_sample[ori_chrm][ori_insertion_pos][sample_id] = []
                m_ins_pos_sample[ori_chrm][ori_insertion_pos][sample_id].append(0)
                m_ins_pos_sample[ori_chrm][ori_insertion_pos][sample_id].append(0)
            if b_left == True:
                m_ins_pos_sample[ori_chrm][ori_insertion_pos][sample_id][0] += 1
            else:
                m_ins_pos_sample[ori_chrm][ori_insertion_pos][sample_id][1] += 1

        samfile.close()
        return m_ins_pos_left, m_ins_pos_right, m_rep_pos_left, m_rep_pos_right, m_ins_pos_sample, m_polyA


    ####
    # parse the re-aligned discordant reads
    # return two dict: 1) the aligned position on the repeat; 2) candidate transductions
    def parse_disc_algnmt(self, sf_disc_alignmt, bmapped_cutoff, mrmsk, i_flank_lenth):
        samfile = pysam.AlignmentFile(sf_disc_alignmt, "r", reference_filename=self.sf_reference)
        m_disc_pos = {}
        m_transduction_disc = {}
        m_polyA = {}
        m_disc_sample = {}
        for algnmt in samfile.fetch():
            if algnmt.is_unmapped == True:  ####skip the unmapped reads
                continue
            mapq = algnmt.mapping_quality
            ####also, for clipped mapped reads, need to check the clipped parts whether can be split to two parts!!!!!!
            #
            read_info = algnmt.query_name  # in format: #read_id~is_first~s_insertion_chrm~s_insertion_pos~sample_id
            read_info_fields = read_info.split(SEPERATOR)
            ins_chrm = read_info_fields[-3]
            ins_pos = int(read_info_fields[-2])
            sample_id = read_info_fields[-1]

            # first check whether read is qualified mapped
            l_cigar = algnmt.cigar
            b_clip_qualified_algned, n_map_bases = self.is_clipped_part_qualified_algnmt(l_cigar, bmapped_cutoff)
            if b_clip_qualified_algned == False:  # skip the unqualified re-aligned parts
                continue

            map_pos_copy = algnmt.reference_start  # map position on repeat copy
            mapped_copy_id = algnmt.reference_name  # map to which repeat copy
            mapped_copy_fields = mapped_copy_id.split(SEPERATOR)
            copy_chrm = mapped_copy_fields[0]
            copy_start = int(mapped_copy_fields[1])
            copy_end = int(mapped_copy_fields[2])

            read_seq = algnmt.query_sequence
            b_polya = self.contain_poly_A_T(read_seq, N_MIN_A_T)
            if b_polya == True:
                if ins_chrm not in m_polyA:
                    m_polyA[ins_chrm] = {}
                if ins_pos not in m_polyA[ins_chrm]:
                    m_polyA[ins_chrm][ins_pos] = 1
                else:
                    m_polyA[ins_chrm][ins_pos] += 1
            # mapped to flank region of the copy
            # transduction situation
            if ((mapq >= MINIMAL_TRANSDUCT_MAPQ) and (map_pos_copy < (i_flank_lenth - 15)) or (
                        map_pos_copy > (copy_end - copy_start + i_flank_lenth - 15))) and b_polya == False:
                if ins_chrm not in m_transduction_disc:
                    m_transduction_disc[ins_chrm] = {}
                if ins_pos not in m_transduction_disc[ins_chrm]:
                    m_transduction_disc[ins_chrm][ins_pos] = {}
                transduct_pos = "{0}~{1}~{2}".format(copy_chrm, copy_start, copy_end)
                if transduct_pos not in m_transduction_disc[ins_chrm][ins_pos]:
                    m_transduction_disc[ins_chrm][ins_pos][transduct_pos] = []
                m_transduction_disc[ins_chrm][ins_pos][transduct_pos].append(map_pos_copy)
                continue

            map_pos_copy+=(len(read_seq)/2) ####use the middle as the map position
            ###not transduction, then find the position on consensus
            if copy_chrm in mrmsk:
                if copy_start in mrmsk[copy_chrm]:
                    copy_record = mrmsk[copy_chrm][copy_start][0]
                    cns_start = copy_record[-2]
                    cns_end = copy_record[-1]
                    b_copy_rc = copy_record[1]

                    i_offset = map_pos_copy - i_flank_lenth
                    if b_copy_rc == True:
                        i_offset = (cns_end - cns_start) - map_pos_copy + i_flank_lenth
                    # convert the map_pos on repeat copy to position on consensus
                    pos_in_consensus = cns_start + i_offset
                    if ins_chrm not in m_disc_pos:
                        m_disc_pos[ins_chrm] = {}
                    if ins_pos not in m_disc_pos[ins_chrm]:
                        m_disc_pos[ins_chrm][ins_pos] = []
                    m_disc_pos[ins_chrm][ins_pos].append(pos_in_consensus)

            ###
            if ins_chrm not in m_disc_sample:
                m_disc_sample[ins_chrm] = {}
            if ins_pos not in m_disc_sample[ins_chrm]:
                m_disc_sample[ins_chrm][ins_pos] = {}
            if sample_id not in m_disc_sample[ins_chrm][ins_pos]:
                m_disc_sample[ins_chrm][ins_pos][sample_id] = 1
            else:
                m_disc_sample[ins_chrm][ins_pos][sample_id] += 1

        samfile.close()
        return m_disc_pos, m_disc_sample, m_transduction_disc, m_polyA


    #For this version, reads are aligned to the consensus sequence
    #Note that the poly-A part in the consensus has been removed, so it's possible reads will be clipped mapped (at end)
    def parse_disc_algnmt_consensus(self, sf_disc_alignmt, bmapped_cutoff):
        samfile = pysam.AlignmentFile(sf_disc_alignmt, "r", reference_filename=self.sf_reference)
        m_disc_pos = {}
        m_polyA = {}
        m_disc_sample = {}
        for algnmt in samfile.fetch():
            # mapq = algnmt.mapping_quality
            # if mapq<MINIMUM_DISC_MAPQ:##############Here should be very careful for SVA and Alu!!!!!!!!!!!!!!!!!!!!!!!#
            #     continue
            ####also, for clipped mapped reads, need to check the clipped parts whether can be split to two parts!!!!!!
            #
            read_info = algnmt.query_name  # in format: #read_id~is_first~s_insertion_chrm~s_insertion_pos~sample_id
            read_info_fields = read_info.split(SEPERATOR)
            ins_chrm = read_info_fields[-3]
            ins_pos = int(read_info_fields[-2])
            sample_id = read_info_fields[-1]

            read_seq = algnmt.query_sequence
            b_polya = self.is_consecutive_polyA_T(read_seq)
            if b_polya == True:
                if ins_chrm not in m_polyA:
                    m_polyA[ins_chrm] = {}
                if ins_pos not in m_polyA[ins_chrm]:
                    m_polyA[ins_chrm][ins_pos] = 1
                else:
                    m_polyA[ins_chrm][ins_pos] += 1

            if algnmt.is_unmapped == True:  ####skip the unmapped reads
                continue
            # first check whether read is qualified mapped
            l_cigar = algnmt.cigar
            b_clip_qualified_algned, n_map_bases = self.is_clipped_part_qualified_algnmt(l_cigar, bmapped_cutoff)
            if b_clip_qualified_algned == False:  # skip the unqualified re-aligned parts
                continue
            # map position on repeat consensus
            pos_in_consensus = algnmt.reference_start + len(read_seq)/2
            if ins_chrm not in m_disc_pos:
                m_disc_pos[ins_chrm] = {}
            if ins_pos not in m_disc_pos[ins_chrm]:
                m_disc_pos[ins_chrm][ins_pos] = []
            m_disc_pos[ins_chrm][ins_pos].append(pos_in_consensus)

            ###
            if ins_chrm not in m_disc_sample:
                m_disc_sample[ins_chrm] = {}
            if ins_pos not in m_disc_sample[ins_chrm]:
                m_disc_sample[ins_chrm][ins_pos] = {}
            if sample_id not in m_disc_sample[ins_chrm][ins_pos]:
                m_disc_sample[ins_chrm][ins_pos][sample_id] = 1
            else:
                m_disc_sample[ins_chrm][ins_pos][sample_id] += 1

        samfile.close()
        return m_disc_pos, m_disc_sample, m_polyA


    # # Given annotation file, flank region length,
    # def re_locate_in_consensus(self, xannotation, i_flank_length, chrm, map_pos):
    #     return
    #####

    def is_poly_A_T(self, seq):  ###for a temp version here
        cnt_A = 0
        cnt_T = 0
        for ch in seq:
            if ch == 'A' or ch == 'a':
                cnt_A += 1
            elif ch == 'T' or ch == 't':
                cnt_T += 1
        cnt_cutoff = (len(seq) * 0.75)
        if (cnt_A > cnt_cutoff) or (cnt_T > cnt_cutoff):
            return True
        return False

    #
    def contain_poly_A_T(self, seq, n_min_cnt):
        max_A = 0
        max_T = 0
        cum_A = 0
        cum_T = 0
        for ch in seq:
            if ch == 'A' or ch == 'a':
                cum_A += 1
            else:
                if cum_A > max_A:
                    max_A = cum_A
                cum_A = 0

            if ch == 'T' or ch == 't':
                cum_T += 1
            else:
                if cum_T > max_T:
                    max_T = cum_T
                cum_T = 0
        if cum_A > max_A:
            max_A = cum_A
        if cum_T > max_T:
            max_T = cum_T

        if max_A >= n_min_cnt or max_T >= n_min_cnt:
            return True
        else:
            return False


    def is_consecutive_polyA_T(self, seq):
        if ("AAAAA" in seq) or ("TTTTT" in seq) or ("AATAA" in seq) or ("TTATT" in seq):
            return True
        else:
            return False
####
    # for given a lits of candidate sites, and a given list of bam files
    # get all the clipped and discordant reads
    def collect_clipped_disc_reads(self, sf_candidate_list, extnd, bin_size, sf_clip_fq, sf_disc_fa):
        with open(sf_clip_fq, "w") as fout_clip_fq, open(sf_disc_fa, "w") as fout_disc_fa:
            # first collect all the clipped and discordant reads
            ####test only####
            sample_cnt = 0
            with open(self.sf_bam_list) as fin_bam_list:
                for line in fin_bam_list:  # for each bam file
                    sf_bam = line.rstrip()
                    xclip_disc = XClipDisc(sf_bam, self.working_folder, self.n_jobs, self.sf_reference)
                    sf_disc_fa_tmp = self.working_folder + "temp_disc.fa" + str(sample_cnt)
                    sf_clip_fa_tmp = self.working_folder + "temp_clip.fq" + str(sample_cnt)

                    ####test only
                    xclip_disc.collect_clipped_disc_reads_of_given_list(sf_candidate_list, extnd, bin_size,
                                                                        sf_clip_fa_tmp, sf_disc_fa_tmp)
                    with open(sf_disc_fa_tmp) as fin_tmp_fa:
                        n_cnt = 0
                        for line in fin_tmp_fa:
                            if n_cnt % 2 == 0:
                                line = line.rstrip() + SEPERATOR + str(sample_cnt) + "\n"  ###here add the sample id
                            fout_disc_fa.write(line)
                            n_cnt += 1
                    with open(sf_clip_fa_tmp) as fin_clip_fq:
                        n_cnt = 0
                        for line in fin_clip_fq:
                            if n_cnt % 4 == 0:
                                line = line.rstrip() + SEPERATOR + str(sample_cnt) + "\n"
                            fout_clip_fq.write(line)
                            n_cnt += 1
                    sample_cnt += 1

    # re-align the collected clipped and discordant reads
    def realign_clipped_reads(self, sf_ref, sf_reads, sf_out_sam):
        cmd = "{0} mem -t {1} -T {2} -k {3} {4} {5} > {6}".format(BWA_PATH, self.n_jobs, BWA_REALIGN_CUTOFF,
                                                                  BWA_REALIGN_CUTOFF, sf_ref, sf_reads, sf_out_sam)
        Popen(cmd, shell=True, stdout=PIPE).communicate()

    # re-align the collected clipped and discordant reads
    def realign_disc_reads(self, sf_ref, sf_reads, sf_out_sam):
        cmd = "{0} mem -t {1} {2} {3} > {4}".format(BWA_PATH, self.n_jobs, sf_ref, sf_reads, sf_out_sam)
        Popen(cmd, shell=True, stdout=PIPE).communicate()

    # check the consistency of clipped and discordant positions on realigned repeat copies
    # this is based on the observation that: reads aligned to repeat copies also follow normal distribution
    def check_clip_disc_pos_consistency_on_rep(self):
        return

####
####