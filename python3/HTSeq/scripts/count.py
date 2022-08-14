# NOTES:
# latest commit changed sequences in chromosomes dictionary to bytearray from bytes (bytestrings) for making N-masking genome more efficient
# should not have broken things hopefully

VCF_FILENAME = "/home/linux/RNAseq/mgp_REL2005_snps_indels.vcf"
BED_FILENAME = "/home/linux/RNAseq/SNPs.bed"

# Some example commandline invokations using python3 to illustrate use and expected results below.
# testrun:
#           allelic_count.py --SNP-allelic-analysis --SNP-vcf-filename=/home/linux/RNAseq/mgp_REL2005_snps_indels.vcf --SNP-bed-filename=/home/linux/RNAseq/SNPs_NEW.bed --SNP-use-bed-file-cache --SNP-ref-sample=129S1_SvImJ --SNP-ref-sample=129S5SvEvBrd --SNP-ref-sample=C57BL_6NJ --SNP-ref-sample=CBA_J --SNP-ref-sample=DBA_2J --SNP-alt-sample=CAST_EiJ --SNP-force-biallelic --additional-attr=gene_name --format=sam --mode=union --stranded=reverse --minaqual=2 --counts_output=HATXCDK8KOnodox_alleleic_counts_NEW /home/linux/RNAseq/Fastq/HATXCDK8KOnodox.sa /home/linux/RNAseq/Mus_musculus.GRCm38.68.gtf
#
# python3 allelic_count.py --SNP-allelic-analysis --SNP-fast-vcf-reader --SNP-vcf-filename=/home/linux/RNAseq/mgp_REL2005_snps_indels.vcf --SNP-bed-filename=/home/linux/RNAseq/SNPs_NEW.bed --SNP-use-bed-file-cache --SNP-ref-sample=129S1_SvImJ --SNP-ref-sample=129S5SvEvBrd --SNP-ref-sample=C57BL_6NJ --SNP-ref-sample=CBA_J --SNP-ref-sample=DBA_2J --SNP-alt-sample=CAST_EiJ --SNP-force-biallelic --additional-attr=gene_name --format=sam --mode=union --stranded=reverse --minaqual=2 --counts_output=HATXCDK8KOnodox_alleleic_counts_NEW /home/linux/RNAseq/Fastq/HATXCDK8KOnodox.sa /home/linux/RNAseq/Mus_musculus.GRCm38.68.gtf
#
# CONSISTENCY TEST
#
# python3 allelic_count.py --check-consistency --num-consistency-tests=100000 --genome-fasta-filename=/home/linux/RNAseq/Mus_musculus.GRCm38.68.dna.toplevel.fa/Mus_musculus.GRCm38.68.dna.chromosome.1.fa --SNP-vcf-filename=/home/linux/RNAseq/mgp_REL2005_snps_indels.vcf --SNP-vcf-file-offset=1 /home/linux/RNAseq/Fastq/HATXCDK8KOnodox.sa /home/linux/RNAseq/Mus_musculus.GRCm38.68.gtf
# for checking VCF file only use:
# python3 allelic_count.py --check-consistency --num-consistency-tests=100000 --genome-fasta-filename=/home/linux/RNAseq/Mus_musculus.GRCm38.68.dna.toplevel.fa/Mus_musculus.GRCm38.68.dna.chromosome.1.fa --SNP-vcf-filename=/home/linux/RNAseq/mgp_REL2005_snps_indels.vcf --SNP-vcf-file-offset=1 None None
#
# 15951727 SNPs found in 92510146 vcf records
# SNP information written to BED cache file /home/linux/RNAseq/SNPs_NEW.bed
# OLD:
# It seem to work but there is a problem with the CAST reads!
#
# 33874599 alignments  processed.
# Polymorphism report: 3397358 assigned: 3392173 reference / 5185 CAST, 5398878 alleles unassigned of 8796236 tested.
#
#                     4375382 ref, 6113 CAST, 2706 unassigned SNPs total
#
# With -1/+1 extended cigarop.ref_iv a few more SNPs are found - is not from double counting as reads are assigned by voting and double votes do not increase assignments
#
# 33874599 alignments  processed.
# Polymorphism report: 3397393 assigned: 3391962 reference / 5431 CAST, 5398843 alleles unassigned of 8796236 tested.
#
#                     4375874 ref, 11047 CAST, 4550 unassigned SNPs total
#
# With also -1/+1 extended cigarop.iv for quick testing if SNPs overlap read alignment interval even a very few more turn up
# - note that the reference does not get any more as SNPs match on ref genome, but CAST does get some as these are mismatches that are not included in M operations
# - also a few more unassigned SNPs
#
# 33874599 alignments  processed.
# Polymorphism report: 3398038 assigned: 3391962 reference / 6076 CAST, 5433209 alleles unassigned of 8831247 tested.
#
#                     4375874 ref, 11692 CAST, 5614 unassigned SNPs total
#
# TXGFPCASTHIRAKOnodox.sa dataset contains CAST genome, analysis shows a small (10.6%, SNPs, 9.3% reads) underrepresentation of CAST (which maybe the missing Xcast chromosome?)
#
# testrun:
#           count_with_alleles.py --allelic-analysis --additional-attr=gene_name --format=sam --mode=union --stranded=reverse --minaqual=2 --counts_output TXGFPCASTHIRAKOnodox_alleleic_counts /home/linux/RNAseq/Fastq/TXGFPCASTHIRAKOnodox.sa /home/linux/RNAseq/Mus_musculus.GRCm38.68.gtf
#
# 31955143 alignments  processed.
# Polymorphism report: 3269872 assigned: 1786378 reference / 1483494 CAST, 5394715 alleles unassigned of 8664587 tested.
#
#                     2298694 ref, 1857067 CAST, 5587 unassigned SNPs total
#
#           count_with_alleles.py --allelic-analysis --additional-attr=gene_name --add-chromosome-information --format=sam --mode=union --stranded=reverse --minaqual=2 --counts_output TXGFPCASTHIRAKOnodox_alleleic_counts /home/linux/RNAseq/Fastq/TXGFPCASTHIRAKOnodox.sa /home/linux/RNAseq/Mus_musculus.GRCm38.68.gtf
#
# ADD FEATURE --SNP-min-base-qual ------------------------------------------------------------------------------------------------------------
#
# SNP min base qual 30
#
# python3 allelic_count.py --SNP-allelic-analysis --SNP-fast-vcf-reader --SNP-min-base-qual=30 --SNP-vcf-filename=/home/linux/RNAseq/mgp_REL2005_snps_indels.vcf --SNP-bed-filename=/home/linux/RNAseq/SNPs_NEW.bed --SNP-use-bed-file-cache --SNP-ref-sample=129S1_SvImJ --SNP-ref-sample=129S5SvEvBrd --SNP-ref-sample=C57BL_6NJ --SNP-ref-sample=CBA_J --SNP-ref-sample=DBA_2J --SNP-alt-sample=CAST_EiJ --SNP-force-biallelic --additional-attr=gene_name --format=sam --mode=union --stranded=reverse --minaqual=2 --counts_output=HATXCDK8KOnodox_alleleic_counts_NEW /home/linux/RNAseq/Fastq/HATXCDK8KOnodox.sa /home/linux/RNAseq/Mus_musculus.GRCm38.68.gtf
#3 3874599 alignments  processed.
#Polymorphism report: 3257794 assigned: 3252199 reference / 5595 alternate, 5573453 alleles unassigned of 8831247 tested.
#
#                     4156672 ref, 10540 alt, 225968 unassigned SNPs total
#
# SNP min base qual 25
#
# python3 allelic_count.py --SNP-allelic-analysis --SNP-fast-vcf-reader --SNP-min-base-qual=25 --SNP-vcf-filename=/home/linux/RNAseq/mgp_REL2005_snps_indels.vcf --SNP-bed-filename=/home/linux/RNAseq/SNPs_NEW.bed --SNP-use-bed-file-cache --SNP-ref-sample=129S1_SvImJ --SNP-ref-sample=129S5SvEvBrd --SNP-ref-sample=C57BL_6NJ --SNP-ref-sample=CBA_J --SNP-ref-sample=DBA_2J --SNP-alt-sample=CAST_EiJ --SNP-force-biallelic --additional-attr=gene_name --format=sam --mode=union --stranded=reverse --minaqual=2 --counts_output=HATXCDK8KOnodox_alleleic_counts_NEW /home/linux/RNAseq/Fastq/HATXCDK8KOnodox.sa /home/linux/RNAseq/Mus_musculus.GRCm38.68.gtf
#
# 33874599 alignments  processed.
# Polymorphism report: 3356331 assigned: 3350600 reference / 5731 alternate, 5474916 alleles unassigned of 8831247 tested.
#
#                     4310315 ref, 10982 alt, 71883 unassigned SNPs total


# HTSeq vcf reader 18:57 - 21:13 15951727 SNPs found in 92510146 vcf records - 2 hours and 15 minutes
# fast vcf reader  21:22 - 21:47 15955367 Polymorphisms selected in 92510146 vcf records. - 25 minutes
# bed read 21:53 - 21:54 - 1 minute only
# GFT read 21:54 - 21:55 - 1 minute only
# BAM proc 21:55 - 22:32 - 37 minutes
# 33874599 alignments  processed.
# Polymorphism report: 3398038 assigned: 3391962 reference / 6076 alternate, 5433209 alleles unassigned of 8831247 tested.

# Reference genome fasts N-masking:
# python3 allelic_count.py --genome-nmask --genome-fasta-filename=/home/linux/RNAseq/Mus_musculus.GRCm38.68.dna.toplevel.fa/Mus_musculus.GRCm38.68.dna.chromosome.1.fa --genome-nmask-filename=Mus_musculus.GRCm38.68.dna.chromosome.N_masked.fa --SNP-vcf-filename=/home/linux/RNAseq/mgp_REL2005_snps_indels.vcf --SNP-vcf-file-offset=1 None None

import sys
import argparse
import operator
import itertools
import warnings
import traceback
import os
import os.path
import multiprocessing
import pysam
import random

import HTSeq


base2shft = {"A": 6, "C": 4, "G": 2, "T": 0}

def check_ref_alt(base, enc_snp):
    if base in base2shft:
        a = (enc_snp >> base2shft[base]) & 3
        return "ref" if a==1 else ("alt" if a==2 else "")
    return ""

def print_vcf_info(vcf_file_name):
    vcf_file = HTSeq.VCF_Reader(vcf_file_name)
    vcf_file.parse_meta()
    vcf_file.make_info_dict()
    print("Reading info from", vcf_file_name)
    print("Reference :", vcf_file.metadata["reference"])
    print("Command   :", vcf_file.metadata["bcftoolsCommand"])
    print("Fileformat:", vcf_file.metadata["fileformat"])
    for morphism in vcf_file:
        print(len(morphism.samples), "Samples :", " ".join(list(morphism.samples.keys())))
        break

def read_vcf(vcf_file_name, ref_samples, alt_samples, omit_ref_genome=False, force_biallelic=False, vcf_file_offset=1, silent=False):  
    SNPs = HTSeq.GenomicArray("auto", stranded=False, typecode='i')
    with open(vcf_file_name) as vcf_file:
        if not silent:
            print("Reading polymorphism information from VCF file:", vcf_file_name)
            print("Using genomic position offset", vcf_file_offset)
        gt_idx=0
        fi_idx=-1
        n=0
        num_polymorphisms=0
        for line in vcf_file:
            line=line.strip()
            if len(line) < 2:
                continue
            elif line[0]=="#":
                if line[1]=="#":            # meta lines
                    if not silent:
                        print(line)
                else:                       # header line
                    header={e:n for n,e in enumerate(line[1:].split("\t"))} # dictionary field name to index
                    if not silent:
                        print("HEADER -> data field index\n","\n".join([k+":"+str(header[k]) for k in header.keys()]))
            else:                           # data lines
                data_fields=line.split("\t")
                if data_fields[header["FORMAT"]][:2]!="GT" or data_fields[header["FORMAT"]][-2:]!="FI":
                    frmt_fields = data_fields[header["FORMAT"]].split(":")
                    gt_idx = frmt_fields.index("GT")
                    fi_idx = frmt_fields.index("FI")
                allele_sequences = [ data_fields[header["REF"]] ] + data_fields[header["ALT"]].split(",")
                ref_alleles=set([allele_sequences[int(gl[0])] for gl in [ sample_data[gt_idx].replace("|","/").split("/") for sample_data in [ data_fields[header[sample]].split(":") for sample in ref_samples ] if sample_data[fi_idx]=="1" ] if len(set(gl))==1 and gl[0]!="." and len(allele_sequences[int(gl[0])]) == 1])
                if not omit_ref_genome and data_fields[header["FILTER"]] == "PASS" and len(data_fields[header["REF"]]) == 1:
                    ref_alleles |= { data_fields[header["REF"]] }
                alt_alleles=set([allele_sequences[int(gl[0])] for gl in [ sample_data[gt_idx].replace("|","/").split("/") for sample_data in [ data_fields[header[sample]].split(":") for sample in alt_samples ] if sample_data[fi_idx]=="1" ] if len(set(gl))==1 and gl[0]!="." and len(allele_sequences[int(gl[0])]) == 1])
                if len(ref_alleles)>0 and len(alt_alleles)>0 and len(ref_alleles.intersection(alt_alleles)) == 0:
                    if not force_biallelic or (len(alt_alleles) == 1 and len(ref_alleles) == 1):
                        num_polymorphisms+=1
                        SNPs[ HTSeq.GenomicPosition(data_fields[header["CHROM"]], int(data_fields[header["POS"]]) -vcf_file_offset, ".") ] = sum( [1<<base2shft[b] for b in ref_alleles] + [2<<base2shft[b] for b in alt_alleles] )
                n+=1
                if not silent:
                    if n % 1000000 == 0:
                        print(n, "vcf records processed,", num_polymorphisms, "SNPs selected.")
    if not silent:
        print(num_polymorphisms,"Polymorphisms selected in", n, "vcf records.") #18965072 Polymorphisms found in 92510146 records / with LABSTRAINS
    return SNPs

def encode_vcf_file(vcf_file_name, ref_list, alt_list, omit_ref_genome=False, force_biallelic=False, vcf_file_offset=1, silent=False):

    def check_genotype(gt):
        a = gt.replace("|","/").split("/")
        if len( set(a) ) == 1 and a[0] != ".":
            return a[0]
        return ""

    vcf_file = HTSeq.VCF_Reader(vcf_file_name)
    SNPs = HTSeq.GenomicArray("auto", stranded=False, typecode='i')
    vcf_file.parse_meta()
    vcf_file.make_info_dict()
    if not silent:
        print("Reading polymorphisms from", vcf_file_name)
        print("Reference :", vcf_file.metadata["reference"])
        print("Command   :", vcf_file.metadata["bcftoolsCommand"])
        print("Fileformat:", vcf_file.metadata["fileformat"])
    num_polymorphisms = 0
    n = 0
    for morphism in vcf_file:
        if not silent:
            if n == 0:
                print(len(morphism.samples), "Samples :", " ".join(list(morphism.samples.keys())))
                print("Ref samples:", " ".join(ref_list))
                print(" -> against:", " ".join(alt_list))
            else:
                if n % 1000000 == 0:
                    print(n, "vcf records processed -", num_polymorphisms, "SNPs selected (", int(100*num_polymorphisms/n), "% )")
        n += 1
        allele_seq = [ e.upper() for e in [ morphism.ref ] + morphism.alt ]  # .ref is string .alt is list
        ref_alleles = set([ allele_seq[int(gt)] for s, gt in [ (sample, check_genotype(morphism.samples[sample]["GT"])) for sample in ref_list if morphism.samples[sample]["FI"] == "1" ] if gt!="" and len(allele_seq[int(gt)])==1 ])
        alt_alleles = set([ allele_seq[int(gt)] for s, gt in [ (sample, check_genotype(morphism.samples[sample]["GT"])) for sample in alt_list if morphism.samples[sample]["FI"] == "1" ] if gt!="" and len(allele_seq[int(gt)])==1 ])
        if not omit_ref_genome and morphism.filter == 'PASS':                # handle reference genome
            if len(morphism.ref) == 1:
                ref_alleles |= {morphism.ref}
        if len(ref_alleles) < 1 or len(alt_alleles) < 1:                     # no valid allele for SNP for either ref or alt list
            continue
        elif force_biallelic:
            if len(ref_alleles) != 1 or len(alt_alleles) != 1:               # only accept biallelic SNPs regarding ref and alt lists
                continue
        if len(ref_alleles.intersection(alt_alleles)) > 0:                   # alleles overlap SNPs are ambigous between ref and alt lists
            continue
        morphism.pos.strand = "."                                            # force non stranded
        morphism.pos.start-=vcf_file_offset                                  # subtract offset to adjust from 1-based to 0-based index <-- potential bug in HTSeq <<<<<<<<<<<<<<< BUG REPORT ???
        SNPs[morphism.pos] = sum( [1<<base2shft[br] for br in ref_alleles] + [2<<base2shft[ba] for ba in alt_alleles] )
        num_polymorphisms += 1
    if not silent:
        print(num_polymorphisms, "SNPs found in", n, "vcf records")
    return SNPs

class UnknownChrom(Exception):
    pass

def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2


def count_reads_single_file(
        isam,
        sam_filename,
        features,
        feature_attr,
        order,
        max_buffer_size,
        stranded,
        overlap_mode,
        multimapped_mode,
        secondary_alignment_mode,
        supplementary_alignment_mode,
        feature_type,
        id_attribute,
        additional_attributes,
        quiet,
        minaqual,
        samout_format,
        samout_filename,
        allele_counts,   # <---------------------------------- allelic analysis
        SNPs,
        snp_min_qual,
        ):

    def write_to_samout(r, assignment, samoutfile, template=None):
        if samoutfile is None:
            return
        if not pe_mode:
            r = (r,)
        for read in r:
            if read is not None:
                read.optional_fields.append(('XF', assignment))
                if samout_format in ('SAM', 'sam'):
                    samoutfile.write(read.get_sam_line() + "\n")
                else:
                    samoutfile.write(read.to_pysam_AlignedSegment(template))

    allelicreads=0
    unreads=0
    totalreads=0
    refreads=0
    altreads=0
    alt_total=0
    ref_total=0
    un_total=0
    try:
        if sam_filename == "-":
            read_seq_file = HTSeq.BAM_Reader(sys.stdin)
        else:
            read_seq_file = HTSeq.BAM_Reader(sam_filename)

        # Get template for output BAM
        if samout_filename is None:
            template = None
            samoutfile = None
        elif samout_format in ('bam', 'BAM'):
            template = read_seq_file.get_template()
            samoutfile = pysam.AlignmentFile(
                    samout_filename, 'wb',
                    template=template,
                    )
        else:
            template = None
            samoutfile = open(samout_filename, 'w')

        read_seq_iter = iter(read_seq_file)
        # Catch empty BAM files
        try:
            first_read = next(read_seq_iter)
            pe_mode = first_read.paired_end
        # FIXME: catchall can hide subtle bugs
        except:
            first_read = None
            pe_mode = False
        if first_read is not None:
            read_seq = itertools.chain([first_read], read_seq_iter)
        else:
            read_seq = []
    except:
        sys.stderr.write(
            "Error occured when reading beginning of SAM/BAM file.\n")
        raise

    # CIGAR match characters (including alignment match, sequence match, and
    # sequence mismatch
    com = ('M', '=', 'X')
    counts = {key: 0 for key in feature_attr}
    # for allele specific analysis
    ref_counts = {key: 0 for key in feature_attr}  # counter for reads assigned to the reference allele  <------------------------------------------------
    alt_counts = {key: 0 for key in feature_attr} # counter for reads assigned to the castaneus allele
    unassigned_counts = {key: 0 for key in feature_attr} # counter for reads with SNPs that could not be assigned
    # either SNP base sequence did not match expectations or no majority of SNPs for either strain

    try:
        if pe_mode:
            if ((supplementary_alignment_mode == 'ignore') and
               (secondary_alignment_mode == 'ignore')):
                primary_only = True
            else:
                primary_only = False
            if order == "name":
                read_seq = HTSeq.pair_SAM_alignments(
                        read_seq,
                        primary_only=primary_only)
            elif order == "pos":
                read_seq = HTSeq.pair_SAM_alignments_with_buffer(
                        read_seq,
                        max_buffer_size=max_buffer_size,
                        primary_only=primary_only)
            else:
                raise ValueError("Illegal order specified.")
        empty = 0
        ambiguous = 0
        notaligned = 0
        lowqual = 0
        nonunique = 0
        i = 0
        for r in read_seq:
            if i > 0 and i % 100000 == 0 and not quiet:
                sys.stderr.write(
                    "%d alignment record%s processed. %d assigned (%d ref / %d alt), %d unassigned of %d tested.\n" %
                    (i, "s" if not pe_mode else " pairs", allelicreads, refreads, altreads, unreads, totalreads))
                sys.stderr.flush()

            i += 1
            if not pe_mode:     # SINGLE END SEQUENCING <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                if not r.aligned:
                    notaligned += 1
                    write_to_samout(
                            r, "__not_aligned", samoutfile,
                            template)
                    continue
                if ((secondary_alignment_mode == 'ignore') and
                   r.not_primary_alignment):
                    continue
                if ((supplementary_alignment_mode == 'ignore') and
                   r.supplementary):
                    continue
                try:
                    if r.optional_field("NH") > 1:
                        nonunique += 1
                        write_to_samout(
                                r,
                                "__alignment_not_unique",
                                samoutfile,
                                template)
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if r.aQual < minaqual:
                    lowqual += 1
                    write_to_samout(
                            r, "__too_low_aQual", samoutfile,
                            template)
                    continue
                # here r is the alignment record - we need the read r.read.SequenceWithQualities and the list of r.cigar operations
                # the problem is that the genomic intervalls are separated from the record and stored separately in iv_seq as list
                # we need more information then iv <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
                if stranded != "reverse":
                    iv_seq = (co.ref_iv for co in r.cigar if co.type in com
                              and co.size > 0)
                else:
                    iv_seq = (invert_strand(co.ref_iv)
                              for co in r.cigar if (co.type in com and
                                                    co.size > 0))
            else:     # PAIRED END SEQUENCING
                if r[0] is not None and r[0].aligned:
                    if stranded != "reverse":
                        iv_seq = (co.ref_iv for co in r[0].cigar
                                  if co.type in com and co.size > 0)
                    else:
                        iv_seq = (invert_strand(co.ref_iv) for co in r[0].cigar
                                  if co.type in com and co.size > 0)
                else:
                    iv_seq = tuple()
                if r[1] is not None and r[1].aligned:
                    if stranded != "reverse":
                        iv_seq = itertools.chain(
                                iv_seq,
                                (invert_strand(co.ref_iv) for co in r[1].cigar
                                if co.type in com and co.size > 0))
                    else:
                        iv_seq = itertools.chain(
                                iv_seq,
                                (co.ref_iv for co in r[1].cigar
                                 if co.type in com and co.size > 0))
                else:
                    if (r[0] is None) or not (r[0].aligned):
                        write_to_samout(
                                r, "__not_aligned", samoutfile,
                                template)
                        notaligned += 1
                        continue
                if secondary_alignment_mode == 'ignore':
                    if (r[0] is not None) and r[0].not_primary_alignment:
                        continue
                    elif (r[1] is not None) and r[1].not_primary_alignment:
                        continue
                if supplementary_alignment_mode == 'ignore':
                    if (r[0] is not None) and r[0].supplementary:
                        continue
                    elif (r[1] is not None) and r[1].supplementary:
                        continue
                try:
                    if ((r[0] is not None and r[0].optional_field("NH") > 1) or
                       (r[1] is not None and r[1].optional_field("NH") > 1)):
                        nonunique += 1
                        write_to_samout(
                                r, "__alignment_not_unique", samoutfile,
                                template)
                        if multimapped_mode == 'none':
                            continue
                except KeyError:
                    pass
                if ((r[0] and r[0].aQual < minaqual) or
                   (r[1] and r[1].aQual < minaqual)):
                    lowqual += 1
                    write_to_samout(
                            r, "__too_low_aQual", samoutfile,
                            template)
                    continue

            try:
                if overlap_mode == "union":
                    fs = set()
                    for iv in iv_seq:
                        if iv.chrom not in features.chrom_vectors:
                            raise UnknownChrom
                        for iv2, fs2 in features[iv].steps():
                            fs = fs.union(fs2)
                elif overlap_mode in ("intersection-strict",
                                      "intersection-nonempty"):
                    fs = None
                    for iv in iv_seq:
                        if iv.chrom not in features.chrom_vectors:
                            raise UnknownChrom
                        for iv2, fs2 in features[iv].steps():
                            if ((len(fs2) > 0) or
                               (overlap_mode == "intersection-strict")):
                                if fs is None:
                                    fs = fs2.copy()
                                else:
                                    fs = fs.intersection(fs2)
                else:
                    sys.exit("Illegal overlap mode.")

                if fs is None or len(fs) == 0:
                    write_to_samout(
                            r, "__no_feature", samoutfile,
                            template)
                    empty += 1
                elif len(fs) > 1:
                    write_to_samout(
                            r, "__ambiguous[" + '+'.join(fs) + "]",
                            samoutfile,
                            template)
                    ambiguous += 1
                else:
                    write_to_samout(
                            r, list(fs)[0], samoutfile,
                            template)

                if fs is not None and len(fs) > 0:
                    if multimapped_mode == 'none':
                        if len(fs) == 1:
                            counts[list(fs)[0]] += 1
                            if allele_counts:
                                ref_snps = 0            # counters for allelic snps in read alignment
                                alt_snps = 0
                                unassigned = 0
                                if r.iv.chrom in SNPs.chrom_vectors:        # make sure that we have the chromosome in our SNPs array
                                    if r.iv.start > 0:
                                        test_iv = HTSeq.GenomicInterval(r.iv.chrom, r.iv.start - 1, r.iv.end + 1, r.iv.strand)
                                    else:
                                        test_iv = HTSeq.GenomicInterval(r.iv.chrom, r.iv.start, r.iv.end + 1, r.iv.strand)
                                    if len( [snp for snp_iv, snp in SNPs[test_iv].steps() if snp != 0] ) > 0:
                                        # iterating over iv_seq_iv list we catch only SNPs that are in fact overlapping read alignments
                                        totalreads += 1
                                        for op in r.cigar:
                                            if op.type in com:              # iterate over all alignment intervals that are matches
                                                if op.size > 0:             # ... and have actual bases on reference - should this always be true for (M,X,=) ?
                                                    # try to expand ref_iv to catch more CAST reads that might not be included in M and possibly also not X cigar ops
                                                    if op.query_from > 0 and op.ref_iv.start > 0:
                                                        op.ref_iv.start -= 1
                                                        op.query_from -= 1
                                                    if op.query_to < len(r.read_as_aligned.seq):
                                                        op.ref_iv.end += 1
                                                        op.query_to += 1
                                                    # ---------------------------------------------------------------------------------------------------------------
                                                    for snp_iv, snp in [ (iv,s) for iv,s in SNPs[op.ref_iv].steps() if s != 0  ]: # steps with 0 at end of snp_iv are omitted
                                                        for snp_pos in snp_iv.range(step=1):                            # iterate over all positions in snp_iv
                                                            if snp_min_qual == 0 or r.read_as_aligned.qual[op.query_from + snp_pos.start - op.ref_iv.start] >= snp_min_qual:
                                                                read_snp = check_ref_alt(chr(r.read_as_aligned.seq[op.query_from + snp_pos.start - op.ref_iv.start]).upper(), snp)
                                                                if read_snp == "ref":
                                                                    ref_snps += 1
                                                                elif read_snp == "alt":
                                                                    alt_snps += 1
                                                                else:
                                                                    print("UNASSIGNED", read_snp, bin(snp), chr(r.read_as_aligned.seq[op.query_from + snp_pos.start - op.ref_iv.start]).upper())
                                                                    unassigned += 1
                                                                # not needed - see above
                                                                # else: handle antisense strand using SequenceWithQualities.get_reverse_complement:
                                                                #     read_seq = r.read.get_reverse_complement() / check if quals need also inverted !!!!!
                                                            else:
                                                                unassigned += 1
                                                else:
                                                    print("WARNING found op.size==0; op.type =", op.type)
                                        alt_total += alt_snps
                                        ref_total += ref_snps
                                        un_total += unassigned
                                        if alt_snps > 0 and alt_snps > ref_snps:
                                            alt_counts[list(fs)[0]] += 1
                                            allelicreads += 1
                                            altreads += 1
                                            if alt_snps < max(alt_snps, ref_snps, unassigned):
                                                print("WARNING: Borderline alleleic assignment ALT:", alt_snps,"(ref:",ref_snps,", unassigned:",unassigned,")")
                                        elif ref_snps > 0 and ref_snps > alt_snps:
                                            ref_counts[list(fs)[0]] += 1
                                            allelicreads += 1
                                            refreads += 1
                                            if ref_snps < max(alt_snps, ref_snps, unassigned):
                                                print("WARNING: Borderline alleleic assignment REF:",ref_snps,"(alt:",alt_snps,", unassigned:",unassigned,")")
                                        else:
                                            unreads += 1
                                            unassigned_counts[list(fs)[0]] += 1
                                            # print("INFO: Unassigned read overlapping SNP position(s) > ref:",ref_snps,"(cast:",cast_snps,", unassigned:",unassigned,")", unreads, totalreads)
                    elif multimapped_mode == 'all':
                        for fsi in list(fs):
                            counts[fsi] += 1
                    elif multimapped_mode == 'fraction':
                        for fsi in list(fs):
                            counts[fsi] += 1.0 / len(fs)
                    elif multimapped_mode == 'random':
                        fsi = random.choice(fs)
                        counts[fsi] += 1
                    else:
                        sys.exit("Illegal multimap mode.")

            except UnknownChrom:
                write_to_samout(
                        r, "__no_feature", samoutfile,
                        template)
                empty += 1

    except:
        sys.stderr.write(
            "Error occured when processing input (%s):\n" %
            (read_seq_file.get_line_number_string()))
        raise

    if not quiet:
        sys.stderr.write(
            "%d %s processed.\n" %
            (i, "alignments " if not pe_mode else "alignment pairs"))
        sys.stderr.flush()

    if samoutfile is not None:
        samoutfile.close()

    if allele_counts:
        print("Polymorphism report: %d assigned: %d reference / %d alternate, %d alleles unassigned of %d tested.\n" % (allelicreads, refreads, altreads, unreads, totalreads))
        print("                     %d ref, %d alt, %d unassigned SNPs total" % (ref_total, alt_total, un_total))

    return {
        'isam': isam,
        'counts': counts,
        'ref_counts': ref_counts,                  # include allelic counts in result dictionary <--------------------------------------------------------------
        'alt_counts': alt_counts,
        'empty': empty,
        'ambiguous': ambiguous,
        'lowqual': lowqual,
        'notaligned': notaligned,
        'nonunique': nonunique,
    }


def count_reads_in_features(
        sam_filenames,
        gff_filename,
        order,
        max_buffer_size,
        stranded,
        overlap_mode,
        multimapped_mode,
        secondary_alignment_mode,
        supplementary_alignment_mode,
        feature_type,
        id_attribute,
        additional_attributes,
        quiet,
        minaqual,
        samouts,
        samout_format,
        output_delimiter,
        output_filename,
        output_append,
        nprocesses,
        feature_query,
        allele_counts,                           # <--------------------------allelic analysis
        SNPs,
        snp_min_qual,
        add_genome_annotation,
        gtf_file_offset,
        gtf_file_end_not_included,
        ):
    '''Count reads in features, parallelizing by file'''

    def parse_feature_query(feature_query):
        if '"' not in feature_query:
            raise ValueError('Invalid feature query')
        if '==' not in feature_query:
            raise ValueError('Invalid feature query')

        idx_quote1 = feature_query.find('"')
        idx_quote2 = feature_query.rfind('"')
        attr_name = feature_query[idx_quote1+1: idx_quote2]

        idx_equal = feature_query[:idx_quote1].find('==')
        attr_cat = feature_query[:idx_equal].strip()

        return {
            'attr_cat': attr_cat,
            'attr_name': attr_name,
            }

    if samouts != []:
        if len(samouts) != len(sam_filenames):
            raise ValueError(
                    'Select the same number of input and output files')
        # Try to open samout files early in case any of them has issues
        if samout_format in ('SAM', 'sam'):
            for samout in samouts:
                with open(samout, 'w'):
                    pass
        else:
            # We don't have a template if the input is stdin
            if (len(sam_filenames) != 1) or (sam_filenames[0] != '-'):
                for sam_filename, samout in zip(sam_filenames, samouts):
                    with pysam.AlignmentFile(sam_filename, 'r') as sf:
                        with pysam.AlignmentFile(samout, 'w', template=sf):
                            pass
    else:
        samouts = [None for x in sam_filenames]

    # Try to open samfiles to fail early in case any of them is not there
    if (len(sam_filenames) != 1) or (sam_filenames[0] != '-'):
        for sam_filename in sam_filenames:
            with pysam.AlignmentFile(sam_filename, 'r') as sf:
                pass

    if feature_query is not None:
        feature_qdic = parse_feature_query(feature_query)
    features = HTSeq.GenomicArrayOfSets("auto", stranded != "no") # ----- features GENOME ANNOTATION of exons ------ <<<<<<<<<<<<<<<<
    if gtf_file_end_not_included:
        gff = HTSeq.GFF_Reader(gff_filename, end_included=False)
    else:
        gff = HTSeq.GFF_Reader(gff_filename)
    attributes = {}
    gene_info = {}   # for annotation of genes -----------------------<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    i = 0
    try:
        for f in gff:
            if f.type == feature_type:
                try:
                    feature_id = f.attr[id_attribute]
                except KeyError:
                    raise ValueError(
                            "Feature %s does not contain a '%s' attribute" %
                            (f.name, id_attribute))
                if stranded != "no" and f.iv.strand == ".":
                    raise ValueError(
                            "Feature %s at %s does not have strand information but you are "
                            "running htseq-count in stranded mode. Use '--stranded=no'." %
                            (f.name, f.iv))

                if feature_query is not None:
                    # Skip the features that don't even have the right attr
                    if feature_qdic['attr_cat'] not in f.attr:
                        continue
                    # Skip the ones with an attribute with a different name
                    # from the query (e.g. other genes)
                    if f.attr[feature_qdic['attr_cat']] != feature_qdic['attr_name']:
                        continue
                if gtf_file_offset == 1:
                    features[f.iv] += feature_id
                else:
                    offset = 1 - gtf_file_offset
                    if f.iv.start >= offset:
                        f.iv.start = f.iv.start - offset
                    if f.iv.end >= offset:
                        f.iv.end = f.iv.end - offset
                attributes[feature_id] = [f.attr[attr] if attr in f.attr else '' for attr in additional_attributes]
                if "chrom" in add_genome_annotation or "gene_length" in add_genome_annotation:
                    # UPDATE START AND END OF GENE
                    if feature_id in gene_info.keys():
                        if gene_info[feature_id].chrom == f.iv.chrom:
                            if gene_info[feature_id].strand == f.iv.strand:
                                if gene_info[feature_id].start > f.iv.start:
                                    gene_info[feature_id].start = f.iv.start
                                if gene_info[feature_id].end < f.iv.end:
                                    gene_info[feature_id].end = f.iv.end
                            else:
                                print("Ambiguous strand information for", feature_id, "encountered. Discarding exon information.")
                        else:
                            print("Ambibugous chromosome definition for", feature_id, "encountered. Discarding exon information.")
                    else:
                        gene_info[feature_id] = f.iv
            i += 1
            if i % 100000 == 0 and not quiet:
                sys.stderr.write("%d GFF lines processed.\n" % i)
                sys.stderr.flush()
        # CALCULATE EXON LENGTH HERE !!!! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        if "exon_length" in add_genome_annotation:
            exon_length = {}
            for iv,feature_ids in features.steps():
                for feature_id in feature_ids:
                    if feature_id in exon_length.keys():
                        exon_lenth[feature_id] += iv.end-iv.start
                    else:
                        exon_length[feature_id] = iv.end-iv.start
        # CALCULATE GENE LENGTH FROM START AND END AND STORE IN ATTRIBUTES
        if len(add_genome_annotation) > 0:
            for feature_id in attributes.keys():
                annotation = []
                if "chrom" in add_genome_annotation:
                    annotation = [gene_info[feature_id].chrom, gene_info[feature_id].start, gene_info[feature_id].end, gene_info[feature_id].strand]
                if "gene_length" in add_genome_annotation:
                    gene_length = gene_info[feature_id].end - gene_info[feature_id].start
                    annotation.append(gene_length)
                if "exon_length" in add_genome_annotation:
                    annotation.append(exon_length[feature_id])
                attributes[feature_id] = annotation + attributes[feature_id]

    except:
        sys.stderr.write(
            "Error occured when processing GFF file (%s):\n" %
            gff.get_line_number_string())
        raise

    feature_attr = sorted(attributes.keys())

    if not quiet:
        sys.stderr.write("%d GFF lines processed.\n" % i)
        sys.stderr.flush()

    if len(feature_attr) == 0:
        sys.stderr.write(
            "Warning: No features of type '%s' found.\n" % feature_type)

    # Prepare arguments for counting function
    args = []
    for isam, (sam_filename, samout_filename) in enumerate(zip(sam_filenames, samouts)):
        args.append((
            isam,
            sam_filename,
            features,
            feature_attr,
            order,
            max_buffer_size,
            stranded,
            overlap_mode,
            multimapped_mode,
            secondary_alignment_mode,
            supplementary_alignment_mode,
            feature_type,
            id_attribute,
            additional_attributes,
            quiet,
            minaqual,
            samout_format,
            samout_filename,
            allele_counts,
            SNPs,
            snp_min_qual,
            ))

    # Count reads
    if nprocesses > 1:
        with multiprocessing.Pool(nprocesses) as pool:
            results = pool.starmap(count_reads_single_file, args)
        results.sort(key=operator.itemgetter('isam'))
    else:
        results = list(itertools.starmap(count_reads_single_file, args))

    # Write output
    other_features = [
        ('__no_feature', 'empty'),
        ('__ambiguous', 'ambiguous'),
        ('__too_low_aQual', 'lowqual'),
        ('__not_aligned', 'notaligned'),
        ('__alignment_not_unique', 'nonunique'),
        ]
    pad = ['' for attr in additional_attributes]
    for ifn, fn in enumerate(feature_attr):
        fields = [fn] + attributes[fn] + [str(r['counts'][fn]) for r in results]
        # for allelic analysis <---------------------------------------------------------- include allele counts in output
        fields = fields + [str(r['ref_counts'][fn]) for r in results]
        fields = fields + [str(r['alt_counts'][fn]) for r in results]
        line = output_delimiter.join(fields)
        if output_filename == '':
            print(line)
        else:
            omode = 'a' if output_append or (ifn > 0) else 'w'
            with open(output_filename, omode) as f:
                f.write(line)
                f.write('\n')

    for title, fn in other_features:
        fields = [title] + pad + [str(r[fn]) for r in results]
        line = output_delimiter.join(fields)
        if output_filename == '':
            print(line)
        else:
            with open(output_filename, 'a') as f:
                f.write(line)
                f.write('\n')

def load_genome_fasta(fasta_filename):
    genome=HTSeq.FastaReader(fasta_filename)
    chromosome={ str(i.name): bytearray(i.seq) for i in genome }
    return chromosome

def write_genome_fasta(chromosome, fasta_filename):
    n = 0
    try:
        with open(fasta_filename, "w") as fasta_file:
            for chrom in chromosome.keys():
                s = HTSeq.Sequence(bytes(chromosome[chrom]), chrom)
                s.write_to_fasta_file(fasta_file)
                n+=1
        print(n, "sequence records written to", fasta_filename)
    except IOError:
        print("Error when attempting to write", n, "sequence record to", fasta_filename)

def vcf_direct(vcf_filename, chromosome):
    f=open(VCF_FILENAME, "r")
    for r in f:
        if r[0]!="#":
            break    # get only first data line
    print("from vcf file directly",r.split('\t')[:6])
    print("Genome pos > Chr", r.split("\t")[0], r.split("\t")[1], ":", chr(chromosome[r.split("\t")[0]][int(r.split("\t")[1])]))
    print("Genome pos - 1 > Chr", r.split("\t")[0], int(r.split("\t")[1])-1, ":", chr(chromosome[r.split("\t")[0]][int(r.split("\t")[1])-1]))

def check_vcf_file_offset(vcf_filename, chromosome, vcf_file_offset=1, test_num=10000):
    vcf_file = HTSeq.VCF_Reader(vcf_filename)
    vcf_file.parse_meta()
    vcf_file.make_info_dict()
    zero_based = 0
    one_based = 0
    offset_based = 0
    n = 0
    for morphism in vcf_file:
        if(chr(chromosome[morphism.pos.chrom][morphism.pos.start]) == morphism.ref):
            zero_based += 1
        if morphism.pos.start > 0:
            if(chr(chromosome[morphism.pos.chrom][morphism.pos.start-1]) == morphism.ref):
                one_based += 1
        if morphism.pos.start >= vcf_file_offset:
            if(chr(chromosome[morphism.pos.chrom][morphism.pos.start-vcf_file_offset]) == morphism.ref):
                offset_based += 1
        n += 1
        if zero_based + one_based + offset_based > test_num:
            break
        else:
            if zero_based + one_based + offset_based < 20:
                print("VCF Position > Chr", morphism.pos.chrom, morphism.pos.start, ":", chr(chromosome[morphism.pos.chrom][morphism.pos.start]))
                if morphism.pos.start > 0:
                    print("VCF Position -1 > Chr", morphism.pos.chrom, morphism.pos.start-1, ":", chr(chromosome[morphism.pos.chrom][morphism.pos.start-1]))
                if vcf_file_offset > 1 and morphism.pos.start >= vcf_file_offset:
                    print("VCF Position", -vcf_file_offset,"> Chr", morphism.pos.chrom, morphism.pos.start-vcf_file_offset, ":", chr(chromosome[morphism.pos.chrom][morphism.pos.start-vcf_file_offset]))
                print("Ref allele =", morphism.ref)
    print("Of", n, "tests performed:")
    print("0-based      :", zero_based)
    print("1-based      :", one_based)
    if vcf_file_offset > 1:
        print("offset-based :", offset_based, "(VCF file offset", vcf_file_offset)
    if offset_based == max(zero_based, one_based, offset_based):
        print("VCF file offset is", vcf_file_offset, ": use --vcf-file-offset=%d" % vcf_file_offset)
        return vcf_file_offset
    elif zero_based > one_based:
        print("VCF Reader generates 0-based positions: use --vcf-file-offset=0")
        return 0
    else:
        print("VCF Reader generates 1-based positions: use --vcf-file-offset=1")
        return 1

def check_gtf_file_stop_codon(gtf_filename, chromosome, test_num=10000):
    zero_based = 0
    one_based = 0
    end_included = 0
    end_not_included = 0
    n = 0
    print("Checking GTF annotations of stop codons")
    gtf_file = HTSeq.GFF_Reader(gtf_filename, end_included=True)
    CDS_stop = ["TAA", "TAG", "TGA"] # UAA, UAG and UGA
    for feature in gtf_file:
        if feature.type=="stop_codon" and feature.iv.strand=="+":
            if feature.iv.end-feature.iv.start == 3:
                end_included += 1
            if feature.iv.end-feature.iv.start == 2:
                end_not_included += 1
            stop_codon = "".join([chr(chromosome[feature.iv.chrom][feature.iv.start]), chr(chromosome[feature.iv.chrom][feature.iv.start+1]), chr(chromosome[feature.iv.chrom][feature.iv.start+2])])
            if stop_codon in CDS_stop:
                zero_based += 1
            stop_codon = "".join([chr(chromosome[feature.iv.chrom][feature.iv.start-1]), chr(chromosome[feature.iv.chrom][feature.iv.start]), chr(chromosome[feature.iv.chrom][feature.iv.start+1])])
            if stop_codon in CDS_stop:
                one_based += 1
            if end_included + end_not_included > test_num and zero_based + one_based > test_num:
                break
            n += 1
    print("Of", n, "tests performed:")
    print("end_included:", end_included)
    print("end_not_included", end_not_included)
    if end_included > end_not_included:
        print("The end is included in the intervals of the GTF file as HTSeq expects as default, everything is fine.")
        print("No action is needed (HTSeq.GFF_Reader( ... end_included=True) is the default.")
        end_included_adjustment = 0
        GFF_READER_END_INCLUDED = True
    else:
        print("WARNING:")
        print("The end is not included in the intervals of the GTF file as HTSeq expects as default.")
        print("action is needed (HTSeq.GFF_Reader( ... end_included=True) is the default.")
        print("Do this by specifying HTSeq.GFF_Reader( ... end_included=False) when using GFF_Reader!")
        end_included_adjustment = 1
        GFF_READER_END_INCLUDED = False
    print("0-based:", zero_based)
    print("1-based", one_based)
    if zero_based > one_based:
        print("GFF Reader generates 0-based positions, everything is fine.")
        GFF_READER_ZERO_BASED = True
    else:
        print("WARNING:")
        print("GFF Reader gernerates 1-based positions that need to be converted to 0-based positions to be consistent witH HTSeq and FastaReader!")
        print("Do this by subtraction 1 from the pos.start of GFF genomic positions!")
        GFF_READER_ZERO_BASED = False
    if n > (zero_based + one_based) * 1.1:
        print("WARNING:")
        print("A large number of stop codons are not matching the consensus. This might indicate a mismatch between the annotation and genome reference.")
        print("Check the genome assembly of the GTF file and the genome FASTA file are corresponding.")
    return GFF_READER_END_INCLUDED, GFF_READER_ZERO_BASED

def check_gtf_file_intron_junctions(gtf_filename, chromosome, end_included=True, test_num=10000):
    gtf_file = HTSeq.GFF_Reader(gtf_filename, end_included=end_included)
    print("Checking intron flanking sequences based on exon annotation in GTF file", gtf_filename)
    print("\texon-end]GT..intron..AG[exon-start")
    zero_based = 0
    one_based = 0
    n = 0
    for feature in gtf_file:
        if feature.type=="exon" and feature.iv.strand=="+":
            intron_end = "".join([chr(chromosome[feature.iv.chrom][feature.iv.start-2]), chr(chromosome[feature.iv.chrom][feature.iv.start-1])])
            intron_start = "".join([chr(chromosome[feature.iv.chrom][feature.iv.end]), chr(chromosome[feature.iv.chrom][feature.iv.end+1])])
            if intron_start == "GT":
                zero_based += 1
            if intron_end == "AG":
                zero_based += 1
            intron_end = "".join([chr(chromosome[feature.iv.chrom][feature.iv.start-3]), chr(chromosome[feature.iv.chrom][feature.iv.start-2])])
            intron_start = "".join([chr(chromosome[feature.iv.chrom][feature.iv.end-1]), chr(chromosome[feature.iv.chrom][feature.iv.end])])
            if intron_start == "GT":
                one_based += 1
            if intron_end == "AG":
                one_based += 1
            if n < 20:
                print("start pos-2", feature.iv.start-2,":","".join([chr(chromosome[feature.iv.chrom][feature.iv.start-2]), chr(chromosome[feature.iv.chrom][feature.iv.start-1]), "[", chr(chromosome[feature.iv.chrom][feature.iv.start])]))
                print("end pos-1", feature.iv.end-1,":","".join([chr(chromosome[feature.iv.chrom][feature.iv.end-1]), "]", chr(chromosome[feature.iv.chrom][feature.iv.end]), chr(chromosome[feature.iv.chrom][feature.iv.end+1])]))
            n += 1
            if n > test_num:
                break
    print("Of", n, "tests performed:")
    print("0-based:", zero_based)
    print("1-based", one_based)
    if zero_based > one_based:
        print("GFF Reader generates 0-based positions, everything is fine.")
        GFF_READER_ZERO_BASED = True
    else:
        print("WARNING:")
        print("GFF Reader gernerates 1-based positions that need to be converted to 0-based positions to be consistent witH HTSeq and FastaReader!")
        print("Do this by subtraction 1 from the pos.start of GFF genomic positions!")
        GFF_READER_ZERO_BASED = False
    if zero_based + one_based < n:
        print("WARNING:")
        print("A large number of slice sites are not matching the consensus. This might indicate a mismatch between the annotation and genome reference.")
        print("Check the genome assembly of the GTF file and the genome FASTA file are corresponding.")        
    return GFF_READER_ZERO_BASED

def check_bam_file(bam_filename, chromosome, test_num=10000):
    bam_file = HTSeq.BAM_Reader(bam_filename)
    zero_based = 0
    one_based = 0
    n = 0
    for r in bam_file:
        if r.aligned:
            if r.iv.strand == "+":
                for op in r.cigar:
                    if op.type == "M":
                        if n < 20:
                            print("Read   >", r.read.seq[op.query_from:op.query_to])
                            print("Genome >", chromosome[op.ref_iv.chrom][op.ref_iv.start:op.ref_iv.end])
                        if r.read.seq[op.query_from:op.query_to] == chromosome[op.ref_iv.chrom][op.ref_iv.start:op.ref_iv.end]:
                            zero_based += 1
                        if r.read.seq[op.query_from:op.query_to] == chromosome[op.ref_iv.chrom][op.ref_iv.start-1:op.ref_iv.end-1]:
                            one_based += 1
                n += 1
                if n > test_num:
                    break
    print("Of", n, "tests performed:")
    print("0-based:", zero_based)
    print("1-based", one_based)
    if zero_based > one_based:
        print("BAM Reader generates 0-based positions, everything is fine.")
        BAM_READER_ZERO_BASED = True
    else:
        print("WARNING:")
        print("BAM Reader gernerates 1-based positions that need to be converted to 0-based positions to be consistent witH HTSeq and FastaReader!")
        print("Do this by subtraction 1 from the pos.start of BAM genomic positions!")
        warned += 1
        BAM_READER_ZERO_BASED = False
    return BAM_READER_ZERO_BASED

def consistency_check(fasta_filename, vcf_filename, vcf_file_offset, gtf_filename, bam_filenames, test_num=10000):
    print("PERFORMING CONSISTENCY CHECK")
    print("Reading genome fasta file", fasta_filename)
    if os.path.exists(fasta_filename):
        chromosome=load_genome_fasta(fasta_filename)
        for c in chromosome.keys():
            print("Chromosome", c, ":", len(chromosome[c]), "basepairs")
    else:
        print("FASTA file for reference genome not found. Use --genome-fasta-filename to specify the path to the FASTA file with the reference genome.")
        sys.exit()
    if vcf_filename is not None and os.path.exists(vcf_filename):
        print()
        print("Checking VCF position index", vcf_filename)
        vcf_direct(vcf_filename, chromosome)
        VCF_FILE_OFFSET = check_vcf_file_offset(vcf_filename, chromosome, vcf_file_offset=vcf_file_offset, test_num=test_num)
    if gtf_filename is not None and os.path.exists(gtf_filename):
        print()
        print("Checking annotation position index in GTF file", gtf_filename)
        GFF_READER_END_INCLUDED, GFF_READER_ZERO_BASED = check_gtf_file_stop_codon(gtf_filename, chromosome, test_num=test_num)
        print()
        if check_gtf_file_intron_junctions(gtf_filename, chromosome, end_included=GFF_READER_END_INCLUDED, test_num=test_num) != GFF_READER_ZERO_BASED:
            print("WARNING: Inconsitent evidence for genomic positions in GTF file are 0-based or 1-based")
    if len(bam_filenames) > 0:
        print()
        print("Checking genome and bam file index")
        for bam_filename in bam_filenames:
            print(bam_filename)
            if os.path.exists(bam_filename):
                check_bam_file(bam_filename, chromosome, test_num=test_num)

# ==============================================================================================================================================================================================================================================================

def nmask_genome(fasta_filename, nmask_filename, vcf_filename, ref_samples, alt_samples, omit_ref_genome=False, vcf_file_offset=1, silent=False):  
    # nmask_genome(args.genome_fasta_filename, args.genome_nmask_filename, args.vcf_filename, args.ref_samples, args.alt_samples, omit_ref_genome=args.omit_ref_genome, vcf_file_offset=args.vcf_file_offset, silent=args.quiet)
    print("N-mask genome")
    print("Reading genome fasta file", fasta_filename)
    if os.path.exists(fasta_filename):
        chromosome=load_genome_fasta(fasta_filename)
        for c in chromosome.keys():
            print("Chromosome", c, ":", len(chromosome[c]), "basepairs")
    else:
        print("FASTA file for reference genome not found. Use --genome-fasta-filename to specify the path to the FASTA file with the reference genome.")
        sys.exit()
    if vcf_filename is not None and os.path.exists(vcf_filename):
        with open(vcf_filename) as vcf_file:
            if not silent:
                print("Reading polymorphism information from VCF file:", vcf_filename)
                print("Using genomic position offset", vcf_file_offset)
            gt_idx=0
            fi_idx=-1
            n=0
            bases_changed=0
            bases_N=0
            for line in vcf_file:
                line=line.strip()
                if len(line) < 2:
                    continue
                elif line[0]=="#":
                    if line[1]=="#":            # meta lines
                        if not silent:
                            print(line)
                    else:                       # header line
                        header={e:n for n,e in enumerate(line[1:].split("\t"))} # dictionary field name to index
                        if not silent:
                            print("HEADER -> data field index\n","\n".join([k+":"+str(header[k]) for k in header.keys()]))
                else:                           # data lines
                    data_fields=line.split("\t")
                    if data_fields[header["FORMAT"]][:2]!="GT" or data_fields[header["FORMAT"]][-2:]!="FI":
                        frmt_fields = data_fields[header["FORMAT"]].split(":")
                        gt_idx = frmt_fields.index("GT")
                        fi_idx = frmt_fields.index("FI")
                    allele_sequences = [ data_fields[header["REF"]] ] + data_fields[header["ALT"]].split(",")
                    ref_alleles=set([allele_sequences[int(gl[0])] for gl in [ sample_data[gt_idx].replace("|","/").split("/") for sample_data in [ data_fields[header[sample]].split(":") for sample in ref_samples ] if sample_data[fi_idx]=="1" ] if len(set(gl))==1 and gl[0]!="." and len(allele_sequences[int(gl[0])]) == 1])
                    if not omit_ref_genome and data_fields[header["FILTER"]] == "PASS" and len(data_fields[header["REF"]]) == 1:
                        ref_alleles |= { data_fields[header["REF"]] }
                    alt_alleles=set([allele_sequences[int(gl[0])] for gl in [ sample_data[gt_idx].replace("|","/").split("/") for sample_data in [ data_fields[header[sample]].split(":") for sample in alt_samples ] if sample_data[fi_idx]=="1" ] if len(set(gl))==1 and gl[0]!="." and len(allele_sequences[int(gl[0])]) == 1])
                    # here we need to mask the genome with ref_alleles and alt_alleles containing sets with the bases
                    # first if there is any difference between ref and alt or ref and the reference genome then the base goes to N
                    # check if the ref_allele is actually correct in the genome else write a warning
                    chrom = data_fields[header["CHROM"]]
                    pos = int(data_fields[header["POS"]]) - vcf_file_offset
                    vcf_ref = data_fields[header["REF"]]
                    if len(vcf_ref)==1 and all([len(e)==1 for e in list(ref_alleles)+list(alt_alleles)]):                    # only use SNPs
                        genome_ref = chr(chromosome[chrom][pos])
                        """if genome_ref != vcf_ref:
                            print("Warning: Chrom %s position %d VCF ref <%s> mismatches genome reference chrom %s position %d <%s>, %d bases changed. %d N-masked" % (chrom, pos+vcf_file_offset, vcf_ref, chrom, pos, genome_ref, bases_changed, bases_N) )
                        else:
                            print("Chrom %s pos %d VCF natches genome ref <%s>." % (chrom, pos, genome_ref) )
                        # if vcf ref and alt alleles are the same and different from genome ref then and only then use the vcf base for the genome position."""
                        if len(ref_alleles)==1 and ref_alleles == alt_alleles:
                            # ref and alt allele are the very same base - so change the genome to this base if different from the genome reference base
                            new_base = ref_alleles.pop()
                            if genome_ref != new_base:
                                print("CHR %s pos %d relacing %s with %s", chrom, pos, genome_ref, new_base)
                                chromosome[chrom][pos] = ord(new_base)
                                bases_changes+=1
                            # else ref and alt allele are the same as the genome reference base and nothing should be done!
                        else:
                            chromosome[chrom][pos] = ord(b"N")  # we use bytearray but this works with bytes strings: chromosomes[chrom] = chromosomes[chrom][:pos] + b"N" + chromosomes[chrom][pos+1:]
                            bases_N+=1
                    n+=1
                    if not silent:
                        if n % 1000000 == 0:
                            print(n, "vcf records processed,", bases_changed, "bases changed,", bases_N, "N-masked.")
        if not silent:
            print(n,"vcf records processed.")
            print(bases_changed, "bases changed in reference genome.")
            print(bases_N, "bases N-masked in reference genome.")
            print(bases_changed+bases_N, "total changes made to reference genome.")
        write_genome_fasta(chromosome, nmask_filename)
    else:
        print("VCF file for polymorphisms not found. Use --SNP-vcf-filename to specify the path to the VCF file with the SNP information.")
        sys.exit()

# ==============================================================================================================================================================================================================================================================

def my_showwarning(message, category, filename, lineno=None, file=None,
                   line=None):
    sys.stderr.write("Warning: %s\n" % message)

def main():

    pa = argparse.ArgumentParser(
        usage="%(prog)s [options] alignment_file gff_file",
        description="This script takes one or more alignment files in SAM/BAM " +
        "format and a feature file in GFF format and calculates for each feature " +
        "the number of reads mapping to it. See " +
        "http://htseq.readthedocs.io/en/master/count.html for details.",
        epilog="Written by Simon Anders (sanders@fs.tum.de), " +
        "European Molecular Biology Laboratory (EMBL) and Fabio Zanini " +
        "(fabio.zanini@unsw.edu.au), UNSW Sydney. (c) 2010-2020. " +
        "Released under the terms of the GNU General Public License v3. " +
        "Part of the 'HTSeq' framework, version %s." % HTSeq.__version__)

    pa.add_argument(
            "samfilenames", nargs='+', type=str,
            help="Path to the SAM/BAM files containing the mapped reads. " +
            "If '-' is selected, read from standard input")

    pa.add_argument(
            "featuresfilename", type=str,
            help="Path to the GTF file containing the features")

    pa.add_argument(
            "-f", "--format", dest="samtype",
            choices=("sam", "bam", "auto"), default="auto",
            help="Type of <alignment_file> data. DEPRECATED: " +
            "file format is detected automatically. This option is ignored.")

    pa.add_argument(
            "-r", "--order", dest="order",
            choices=("pos", "name"), default="name",
            help="'pos' or 'name'. Sorting order of <alignment_file> (default: name). Paired-end sequencing " +
            "data must be sorted either by position or by read name, and the sorting order " +
            "must be specified. Ignored for single-end data.")

    pa.add_argument(
            "--max-reads-in-buffer", dest="max_buffer_size", type=int,
            default=30000000,
            help="When <alignment_file> is paired end sorted by position, " +
            "allow only so many reads to stay in memory until the mates are " +
            "found (raising this number will use more memory). Has no effect " +
            "for single end or paired end sorted by name")

    pa.add_argument(
            "-s", "--stranded", dest="stranded",
            choices=("yes", "no", "reverse"), default="yes",
            help="Whether the data is from a strand-specific assay. Specify 'yes', " +
            "'no', or 'reverse' (default: yes). " +
            "'reverse' means 'yes' with reversed strand interpretation")

    pa.add_argument(
            "-a", "--minaqual", type=int, dest="minaqual",
            default=10,
            help="Skip all reads with MAPQ alignment quality lower than the given " +
            "minimum value (default: 10). MAPQ is the 5th column of a SAM/BAM " +
            "file and its usage depends on the software used to map the reads.")

    pa.add_argument(
            "-t", "--type", type=str, dest="featuretype",
            default="exon",
            help="Feature type (3rd column in GTF file) to be used, " +
            "all features of other type are ignored (default, suitable for Ensembl " +
            "GTF files: exon)")

    pa.add_argument(
            "-i", "--idattr", type=str, dest="idattr",
            default="gene_id",
            help="GTF attribute to be used as feature ID (default, " +
            "suitable for Ensembl GTF files: gene_id). All feature of the " +
            "right type (see -t option) within the same GTF attribute will " +
            "be added together. The typical way of using this option is to " +
            "count all exonic reads from each gene and add the exons " +
            "but other uses are possible as well.")

    pa.add_argument(
            "--additional-attr", type=str,
            action='append',
            default=[],
            help="Additional feature attributes (default: none, " +
            "suitable for Ensembl GTF files: gene_name). Use multiple times " +
            "for more than one additional attribute. These attributes are " +
            "only used as annotations in the output, while the determination " +
            "of how the counts are added together is done based on option -i.")

    pa.add_argument(
            "-m", "--mode", dest="mode",
            choices=("union", "intersection-strict", "intersection-nonempty"),
            default="union",
            help="Mode to handle reads overlapping more than one feature " +
            "(choices: union, intersection-strict, intersection-nonempty; default: union)")

    pa.add_argument(
            "--nonunique", dest="nonunique", type=str,
            choices=("none", "all", "fraction", "random"), default="none",
            help="Whether and how to score reads that are not uniquely aligned " +
            "or ambiguously assigned to features " +
            "(choices: none, all, fraction, random; default: none)")

    pa.add_argument(
            "--secondary-alignments", dest="secondary_alignments", type=str,
            choices=("score", "ignore"), default="ignore",
            help="Whether to score secondary alignments (0x100 flag)")

    pa.add_argument(
            "--supplementary-alignments", dest="supplementary_alignments", type=str,
            choices=("score", "ignore"), default="ignore",
            help="Whether to score supplementary alignments (0x800 flag)")

    pa.add_argument(
            "-o", "--samout", type=str, dest="samouts",
            action='append',
            default=[],
            help="Write out all SAM alignment records into " +
            "SAM/BAM files (one per input file needed), annotating each line " +
            "with its feature assignment (as an optional field with tag 'XF')" +
            ". See the -p option to use BAM instead of SAM.")

    pa.add_argument(
            "-p", '--samout-format', type=str, dest='samout_format',
            choices=('SAM', 'BAM', 'sam', 'bam'), default='SAM',
            help="Format to use with the --samout option."
            )

    pa.add_argument(
            "-d", '--delimiter', type=str, dest='output_delimiter',
            default='\t',
            help="Column delimiter in output (default: TAB)."
            )
    pa.add_argument(
            "-c", '--counts_output', type=str, dest='output_filename',
            default='',
            help="Filename to output the counts to instead of stdout."
            )

    pa.add_argument(
            '--append-output', action='store_true', dest='output_append',
            help='Append counts output. This option is useful if you have ' +
            'already creates a TSV/CSV/similar file with a header for your ' +
            'samples (with additional columns for the feature name and any ' +
            'additionl attributes) and want to fill in the rest of the file.'
            )

    pa.add_argument(
            "-n", '--nprocesses', type=int, dest='nprocesses',
            default=1,
            help="Number of parallel CPU processes to use (default: 1)."
            )

    pa.add_argument(
            '--feature-query', type=str, dest='feature_query',
            default=None,
            help='Restrict to features descibed in this expression. Currently ' +
            'supports a single kind of expression: attribute == "one attr" to ' +
            'restrict the GFF to a single gene or transcript, e.g. ' +
            '--feature-query \'gene_name == "ACTB"\' - notice the single ' +
            'quotes around the argument of this option and the double ' +
            'quotes around the gene name. Broader queries might become ' +
            'available in the future.',
            )

    pa.add_argument(
            "-q", "--quiet", action="store_true", dest="quiet",
            help="Suppress progress report")  # and warnings" )

    pa.add_argument(
            "--version", action="store_true",
            help='Show software version and exit')

    pa.add_argument(
            "-z", "--SNP-allelic-analysis", action="store_true", dest="allele_counts",
            help="Perform allele specific analysis")                                         # <----- allelic analysis

    pa.add_argument(
            "--add-chromosome-information", action="store_true", dest="add_chromosome_information",
            help="include chromosome information (chromosome, start, and end genomic position, and strand) in output")

    pa.add_argument(
            "--add-genelength-information", action="store_true", dest="add_genelength_information",
            help="include gene locus length information in output (basepairs from start of first to end of last exon)")

    pa.add_argument(
            "--add-exonlength-information", action="store_true", dest="add_exonlength_information",
            help="include information of the length of all exons of a gene superimposed in output for calculating RKM and TPM")

    pa.add_argument(
            "--SNP-min-base-qual", type=int, dest="snp_min_qual",
            default=0,
            help="Select the minimum quality of the read sequence at the position of the polymorphic basepair." +
            "Default is --SNP-min-base-qual=0 and will accept any quality. If a higher quality is specified " +
            "than the read base quality the SNP will not be used for allelic count. The total (non allelic " +
            "read count is unaffected by this option."
            )

    pa.add_argument(
            "--SNP-omit-ref-genome", action="store_true", dest="omit_ref_genome",
            help="Do not include the reference genome sequence as a reference sample")

    pa.add_argument(
            "--SNP-force-biallelic", action="store_true", dest="force_biallelic",
            help="Only consider biallelic SNPs, discard SNPs with multiallelic reference or alternate samples")

    pa.add_argument(
            "--SNP-fast-vcf-reader", action="store_true", dest="use_fast_vcf",
            help="Read VCF file directly using a fast implementation.")

    pa.add_argument(
            '--SNP-vcf-filename', type=str, dest='vcf_filename',
            default='',
            help="Filename for VCF file with SNP information."
            )
    
    pa.add_argument(
            '--SNP-bed-filename', type=str, dest='bed_filename',
            default='',
            help="Filename for BED file with SNP information."
            )

    pa.add_argument(
            "--SNP-use-bed-file-cache", action="store_true", dest="use_bed_file_cache",
            help="Generate and use a bed file cache for SNP information. Filename must be given using --SNP_bed_filename.")

    pa.add_argument(
            "--SNP-force-generate-bed-file-cache", action="store_true", dest="force_generate_bed_file_cache",
            help="Generate a new bed file cache for SNP information, If the file already exists overwrite. Filename must be given using --SNP_bed_filename.")

    pa.add_argument(
            "--SNP-ref-sample", type=str, dest='ref_samples',
            action='append',
            default=[],
            help="Name of samples in VCF file that are considered ref. " +
            "Use multiple times for more than one reference sample. " +
            "Note that the exact spelling of sample names as in the VCF file is required. " +
            "Filename must be given using --SNP_vcf_filename. ")

    pa.add_argument(
            "--SNP-alt-sample", type=str, dest='alt_samples',
            action='append',
            default=[],
            help="Name of samples in VCF file that are considered alt. " +
            "Use multiple times for more than one alternate sample. " +
            "Note that the exact spelling of sample names as in the VCF file is required. " +
            "Filename must be given using --SNP_vcf_filename. ")

    pa.add_argument(
            "--SNP-list-vcf-file-samples", action="store_true", dest="list_vcf_file_samples",
            help="List the names of all samples in the VCF file and exit. Filename must be given using --SNP_vcf_filename.")

    pa.add_argument(
            "--SNP-vcf-file-offset", type=int, dest="vcf_file_offset",
            default=1,
            help="Force a specific offset for the genomic position of the first basepair for the VCF file " +
            "specified with --SNP_vcf_filename. The default is 1 (--SNP-vcf_file-offset=1). If in doubt use " +
            "--check-consistency for perfoming a validation of genomic coordinates for your VCF file."
            )

    pa.add_argument(
            "--check-consistency", action="store_true", dest="consistency_check",
            help="Perform a verification of genome coordinates for your annotation GTF, VCF, and SAM/BAM files. " +
            "Filenames must be provided as arguments, whereby the last argument is used as the annotation GTF file. " +
            "Use --SNP-vcf-filename to provide a VCF file with polyomorphism information. In case a VCF file is provided the reference allele " +
            "will be checked against the reference genome FASTA, --SNP-ref-sample and --SNP-omit-ref-genome are ignored. " +
            "For GTF and BAM files reasonable validation of the genomic corrdinates between the reference genome FASTA and " +
            "the files is performed. Warnings will be shown in case suspected errors are identified and solutions are suggested. " +
            "It is recommended to perform this check for allelic analysis with all filetyes as a mismatch in the genomic" +
            "position of just one basepair will miss any SNPs!\n" +
            "From the definition of different filetypes the following can be expected:\n"+
            "BAM files have a 0-based genomic position index\n" +
            "SAM files have a 1-based genomic position index\n" +
            "GTF files have a 1-based genomic position index\n" +
            "BED files have a 0-based genomic position index\n" +
            "VCF files have a 1-based genomic position index\n" +
            "BCF files have a 0-based genomic position index\n" +
            "FASTA files have a 0-based genomic position index\n" +
            "FASTQ files have a 0-based genomic position index\n" +
            "These files have positions for a specific genome assembly of a specific species. Erronous mixing genome assemblies," +
            "(or species) obviously will lead to useless results, but might occur due to filename mistakes or confusion between " +
            "genome assemblies from different providers. The latter is a real concern as often prebuilt indices are used for " +
            "NGS aligners. Using --check-consistency an independent verification for correctness of the analysis strategy can be obtained.")

    pa.add_argument(
            '--genome-fasta-filename', type=str, dest='genome_fasta_filename',
            default='',
            help="Filename of a FASTA file with the reference genome. Note that chromosome names and genomic positions must match with that of BAM, GTF, and VCF files. " +
            "This file is  used if a consistency check is performed using option --check-consistency or if an N-masked genome is generated with --genome-nmask."
            )

    pa.add_argument(
            "--num-consistency-tests", type=int, dest="test_num",
            default=10000,
            help="Maximal number of tests sampled from each file for performing a consistency test. " +
            "The default is 10000. If a large number is specified the actual tests might be limited by the data. " +
            "Use with --check-consistency for validation of genomic coordinates for your VCF, GTF, and BAM files."
            )

    # N-MASK GENOME IMPLEMENTED ======================================================================================
    
    pa.add_argument(
            "--genome-nmask", action="store_true", dest="genome_nmask",
            help="Produce an N-masked genome for the input genme FASTA file given by --genome-fasta-filename usinf either " +
            "a VCF file given by --SNP-vcf-filename or a BED file given by --SNP-bed-filename. The resulting N-masked " +
            "genome is stored output to a FASTA file given by --genome-nmask-filename. " +
            "Note that chromosome names and genomic positions must match in the genome FASTA and VCF input files."
            )

    pa.add_argument(
            '--genome-nmask-filename', type=str, dest='genome_nmask_filename',
            default='',
            help="Filename of a FASTA file to which the N-masked geneome will be stored.with the reference genome. " +
            "Note that chromosome names and genomic positions must match in the genome FASTA and VCF input files."
            )


    # ================================================================================================================

    pa.add_argument(
            "--GTF-file-offset", type=int, dest="gtf_file_offset",
            default=1,
            help="Force a specific offset for the genomic position of the first basepair for the GTF genomic feature file " +
            "specified as the last paramter on the command line. The default is 1 (--GTF_file-offset=1). If in doubt use " +
            "--check-consistency for perfoming a validation of genomic coordinates for your GTF genomic feature file."
            )

    pa.add_argument(
            "--GTF-file-end-not-included", action="store_true", dest="gtf_file_end_not_included",
            help="This is a flag that can be used for a GFF/GTF file with intervals where the end is not included. " +
            "The default is a ( ] interval where the end position is the index after the last position of the interval. " +
            "This corresponds to end_included and can be overridden using --GTF-file-end-not-included. If in doubt use " +
            "--check-consistency for perfoming a validation of genomic intervals for your GTF genomic feature file."
            )

    args = pa.parse_args()

    if args.version:
        print(HTSeq.__version__)
        sys.exit()

    if args.consistency_check:   # perform checks depending on which files are provided: BAM, GTF, VCF, requires genome FASTA
        consistency_check(args.genome_fasta_filename, args.vcf_filename, args.vcf_file_offset, args.featuresfilename, args.samfilenames, test_num=args.test_num)
        sys.exit()

    if args.allele_counts:
        if args.mode != "union":
            print("WARNING: For allelic analysis count mode=union is required!")
            sys.exist()

    if args.list_vcf_file_samples:
        if os.path.exists(args.vcf_filename):
            print_vcf_info(args.vcf_filename)
        else:
            print("VCF file not found. Use --SNP-vcf-filename to specify the path to a VCF file.")
        sys.exit()
        
    if args.allele_counts:                              # <=============================== load SNPs
        if not args.quiet:
            print(args)
        if args.force_generate_bed_file_cache or not os.path.exists(args.bed_filename):
            if os.path.exists(args.vcf_filename):
                if args.use_fast_vcf:
                    SNPs = read_vcf(args.vcf_filename, args.ref_samples, args.alt_samples, omit_ref_genome=args.omit_ref_genome, force_biallelic=args.force_biallelic, vcf_file_offset=args.vcf_file_offset, silent=args.quiet)
                else:
                    SNPs = encode_vcf_file(args.vcf_filename, args.ref_samples, args.alt_samples, omit_ref_genome=args.omit_ref_genome, force_biallelic=args.force_biallelic, vcf_file_offset=args.vcf_file_offset, silent=args.quiet)
                if args.use_bed_file_cache or args.force_generate_bed_file_cache:
                    if os.access(args.bed_filename, os.W_OK) or os.access(os.path.dirname(args.bed_filename), os.W_OK):
                        SNPs.write_bedgraph_file(args.bed_filename, strand=".")
                        if not args.quiet:
                            print("SNP information written to BED cache file", args.bed_filename)
                    else:
                        if not args.quiet:
                            print("BED file path", args.bed_filename, "not writeable. Use --SNP-bed-filename to specify a valid path to write a BED file.")
            else:
                print("VCF file not found. Use --SNP-vcf-filename to specify the path to a VCF file.")
                sys.exit()
        elif args.use_bed_file_cache and os.path.exists(args.bed_filename):
            SNPs = HTSeq.GenomicArray.from_bedgraph_file( args.bed_filename, strand=".", typecode="i" )
            for chrom in SNPs.chrom_vectors:
                SNPs[chrom]["."].iv.end = sys.maxsize   # needed to avoid error as GenomicArray can be queried at any length
        else:
            print("No SNP information found. Use --SNP-vcf-filename or --SNP-bed-filename to specify the path to a VCF file or BED file.")
            sys.exit()
    else:
        SNPs = 0                                        # this is just a dummy value to pass as parameter

    if args.genome_nmask:
        nmask_genome(args.genome_fasta_filename, args.genome_nmask_filename, args.vcf_filename, args.ref_samples, args.alt_samples, omit_ref_genome=args.omit_ref_genome, vcf_file_offset=args.vcf_file_offset, silent=args.quiet)
        sys.exit()
        
    add_genome_annotation = []
    if args.add_chromosome_information:
        add_genome_annotation.append("chrom")
    if args.add_genelength_information:
        add_genome_annotation.append("gene_length")
    if args.add_exonlength_information:
        add_genome_annotation.append("exon_length")

    warnings.showwarning = my_showwarning
    try:
        count_reads_in_features(
            args.samfilenames,
            args.featuresfilename,
            args.order,
            args.max_buffer_size,
            args.stranded,
            args.mode,
            args.nonunique,
            args.secondary_alignments,
            args.supplementary_alignments,
            args.featuretype,
            args.idattr,
            args.additional_attr,
            args.quiet,
            args.minaqual,
            args.samouts,
            args.samout_format,
            args.output_delimiter,
            args.output_filename,
            args.output_append,
            args.nprocesses,
            args.feature_query,
            args.allele_counts,    # <------------------------------------------------- allelic analysis
            SNPs,
            args.snp_min_qual,
            args.add_genome_annotation,
            args.gtf_file_offset,
            args.gtf_file_end_not_included,
            )
    except:
        sys.stderr.write("  %s\n" % str(sys.exc_info()[1]))
        sys.stderr.write("  [Exception type: %s, raised in %s:%d]\n" %
                         (sys.exc_info()[1].__class__.__name__,
                          os.path.basename(traceback.extract_tb(
                              sys.exc_info()[2])[-1][0]),
                          traceback.extract_tb(sys.exc_info()[2])[-1][1]))
        sys.exit(1)


if __name__ == "__main__":
    main()
