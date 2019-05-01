import warnings
import numpy as np
import pysam

def aln_ok(aln, ignore_secondary=True, ignore_supplementary=True,
           ignore_qc_failed=True, mapq=1):
    if aln.is_unmapped:
        return False
    if ignore_secondary and aln.is_secondary:
        return False
    if ignore_supplementary and aln.is_supplementary:
        return False
    if ignore_qc_failed and aln.is_qcfail:
        return False
    if aln.mapq < mapq:
        return False
    return True

def is_paired(bamf):
    aln = next(pysam.AlignmentFile(bamf, 'r'))
    return aln.is_paired

# Get gaps in alignment
# These usually correspond to introns
def get_skipped(aln):
    skipped = []
    prev = None
    for blk in aln.get_blocks():
        if (not prev is None) and (prev != blk[0]):
            skipped.append((prev, blk[0]))
        prev = blk[1]
    return skipped

# Get linker size for a mate pair of reads
# Linker is insert size subtracting leangth of reads
# It can be negative
def get_lsize(aln1, aln2):
    a1, a2 = sorted([aln1, aln2], key=lambda a: a.reference_start)
    a1_end = a1.reference_end - a1.get_cigar_stats()[0][3]
    a2_start = a2.reference_start
    return a2_start - a1_end

def proper_pair(a1, a2):
    if a1 and a2 \
       and not a1.is_unmapped and not a2.is_unmapped \
       and a1.reference_name == a2.reference_name:
           return True
    else:
        return False

# Get sample of linker size distribution
# TODO: peek into bam more wisely
def bam_lsizes(bamf, N=100000):
    sz = np.zeros(N, dtype=np.int)
    cnt = 0
    for a1, a2 in read_pairs(bamf, quiet=True):
        if proper_pair(a1, a2):
            sz[cnt] = get_lsize(a1, a2)
            cnt += 1
            if cnt >= N:
                break
    if cnt<N:
        sz = sz[:cnt]

    return sz

def lib_param(bamf):
    bam = pysam.AlignmentFile(bamf, 'r')
    _rl = []
    paired = set()

    cnt = 0
    for aln in bam:
        cnt += 1
        _rl.append(aln.infer_read_length())
        paired.add(aln.is_paired)
        if cnt > 10000:
            break

    if len(paired)>1:
        raise ValueError('Mix of paired and unpaired alignments detected')
    elif len(paired)==1:
        paired = list(paired)[0]
    else:
        print('Should not be here')

    if paired:
        if (not 'SO' in bam.header['HD']) or \
           (not bam.header['HD']['SO'] in ('coordinate', 'queryname')):
            raise ValueError('BAM file is not sorted. '\
                             +'To analyze paired-end data BAM should sorted.')

    rl = np.int(np.median(_rl))
    return (rl, paired)

def aln_mate_hash(aln):
    return (aln.query_name,
            "second" if aln.is_read1 else "first",
            aln.next_reference_name if not aln.mate_is_unmapped else None,
            aln.next_reference_start if not aln.mate_is_unmapped else None,
            aln.reference_name if not aln.is_unmapped else None,
            aln.reference_start if not aln.is_unmapped else None,
            -aln.template_length if not aln.is_unmapped and not aln.mate_is_unmapped else None)

def aln_hash(aln):
    return (aln.query_name,
            "first" if aln.is_read1 else "second",
            aln.reference_name if not aln.is_unmapped else None,
            aln.reference_start if not aln.is_unmapped else None,
            aln.next_reference_name if not aln.mate_is_unmapped else None,
            aln.next_reference_start if not aln.mate_is_unmapped else None,
            aln.template_length if not aln.is_unmapped and not aln.mate_is_unmapped else None)

# stolen from HTSeq library
def read_pairs(fl, max_buffer_size=1000000, quiet=False):
    bam = pysam.AlignmentFile(fl, 'r')
    aln_buffer = {}
    ambiguous_pairing_counter = 0
    for aln in bam:
        matekey = aln_mate_hash(aln)
        if matekey in aln_buffer:
            if len(aln_buffer[matekey]) == 1:
                mate = aln_buffer[matekey][0]
                del aln_buffer[matekey]
            else:
                mate = aln_buffer[matekey].pop(0)
                if ambiguous_pairing_counter == 0:
                    ambiguous_pairing_first_occurance = matekey
                ambiguous_pairing_counter += 1
            if aln.is_read1:
                yield (aln, mate)
            else:
                yield (mate, aln)
        else:
            alnkey = aln_hash(aln)
            if alnkey not in aln_buffer:
                aln_buffer[alnkey] = [aln]
            else:
                aln_buffer[alnkey].append(aln)
            if len(aln_buffer) > max_buffer_size:
                raise ValueError(
                    "Maximum alignment buffer size exceeded while pairing SAM alignments.")

    if len(aln_buffer) > 0:
        if not quiet:
            warnings.warn("Mate records missing for %d records; first such record: %s." %
                          (len(aln_buffer), str(list(aln_buffer.values())[0][0])))
        for aln_list in list(aln_buffer.values()):
            for aln in aln_list:
                if aln.is_read1:
                    yield (aln, None)
                else:
                    yield (None, aln)

    if ambiguous_pairing_counter > 0 and not quiet:
        warnings.warn("Mate pairing was ambiguous for %d records; mate key for first such record: %s." %
                      (ambiguous_pairing_counter, str(ambiguous_pairing_first_occurance)))
