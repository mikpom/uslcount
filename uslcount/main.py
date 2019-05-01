from __future__ import print_function, division, unicode_literals
import warnings
import math
import pysam
import numpy as np
from .bamutils import read_pairs, lib_param, \
    aln_ok, get_skipped, bam_lsizes
from .genomic import steps
from .gdata import GenomicData
from .trees import pred_feats, dtrees, minprob, maxprob, predict_proba

# Defining fields for data structures
cnt_cols = ['ustno', 'ust', 'st', 'stno', 'intr', 'jctn']
cnt_fields = [('gid', 'U40')]
for col in cnt_cols:
    cnt_fields.append((col, np.uint32))

stats_fields = [('unq_frac', np.float64),
                ('jctF', np.float64), ('nearbyRPK', np.float64),
                ('RPK', np.float64), ('exo2intr', np.float64),
                ('prob', np.float64)]

# Collect per-gene counts from BAM file
def count_bam(bam_file, gdat, paired=False, strand='N', ignore_intronic=True,
              ignore_secondary=True, ignore_supplementary=True,
              ignore_qc_failed=True, mapq=10, quiet=False, mode='count'):

    def _print_progress(cnt, paired=False, quiet=False):
        if (cnt % 1000000 == 0) and (not quiet):
            if paired:
                print(str(int(cnt/1000000))+' mln alignment pairs processed')
            else:
                print(str(int(cnt/1000000))+' mln alignments processed')

    # mode is one of count, analyze, dev
    modes = ['count', 'analyze', 'dev']
    if mode not in modes:
        raise ValueError('mode should one of '+','.join(modes))
    if mode=='dev' and strand=='N':
        raise ValueError('Strand should be indicated for dev mode')

    # Get proper alignment iterator
    if paired:
        aln_iter = read_pairs(bam_file, quiet=quiet)
    else:
        aln_iter = pysam.AlignmentFile(bam_file, 'r')

    seen_chroms = set()

    # Initialize count dictionary
    knt = {}
    for gn in gdat.genes:
        knt[gn] = {c:0 for c in cnt_cols}

    # Iterate over alignments
    if not quiet:
        print('Starting processing alignments')
    aln_n = 0
    for aln in aln_iter:
        aln_n += 1
        _print_progress(aln_n, paired=paired, quiet=quiet)

        # collect alignment blocks and junctions
        jctns = [] # list of splice junction as (chrom, start, end) tuples
        blocks = [] # list of alignemnts blocks
        if paired:
            chroms = set()
            blocks = []
            for seg in aln:
                if not seg is None:
                    if aln_ok(seg, ignore_secondary=ignore_secondary,
                              ignore_supplementary=ignore_supplementary,
                              ignore_qc_failed=ignore_qc_failed,
                              mapq=mapq):
                        blocks += seg.get_blocks()
                        chroms.add(seg.reference_name)
                        if ('N' in seg.cigarstring and mode=='analyze') or (mode=='dev'):
                            jctns += get_skipped(seg)
            if len(chroms) > 1: # mates mapped to different chromosomes
                continue
            elif len(chroms) == 1:
                chrom = list(chroms)[0]
        else:
            if aln_ok(aln, ignore_secondary=ignore_secondary,
                      ignore_supplementary=ignore_supplementary,
                      ignore_qc_failed=ignore_qc_failed,
                      mapq=mapq):
                blocks = aln.get_blocks()
                chrom = aln.reference_name
                if ('N' in aln.cigarstring and mode=='analyze') or (mode=='dev'):
                    jctns = get_skipped(aln)

        if len(blocks)==0:
            continue
        if not chrom in gdat.chroms:
            if not chrom in seen_chroms:
                msg = 'Reference sequence "{:s}" encountered in BAM alignment'+\
                      ' record but not found in genomic data.'
                warnings.warn(msg.format(chrom))
                seen_chroms.add(chrom)
            continue

        # Collect intronic features
        fs_intr = set()
        for blk in blocks:
            for iv, fs in steps(gdat.introns[chrom], blk[0], blk[1]).items():
                fs_intr = fs_intr.union(fs)

        # Infer alignment strand
        if strand != 'N':
            _aln = aln[0] if paired else aln
            if strand == "F":
                aln_strand = '-' if _aln.is_reverse else '+'
            elif strand == 'R':
                aln_strand = '+' if _aln.is_reverse else '-'

        # Assign unstranded and stranded-no-overlap counts
        # We do so unless library is stranded and mode is 'count'
        if not (mode=='count' and strand!='N'):
            fs_ust = set()
            for blk in blocks:
                for iv, fs in steps(gdat.exons_ust[chrom], blk[0], blk[1]).items():
                    fs_ust = fs_ust.union(fs)

            if len(fs_ust) > 1:
                for gn in fs_ust:
                    if (not gn in fs_intr) or (not ignore_intronic):
                        knt[gn]['ust'] += 1
            elif len(fs_ust) == 1:
                gn = list(fs_ust)[0]
                if (not gn in fs_intr) or (not ignore_intronic):
                    knt[gn]['ustno'] += 1
                    knt[gn]['ust'] += 1
                if strand != 'N':
                    if (gdat.strands[gn] == aln_strand) and \
                       ((not gn in fs_intr) or (not ignore_intronic)):
                        knt[gn]['stno'] += 1

        # Assign intronic counts
        if mode in ('analyze', 'dev'):
            if len(fs_intr) >= 1:
                for fs in fs_intr:
                    knt[fs]['intr'] += 1

        # Assign junction counts
        # strand is completely ignored here
        if mode in ('analyze', 'dev'):
            counted = set()
            if len(jctns) >= 1:
                for jctn in jctns:
                    try:
                        # sometimes genes can share a splice junction
                        # #sadButTrue
                        gg = gdat.intr2gn[(chrom, jctn[0], jctn[1])]
                    except KeyError:
                        continue
                    for gn in gg:
                        if not gn in counted:
                            knt[gn]['jctn']+= 1
                            counted.add(gn)

        # collect stranded counts
        # we do so if mode is 'dev' or mode is 'count' and library is stranded
        if (mode=='count' and strand != 'N') or (mode=='dev'):
            fs_st = set()
            for blk in blocks:
                for iv, fs in steps(gdat.exons_st[chrom][aln_strand], blk[0], blk[1]).items():
                    fs_st = fs_st.union(fs)
            if len(fs_st)==1:
                gn = list(fs_st)[0]
                if (not gn in fs_intr) or (not ignore_intronic):
                    knt[gn]['st'] += 1

    # Prepare numpy structured array for output
    arr = np.zeros(shape=(len(gdat.genes), ), dtype=cnt_fields)
    gg = list(gdat.genes)
    arr['gid'] = gg
    for col in cnt_cols:
        arr[col] = [knt[g][col] for g in gg]
    return arr

# Collect alignment statistics for analyze task
def collect_aln_stats(cnt, gdat, jct_freq):
    if set(cnt['gid']) != set(gdat.genes):
        raise ValueError('Counts should be provided for exactly '+\
                         'the genes in GenomicData object')
    dtype = [(name, cnt.dtype.fields[name][0]) for name in cnt.dtype.names]
    stats = np.zeros(shape=cnt.shape, dtype=dtype+stats_fields)
    for fld in cnt.dtype.names:
        stats[fld] = cnt[fld]

    # Reads per kilobase of gene span
    rpk = {}
    for rec in stats:
        gid = rec['gid']
        if gdat.len_ust_unq[gid] > 0:
           _rpk  = rec['ustno'] * 1000 / gdat.len_ust_unq[gid]
        else:
            _rpk = 0
        rec['RPK'] = _rpk
        rpk[gid] = _rpk

    # Max RPK of nearby genes
    dlt = 2000
    for rec in stats:
        gid = rec['gid']
        gn = gdat.genes[gid]
        gn_in_rgn = gdat.genes_rgn(gn.chrom, max(0, gn.start-dlt), gn.end+dlt)
        gn_in_rgn.remove(gid)
        if len(gn_in_rgn) > 0:
            rec['nearbyRPK'] = max([rpk[g] for g in gn_in_rgn])
        else:
            rec['nearbyRPK'] = 0

    # Exonic to intronic ratio
    for rec in stats:
        gid = rec['gid']
        if (gdat.len_ust_intr[gid] == 0) or (gdat.len_ust_unq[gid] == 0):
            rec['exo2intr'] = np.nan
        else:
            intr_x = rec['intr'] / gdat.len_ust_intr[gid]
            exo_x = rec['ust'] / gdat.len_ust_unq[gid]
            if (intr_x != 0) and (exo_x !=0):
                rec['exo2intr'] = np.log2(exo_x / intr_x)
            elif (intr_x==0) and (exo_x==0):
                rec['exo2intr'] = np.nan
            elif (intr_x!=0) and (exo_x==0):
                rec['exo2intr'] = -np.inf
            elif (intr_x==0) and (exo_x!=0):
                rec['exo2intr'] = np.inf
            else:
                print('Should not be here')

    # Junction factor
    for rec in stats:
        gid = rec['gid']
        expect_jct = jct_freq[gid] * rec['ust']
        obs_jct = rec['jctn']
        if (expect_jct != 0) and (obs_jct !=0):
            rec['jctF'] = np.log2(obs_jct / expect_jct)
        elif (expect_jct==0) and (obs_jct==0):
            rec['jctF'] = np.nan
        elif (expect_jct!=0) and (obs_jct==0):
            rec['jctF'] = -np.inf
        elif (expect_jct==0) and (obs_jct!=0):
            rec['jctF'] = np.inf
        else:
            print('Should not be here')

    # Fraction of unique unstranded gene span
    stats['unq_frac'] = [gdat.len_ust_unq[g] / gdat.len_ust_tot[g] for g in stats['gid']]

    # Calculate decision probs
    maxfloat = np.finfo(np.float64).max
    minfloat = np.finfo(np.float64).min
    for rec in stats:
        gid = rec['gid']
        gf = gdat.genes[gid]
        if rec['ustno']>0:
            dt = np.asarray(rec[pred_feats[gf.struct]].item())
            np.clip(dt, minfloat, maxfloat, out=dt)
            pr = predict_proba(dt, dtrees[gf.struct])
            rec['prob'] = pr['err']
        else:
            rec['prob'] = np.nan

    return stats

# Run analysis of BAM from unstranded data
def analyze_bam(bam, gdat, strand='N', ignore_intronic=True,
                ignore_secondary=True, ignore_supplementary=True,
                ignore_qc_failed=True, mapq=10, 
                quiet=False):

    rl, paired = lib_param(bam)

    cnt = count_bam(bam, gdat, paired=paired, mode='analyze',
                    strand=strand,
                    ignore_intronic=ignore_intronic,
                    ignore_secondary=ignore_secondary,
                    ignore_supplementary=ignore_supplementary,
                    ignore_qc_failed=ignore_qc_failed,
                    mapq=mapq, quiet=quiet)

    if not quiet:
        print('Parsing BAM is done!')
        print('Calculating junction probabilities...')
    if paired:
        lsizes = bam_lsizes(bam)
        lsizes = lsizes[lsizes>=-rl]
        freqs = gdat.jct_freq(rl, lsizes, quiet=quiet)
        stats = collect_aln_stats(cnt, gdat, freqs)
    else:
        freqs = gdat.jct_freq(rl)
        stats = collect_aln_stats(cnt, gdat, freqs)
    return stats

def write_count_res(cnt, fh, stranded=False):
    cnt_col = 'st' if stranded else 'ustno'
    fh.write('gene_id\tcount\n')
    for rec in cnt:
        fh.write(rec['gid']+'\t'+str(rec[cnt_col])+'\n')

def write_all_counts(cnt, fh):
    fh.write('gene_id\t')
    fh.write('\t'.join(cnt_cols))
    fh.write('\n')
    for rec in cnt:
        fields = [rec['gid']]
        fields += [str(rec[c]) for c in cnt_cols]
        fh.write('\t'.join(fields))
        fh.write('\n')

field2col = {'ustno':'count', 'jctn':'junctions'}
def write_analyze_res(stats, fh):
    cols = ['ustno', 'jctn']
    fh.write('gene_id\t')
    for col in cols:
        fh.write(field2col.get(col, col)+'\t')
    fh.write('confidence_score\n')

    for rec in stats:
        fields = []
        fields.append(rec['gid'])
        for col in cols:
            fields.append(str(rec[col]))
            
        # probability that gene is erroneously quantitated
        pr = rec['prob']
        
        if not np.isnan(pr):
            score = math.ceil(1000 * (1-pr) + 0.5) / 1000
            fields.append(str(score))
        else:
            fields.append('NaN')
        fh.write('\t'.join(fields)+'\n')
