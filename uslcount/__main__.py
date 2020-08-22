import sys
import argparse
from .bamutils import is_paired, EmptyBamError
from .main import count_bam, analyze_bam, write_count_res, \
    write_all_counts, write_analyze_res
from .gdata import GenomicData

tasks = ['build', 'analyze', 'count']
if len(sys.argv)<2 or (sys.argv[1] not in tasks):
    usage = """Usage: python3 -m uslcount <task> [options]

    available tasks are
        {:s}
        {:s}
        {:s}
    For task specific options:

    > python3 -m uslcount <task> -h (--help)
    """.format(*tasks)
    print(usage)
    sys.exit()

elif sys.argv[1] == 'build':
    
    prsr=argparse.ArgumentParser(prog='build',
        description='Build genomic data DB for analysis')
    prsr.add_argument('gtf', metavar='GTF', type=str,
                      help='Path to file GTF annotation file')
    prsr.add_argument('out_dir', metavar='DIR', type=str, help='Output directory')
    args = prsr.parse_args(sys.argv[2:])
    
    g = GenomicData.build(args.gtf, args.out_dir, sj_overhang=3)

elif sys.argv[1] in ('analyze', 'count'):

    # First add arguments common for both tasks
    prsr=argparse.ArgumentParser()
    prsr.add_argument('-out',  metavar='OUT', help='path to outfile',
                      type=str, default='stdout')
    prsr.add_argument('-count_intron_overlapping', choices=['yes', 'no'], metavar='yes|no',
                      help='count alignments overlapping gene introns ("no" by default)',
                      default='no')
    prsr.add_argument('-count_secondary', choices=['yes', 'no'], metavar='yes|no',
                      help='count secondary alignments ("no" by default)',
                      default='no')
    prsr.add_argument('-count_supplementary', choices=['yes', 'no'], metavar='yes|no',
                      help='count supplementary alignments ("no" by default)',
                      default='no')
    prsr.add_argument('-count_qc_failed', choices=['yes', 'no'], metavar='yes|no',
                      help='count qc failed (0x200 bit flag set) alignments ("yes" by default)',
                      default='yes')
    prsr.add_argument('-mapq',  metavar='INT',
                      help='minimal MAPQ considered (default 10)',
                      type=int, default=10)
    prsr.add_argument('-quiet',  help='do not print anything',
                      action='store_true', default=False)
    prsr.add_argument('gdir', metavar='DATA_DIR', type=str,
                      help='Path to genomic data directory')
    prsr.add_argument('bam',  metavar='BAM', help='Path to bam/sam file to analyze',
                      type=str)

    # Add task specific arguments and description
    if sys.argv[1] == 'analyze':
        prsr.prog = 'analyze'
        prsr.description='Gather counts from BAM/SAM alignment file prepared '+\
                          'with unstranded library and Confidence Score'
    elif sys.argv[1] == 'count':
        prsr.add_argument('-strand',  metavar='{F, R, N}',
                          help='Library strand. N - unstranded library (default), F - forward, R - reverse.',
                          choices=['F', 'R', 'N'], type=str, default='N')
        prsr.add_argument('-dev', help='Switch to full count output mode for development purposes.',
                          action='store_true', default=False)
        prsr.prog = 'count'
        prsr.description='Gather counts from BAM/SAM alignment file.'

    # Parse arguments
    args = prsr.parse_args(sys.argv[2:])
    ignore_intronic = False if args.count_intron_overlapping=='yes'else True
    ignore_secondary = False if args.count_secondary=='yes' else True
    ignore_supplementary = False if args.count_supplementary=='yes' else True
    ignore_qc_failed = False if args.count_intron_overlapping=='yes' else True

    # Load Genome
    gdat = GenomicData.load(args.gdir, quiet=args.quiet)

    # Run analysis
    if sys.argv[1] == 'analyze':
        try:
            res = analyze_bam(args.bam, gdat,
                              ignore_intronic=ignore_intronic,
                              ignore_secondary=ignore_secondary,
                              ignore_supplementary=ignore_supplementary,
                              ignore_qc_failed=ignore_qc_failed,
                              mapq=args.mapq, 
                              quiet=args.quiet)
        except EmptyBamError:
            print('No alignments')
            sys.exit()
        
    elif sys.argv[1] == 'count':
        if args.strand == 'N':
            stranded = False
        else:
            stranded = True
        paired = is_paired(args.bam)
        if args.dev:
            mode='dev'
        else:
            mode='count'
        res = count_bam(args.bam, gdat, mode=mode, paired=paired,
                        strand=args.strand,
                        ignore_intronic=ignore_intronic,
                        ignore_secondary=ignore_secondary,
                        ignore_supplementary=ignore_supplementary,
                        ignore_qc_failed=ignore_qc_failed,
                        mapq=args.mapq, 
                        quiet=args.quiet)

    # Write results
    if args.out == 'stdout':
        out_fh = sys.stdout
        opened_file = False
    else:
        out_fh = open(args.out, 'wt')
        opened_file = True

    if sys.argv[1] == 'analyze':
        write_analyze_res(res, out_fh)
    elif sys.argv[1] == 'count':
        if mode=='count':
            write_count_res(res, out_fh, stranded=stranded)
        elif mode=='dev':
            write_all_counts(res, out_fh)

    if opened_file:
        out_fh.close()
