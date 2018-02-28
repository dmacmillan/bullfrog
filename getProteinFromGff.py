import argparse
import os
import sys
import logging
import pysam
from itertools import groupby
from customFunctions import *

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Given a gff file, extract protein sequences based on transcript sequences'
    )
    
    parser.add_argument('gff', help='gff file input')
    parser.add_argument('transcripts', help='Fasta file of transcripts from which to fetch')
    parser.add_argument('-ts', '--trim_start', action='store_true', help='If enabled, trim the start of the protein to the nearest methionine amino acid')
    parser.add_argument("-l", "--log", dest="logLevel", default='warning', choices=['debug', 'info', 'warning', 'error', 'critical'], help='set the logging level. default = "warning"')
    parser.add_argument('-n', '--name', default='protein.tsv', help='Name for the file output. Default is "protein.tsv"')
    parser.add_argument('-o', '--outdir', default=os.getcwd(), help='Path to output to. Default is {}'.format(os.getcwd()))
    
    args = parser.parse_args()
    
    # Create output directory if it does not exist
    if not os.path.isdir(args.outdir):
        try:
            os.makedirs(args.outdir)
        except OSError:
            pass
    
    # Logging
    logging.basicConfig(
        filename = os.path.join(args.outdir, '{}.log'.format(os.path.basename(__file__))),
        level = getattr(logging, args.logLevel.upper()),
        filemode = 'w'
    )
    
    # Open gff file
    gff_file_object = open(args.gff, 'r')
    
    # Open transcripts fasta
    transcripts = pysam.FastaFile(args.transcripts)
    
    # Load gff into Pysam object
    gff = pysam.tabix_iterator(
        gff_file_object,
        parser = pysam.asGTF()
    )
    
    # Sorting GFF is annoying in bash alone,
    # so I opt to store it in memory and then spit it back
    # out per transcript
    txts = {}

    with open(os.path.join(args.outdir, args.name), 'w') as o:
        for line in gff:
            attributes = dict([x.split('=') for x in line.attributes.split('; ')])
            name = attributes['Name']
            if name not in txts:
                txts[name] = [line] 
            else:
                txts[name].append(line)
        for txt in txts:
            logging.debug('txt: {}'.format(txt))
            txts[txt] = sorted(txts[txt], key = lambda x: x.start)
            for x in txts[txt]:
                logging.debug('x: {}'.format(x))
            strand = txts[txt][0].strand
            cds = [x for x in txts[txt] if x.feature == 'CDS']
            if not cds:
                logging.warning('Transcript "{}"" has no CDS!'.format(txt))
                continue
            seq = ('').join([transcripts.fetch(x.contig, x.start, x.end) for x in cds])
            if strand == '-':
                seq = rev_comp(seq)
            aa = trans(seq)
            # Strip '-' suffix
            aa = aa.strip('-')
            # Trim sequence to nearest methionine
            if args.trim_start:
                try:
                    aa = aa[aa.index('M'):]
                except ValueError:
                    logging.warning(
                        'Sequence "{}" has no methionine. {}'.format(
                            aa,
                            txt
                        )
                    )
            logging.debug('aa: {}'.format(aa))
            o.write(
                ('\t').join((
                    cds[0].contig,
                    txt,
                    aa
                )) + '\n'
            )

#    with open(os.path.join(args.outdir, args.name), 'w') as o:
#        for key, group in groupby(
#            gff,
#            #key = lambda entry: dict([x.split('=') for x in entry.attributes.split(';')])['Name']
#            key = lambda entry: dict([x.split('=') for x in entry.attributes.split('; ')])['Name']
#        ):
#            logging.debug('Group key: "{}"'.format(key))
#            temp = [x for x in group]
#            for x in temp:
#                logging.debug('\tx: {}'.format(x))
#            cds = [x for x in temp if x.feature == 'CDS']
#            #cds = [x for x in group if x.feature == 'CDS']
#            logging.debug('cds: {}'.format(cds))
#            strand = cds[0].strand
#            logging.debug('strand: {}'.format(strand))
#            if strand == '+':
#                seq = ('').join([transcripts.fetch(x.contig, x.start, x.end) for x in cds])
#            else:
#                if len(cds) > 1:
#                    seq = ('').join([transcripts.fetch(x.contig, x.start, x.end)[::-1] for x in cds])
#                    seq = comp(seq)
#                else:
#                    seq = ('').join([transcripts.fetch(x.contig, x.start, x.end) for x in cds])
#                    seq = rev_comp(seq)
#            aa = trans(seq)
#            # Strip '-' suffix
#            aa = aa.strip('-')
#            logging.debug('key: {}'.format(key))
#            logging.debug('seq: {}'.format(seq))
#            # Trim sequence to nearest methionine
#            if args.trim_start:
#                try:
#                    aa = aa[aa.index('M'):]
#                except ValueError:
#                    logging.warning(
#                        'Sequence "{}" has no methionine. {}'.format(
#                            aa,
#                            key
#                        )
#                    )
#            logging.debug('aa: {}'.format(aa))
#            o.write(
#                ('\t').join((
#                    cds[0].contig,
#                    key,
#                    aa
#                )) + '\n'
#            )
