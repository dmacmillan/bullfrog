import argparse
import os
import sys
import logging
from itertools import groupby
import time
from customFunctions import *

# kwargs recognized:
    # gffs - A dictionary, where the key 
        # is the name and where the value
        # is a dictionary of mRNA names each with a 
        # list of (start, stop) values
    # exons - A dictionary of gene name -> scaffold name,
        # mRNA name -> start, end, and strand
def parse_psl(psl, feature, **kwargs):
    stime = time.time()
    # Keep track of Qnames so as to avoid looking at
    # secondary alignments
    to_analyze = {}
    with open(psl, 'r') as f:
        logging.debug('Reading psl: {}'.format(psl))
        # Read 5 lines to skip the header
        for i in range(5):
            f.readline()
        for line in f:
            cols = line.strip().split('\t')
            match = int(cols[0])
            cols[0] = match
            # This mRNA name is derived from
            # the maker 
            mrna_name = cols[9]
            scaffold_name = cols[13]
            blocksizes = [int(x) for x in cols[18].strip(',').split(',')]
            nblocks = len(blocksizes)
            # This stuff just ensures that only the alignment
            # with the most matches is looked at
            if mrna_name in to_analyze:
                if match > to_analyze[mrna_name]['cols'][0]:
                    to_analyze[mrna_name] = {
                        'cols': cols,
                        'scaffold_name': scaffold_name
                    }
            else:
                to_analyze[mrna_name] = {
                    'cols': cols,
                    'scaffold_name': scaffold_name
                }
        for mrna_name in to_analyze:
            cols = to_analyze[mrna_name]['cols']
            scaffold_name = to_analyze[mrna_name]['scaffold_name']
            # psl uses a 0-based coordinate system
            # while gff uses a 1-based
            tstarts = [int(x)+1 for x in cols[20].strip(',').split(',')]
            blocksizes = [int(x) for x in cols[18].strip(',').split(',')]
            logging.debug('mrna_name: {}'.format(mrna_name))
            logging.debug('cols: {}'.format(cols))
            logging.debug('tstarts: {}'.format(tstarts))
            logging.debug('blocksizes: {}'.format(blocksizes))
            strand = cols[8]
            # Why is this condition here? Removing...
            # if exons or exons == {}:
                # Get the gene name
            gene_name = mrna_name.rsplit('-',2)[0]
            if 'exons' not in kwargs:
                kwargs['exons'] = {
                    gene_name: {
                        'scaffold_name': scaffold_name,
                        'mrna_name': {
                            mrna_name: {
                                'start': tstarts[0],
                                'end': tstarts[-1] + blocksizes[-1],
                                'strand': strand
                            }
                        }
                    }
                }
            elif gene_name not in kwargs['exons']:
                kwargs['exons'][gene_name] = {
                    'scaffold_name': scaffold_name,
                    'mrna_name': {
                        mrna_name: {
                            'start': tstarts[0],
                            'end': tstarts[-1] + blocksizes[-1],
                            'strand': strand
                        }
                    }
                }
            elif mrna_name not in kwargs['exons'][gene_name]['mrna_name']:
                kwargs['exons'][gene_name]['mrna_name'][mrna_name] = {
                    'start': tstarts[0],
                    'end': tstarts[-1] + blocksizes[-1],
                    'strand': strand
                }
            else:
                kwargs['exons'][gene_name]['mrna_name'][mrna_name]['start'] = min(
                    kwargs['exons'][gene_name]['mrna_name'][mrna_name]['start'],
                    tstarts[0]
                )
                kwargs['exons'][gene_name]['mrna_name'][mrna_name]['end'] = max(
                    kwargs['exons'][gene_name]['mrna_name'][mrna_name]['end'],
                    tstarts[-1] + blocksizes[-1]
                )
            compare_blocks = True
            # Check number of blocks is equal
            try:
                maker_blocksizes = set([x[1]-x[0] for x in kwargs['gffs']['maker'][mrna_name][feature]])
                maker_num_blocks = len(kwargs['gffs']['maker'][mrna_name][feature])
            except KeyError:
                logging.warning(
                    '"{}" not found in reference gff "{}"! ' \
                    'Cannot check block sizes!'.format(mrna_name, 'maker')
                )
                compare_blocks = False
            try:
                gmap_blocksizes = set([x[1]-x[0] for x in kwargs['gffs']['gmap'][mrna_name][feature]])
                gmap_num_blocks = len(kwargs['gffs']['gmap'][mrna_name][feature])
            except KeyError:
                logging.warning(
                    '"{}" not found in reference gff "{}"! ' \
                    'Cannot check block sizes!'.format(mrna_name, 'gmap')
                )
                gmap_blocksizes = set()
                gmap_num_blocks = 0
            blat_blocksizes = set([x-1 for x in blocksizes])
            blat_num_blocks = len(blocksizes)
            if compare_blocks:
                if ((gmap_num_blocks == maker_num_blocks) and (gmap_blocksizes == maker_blocksizes)):
                    pass
                elif ((blat_num_blocks == maker_num_blocks) and (blat_blocksizes == maker_blocksizes)):
                    pass
                else:
                    in_maker_not_gmap = maker_blocksizes - gmap_blocksizes
                    in_maker_not_blat = maker_blocksizes - blat_blocksizes
                    in_gmap_not_maker = gmap_blocksizes - maker_blocksizes
                    in_gmap_not_blat = gmap_blocksizes - blat_blocksizes
                    in_blat_not_maker = blat_blocksizes - maker_blocksizes
                    in_blat_not_gmap = blat_blocksizes - gmap_blocksizes
                    output = [
                        scaffold_name,
                        mrna_name,
                        str(maker_num_blocks),
                        str(gmap_num_blocks),
                        str(blat_num_blocks),
                        (',').join([str(x) for x in in_maker_not_gmap]) or 'NA',
                        (',').join([str(x) for x in in_maker_not_blat]) or 'NA',
                        (',').join([str(x) for x in in_gmap_not_maker]) or 'NA',
                        (',').join([str(x) for x in in_gmap_not_blat]) or 'NA',
                        (',').join([str(x) for x in in_blat_not_maker]) or 'NA',
                        (',').join([str(x) for x in in_blat_not_gmap]) or 'NA',
                    ]
                    logging.warning(
                        ('\t').join(output)
                    )
            for i,start in enumerate(tstarts):
                out = Gff()
                out.source = 'pslToGff'
                out.feature = feature
                out.start = tstarts[i]
                logging.debug('blocksizes: {}'.format(blocksizes))
                logging.debug('tstarts: {}'.format(tstarts))
                logging.debug('i: {}'.format(i))
                out.end = tstarts[i] + blocksizes[i]
                out.score = '0'
                # These boundaries can be fixed in this case
                out.seqname = scaffold_name
                out.strand = strand
                out.frame = '0'
                attributes = {
                    'ID': '{}:{}:{}'.format(mrna_name, feature, i),
                    'Parent': mrna_name,
                    'Name': gene_name
                }
                out.attribute = attributes
                print(out)
    logging.info('Parsed "{}" in {}s'.format(psl, time.time()-stime))
    return exons

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Given psl files, return a gff file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('maker_gff', help='reference gff file with the "ID" attribute set to the Qname of the psl')
    parser.add_argument('gmap_gff', help='Another gff file')
    parser.add_argument('maker_proteins', help='Maker protein sequences in tsv format')
    parser.add_argument('gmap_proteins', help='Gmap protein sequences in tsv format')
    parser.add_argument('blat_proteins', help='Blat protein sequences in tsv format')
    parser.add_argument('-ep', '--exon_psls', nargs='+', help='exon psl file input')
    parser.add_argument('-cp', '--cds_psls', nargs='+', help='CDS psl file input')
    #parser.add_argument('-f', '--feature', default='exon', help='A feature name for the gff output, must match the feature names from the gff input')
    parser.add_argument("-l", "--log", dest="logLevel", default='warning', choices=['debug', 'info', 'warning', 'error', 'critical'], help='set the logging level. default = "warning"')
    
    args = parser.parse_args()
    
    # Logging
    logging.basicConfig(
        filename = os.path.join(os.getcwd(), '{}.log'.format(os.path.basename(__file__))),
        level = getattr(logging, args.logLevel.upper()),
        filemode = 'w'
    )

    # Log columns
    # logging.warning(
    #     'Columns = v3(BLAT PSL) scaffold, v2(Maker GFF) transcript, ' \
    #     'num v3 blocks, num v2 blocks, ' \
    #     'blocksizes in v3 not v2, blocksizes in v2 not v3, ' \
    #     'block size differences v3-v2 for each'
    # )

    # Keep track of gene sizes
    genes = {}

    # Load gff into memory
    maker_gff = {}
    gmap_gff = {}

    for gff in Gff.read(args.maker_gff):
        if gff.feature not in ('CDS', 'exon'):
            continue
        names = gff.attribute['Parent'].split(',')
        for name in names:
            if name not in maker_gff:
                maker_gff[name] = {}
            if gff.feature not in maker_gff[name]:
                maker_gff[name][gff.feature] = []
            maker_gff[name][gff.feature].append((gff.start, gff.end))

    for gff in Gff.read(args.gmap_gff):
        if gff.feature not in ('CDS', 'exon'):
            continue
        name = gff.attribute['Name']
        if name not in gmap_gff:
            gmap_gff[name] = {}
        if gff.feature not in gmap_gff[name]:
            gmap_gff[name][gff.feature] = []
        gmap_gff[name][gff.feature].append((gff.start, gff.end))

    # Keep track of exon boundaries
    exons = {}
    if args.exon_psls:
        for psl in args.exon_psls:
            exons = parse_psl(psl,
                'exon',
                exons = exons,
                gffs = {
                    'maker': maker_gff,
                    'gmap': gmap_gff
                }
            )
    if args.cds_psls:
        for psl in args.cds_psls:
            parse_psl(psl,
                'CDS',
                gffs = {
                    'maker': maker_gff,
                    'gmap': gmap_gff
                }
            )

    for gene in exons:
        out = Gff()
        out.source = 'pslToGff'
        out.feature = 'gene'
        out.start = min([exons[gene]['mrna_name'][x]['start'] for x in exons[gene]['mrna_name']])
        out.end = max([exons[gene]['mrna_name'][x]['end'] for x in exons[gene]['mrna_name']])
        out.score = '0'
        # These boundaries can be fixed in this case
        out.seqname = exons[gene]['scaffold_name']
        out.strand = exons[gene]['mrna_name'][list(exons[gene]['mrna_name'])[0]]['strand']
        out.frame = '0'
        # Added Name attribute to allow for grouping
        attributes = {
            'ID': gene,
            'Name': gene
        }
        out.attribute = attributes
        print(out)
        for mrna in exons[gene]['mrna_name']:
            out.feature = 'mRNA'
            out.start = exons[gene]['mrna_name'][mrna]['start']
            out.end = exons[gene]['mrna_name'][mrna]['end']
            # These boundaries can be fixed in this case
            attributes = {
                'ID': mrna,
                'Parent': gene,
                'Name': gene
            }
            out.attribute = attributes
            print(out)