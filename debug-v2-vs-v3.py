if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Given psl files, return a gff file',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('gff', help='reference gff file with the "ID" attribute set to the Qname of the psl')
    parser.add_argument('-ep', '--exon_psls', nargs='+', help='exon psl file input')
    parser.add_argument('-cp', '--cds_psls', nargs='+', help='CDS psl file input')
    parser.add_argument("-l", "--log", dest="logLevel", default='warning', choices=['debug', 'info', 'warning', 'error', 'critical'], help='set the logging level. default = "warning"')
    
    args = parser.parse_args()
    
    # Logging
    logging.basicConfig(
        filename = os.path.join(os.getcwd(), '{}.log'.format(os.path.basename(__file__))),
        level = getattr(logging, args.logLevel.upper()),
        filemode = 'w'
    )

    # Keep track of gene sizes
    genes = {}

    # GFF format columns
    gff_cols = (
        'seqname',
        'source',
        'feature',
        'start',
        'end',
        'score',
        'strand',
        'frame',
        'attribute'
    )

    # Load gff into memory
    gff = {}

    with open(args.gff, 'r') as f:
        for line in f:
            if line[0] == '#':
                continue
            cols = line.strip().split('\t')
            attributes = cols[8]
            start = int(cols[3])
            end = int(cols[4])
            feature = cols[2]
            # Since this script is very specific I can 
            # hardcode the CDS and exon
            if feature not in ('CDS', 'exon'):
                continue
            attributes = dict([x.split('=') for x in attributes.split(';')])
            names = attributes['Parent'].split(',')
            for name in names:
                if name not in gff:
                    gff[name] = {}
                if feature not in gff[name]:
                    gff[name][feature] = []
                gff[name][feature].append((start, end))

    # Keep track of exon boundaries
    exons = {}
    if args.exon_psls:
        for psl in args.exon_psls:
            exons = parse_psl(psl, 'exon', gff, gff_cols, exons=exons)
    if args.cds_psls:
        for psl in args.cds_psls:
            parse_psl(psl, 'CDS', gff, gff_cols)