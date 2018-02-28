import pysam
import argparse
import logging
import os
import sys
import re

def gmap_attr_to_dict(attributes):
    return dict([x.split('=') for x in attributes.split(';')])

def group_gmap_gff(gff, features=None, use_name=False):
    groups = {}
    for entry in gff:
        if features and entry.feature not in features:
            continue
        attributes = gmap_attr_to_dict(entry.attributes)
        try:
            name = attributes['Parent'].split(',')
        except KeyError:
            print(entry)
            print(attributes)
            sys.exit()
        if use_name:
            name = [attributes['Name']]
        for parent in name:
            if parent not in groups:
                groups[parent] = []
            groups[parent].append(entry)
    return groups

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Given two gff files, compare blocks and block sizes',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument('pslToGff', help='Gff file generated from pslToGff')
    parser.add_argument('makerGff', help='Gff file generated by Maker')
    parser.add_argument('-f', '--features', nargs='+', default=('exon', 'CDS'), help='Only look at these features')
    parser.add_argument("-l", "--log", dest="logLevel", default='warning', choices=['debug', 'info', 'warning', 'error', 'critical'], help='set the logging level. default = "warning"')
    
    args = parser.parse_args()
    
    # Logging
    logging.basicConfig(
        filename = os.path.join(os.getcwd(), 'log-{}'.format(os.path.basename(__file__))),
        level = getattr(logging, args.logLevel.upper()),
        filemode = 'w'
    )

    # Load GFF files
    pslToGff_fo = open(args.pslToGff, 'r')
    pslToGff = pysam.tabix_iterator(pslToGff_fo, pysam.asGTF())
    makerGff_fo = open(args.makerGff, 'r')
    makerGff = pysam.tabix_iterator(makerGff_fo, pysam.asGTF())

    # Group gffs by name
    pslToGff = group_gmap_gff(pslToGff, features=args.features, use_name=True)
    makerGff = group_gmap_gff(makerGff, features=args.features)

    names_in_pslToGff_not_makerGff = pslToGff.keys() - makerGff.keys()
    names_in_makerGff_not_pslToGff = makerGff.keys() - pslToGff.keys()

    cols = [
        'Transcripts in pslToGff not makerGff',
        'Transcripts in makerGff not pslToGff',
    ]

    with open(os.path.join(os.getcwd(), 'names_in_pslToGff_not_makerGff'), 'w') as o:
        for name in names_in_pslToGff_not_makerGff:
            o.write(name + '\n')

    with open(os.path.join(os.getcwd(), 'names_in_makerGff_not_pslToGff'), 'w') as o:
        for name in names_in_makerGff_not_pslToGff:
            o.write(name + '\n')