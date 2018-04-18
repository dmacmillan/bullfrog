import argparse
import os
import sys
from itertools import groupby
from customFunctions import *
from gff import *
from psl import *
import re
import time
import pickle
import pysam
import logging
import copy

def groupGff(
    gff_file,
    is_gmap = False
):
    def process(string, num=1):
        return string.rsplit('.', num)[0]
    def blank(string, *args):
        return string
    do = blank
    if is_gmap:
        do = process
    genes = {}
    mrnas = {}
    features = {}
    for entry in Gff.read(gff_file):
        if entry.feature == 'gene':
            gene_id = entry.attribute['ID']
            if is_gmap:
                gene_id = gene_id.rsplit('-', 2)[0]
            genes[gene_id] = {
                'gene': entry,
                'mRNAs': {}
            }
        elif entry.feature == 'mRNA':
            mrna_id = do(entry.attribute['ID'])
            mrnas[mrna_id] = {
                'mRNA': entry
            }
        else:
            for parent in entry.attribute['Parent'].split(','):
                parent = do(parent)
                if parent not in features:
                    features[parent] = [entry]
                else:
                    features[parent].append(entry)
    for mrna_id in features:
        try:
            mrnas[mrna_id]['features'] = features[mrna_id]
        except KeyError:
            continue
    for mrna_id in mrnas:
        try:
            gene_id = mrnas[mrna_id]['mRNA'].attribute['Parent']
            if is_gmap:
                gene_id = gene_id.rsplit('-', 2)[0]
            genes[gene_id]['mRNAs'][mrna_id] = mrnas[mrna_id]
        except KeyError:
            continue
    return genes

def getProtein(cdss, fasta, rev=False, trunc=True):
    parts = []
    for x in cdss:
        seq = fasta.fetch(x.seqname, x.start, x.end)
        parts.append(seq)
    # seq = ''.join([fasta.fetch(x[0], x[1], x[2]) for x in cdss])
    seq = ''.join(parts)
    if rev:
        return parts, trans(rev_comp(seq), truncate=trunc)
    return parts, trans(seq, truncate=trunc)

def getExonSeqs(exons, fasta):
    parts = []
    for x in exons:
        seq = fasta.fetch(x.seqname, x.start, x.end)
        parts.append(seq)
    return parts

def mapBlocks(scaffold, scaffold_sequence, features, feature_seqs):
    new_features = []
    for i, feature in enumerate(features):
        match = None
        try:
            match = scaffold_sequence.index(feature_seqs[i])
        except ValueError:
            logging.warning('Sequence {} not found in {}. Feature is {}'.format(feature_seqs[i], scaffold, feature))
            # try:
            #     match = scaffold_sequence.index(rev_comp(feature_seqs[i]))
            # except ValueError:
            #     logging.warning('  Reverse comp not found either')
            return None
        new_feature = copy.deepcopy(features[i])
        new_feature.start = match + 1
        new_feature.end = new_feature.start + feature.end - feature.start
        new_feature.seqname = scaffold
        new_features.append(new_feature)
    return new_features

# block and truth_block are tuples of the form: (scaffold_name (string), cds_start (int), cds_end (int))
# block_seq and truth_block_seq are simply strings of nucleotides
# Will attempt to make block equal to truth_block in terms of size
# Cannot determine if the block should be extended or shortened on the right or left
# without comparing the sequence of the two blocks
# Will return the amended block if they are similar, otherwise returns false
def areBlocksSimilar(block, block_seq, truth_block, truth_block_seq):
    logging.debug('block: "{}"'.format(block))
    logging.debug('block_seq: "{}"'.format(block_seq))
    logging.debug('truth_block: "{}"'.format(truth_block))
    logging.debug('truth_block_seq: "{}"'.format(truth_block_seq))
    truth_block_size = truth_block[2] - truth_block[1]
    block_size = block[2] - block[1]
    diff = truth_block_size - block_size
    # Three amino acids
    if abs(diff) <= 9:
        # There are 3 cases to consider
        # Case 1: Truth sequence has extra nucleotides
            # Two subcases for case 1
            # Case 1A: Extra sequence on the left-most end
            # Case 1B: Extra sequence on the right-most end
        # Case 2: Truth sequence has missing nucleotides
            # Two subcases for case 2
            # Case 2A: Missing sequence on the left-most end
            # Case 2B: Missing sequence on the right-most end
        # Case 3: Sequences are of the same length
        # Case 1
        if truth_block_size > block_size:
            logging.debug('Case 1')
            # Case 1A
            if truth_block_seq[diff:] == block_seq:
                logging.debug('Case 1A')
                return (block[0], block[1] - diff, block[2])
            # Case 1B
            else:
                logging.debug('Case 1B')
                return (block[0], block[1], block[2] + diff)
        # Case 2
        elif block_size > truth_block_size:
            logging.debug('Case 2')
            # Case 2A
            if truth_block_seq == block_seq[diff:]:
                logging.debug('Case 2A')
                return (block[0], block[1] + diff, block[2])
            # Case 2B
            else:
                logging.debug('Case 2B')
                return (block[0], block[1], block[2] - diff)
        # Case 3
        else:
            logging.debug('Case 3')
            if block_seq != truth_block_seq:
                return False
            return truth_block
    return False

if __name__ == '__main__':

    parser = argparse.ArgumentParser(
        description='Compare two or more GFF files',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # This could be done better by taking in any number of gffs and providing a labels argument to label each one
    parser.add_argument('-t', '--truth', help='All other GFF files will only be compared to this GFF file')
    parser.add_argument('-tf', '--truth_fasta', help='Fasta file for which the truth GFF describes')
    parser.add_argument('-g', '--gffs', nargs='+', help='GFF files to compare')
    parser.add_argument('-f', '--fastas', nargs='+', help='FASTA files for which the GFF files describe')
    parser.add_argument('-a', '--alignments', nargs='+', help='PSL alignment files for each GFF transcript CDS. Files must have specific names according to labels. I.e. two gff files with labels: "foo" and "bar" would need an alignment file beginning with "foo_bar" if "foo" were the target and "bar" were the query, in other words "bar" is aligned to "foo"')
    parser.add_argument('-l', '--labels', nargs='+', help='A name for each respective GFF file')
    parser.add_argument('-tl', '--transcripts_list', help='A list of transcripts, only these transcripts will be compared')
    # parser.add_argument('-if', '--include_features', nargs='+', help='Only these features will be analyzed')
    # parser.add_argument('-ef', '--exclude_features', nargs='+', help='All features except for these will be analyzed')
    parser.add_argument("--loglevel", default='warning', choices=('debug', 'info', 'warning', 'error', 'critical'), help='set the logging level. default = "warning"')

    args = parser.parse_args()

    try:
        truth_fasta = pysam.FastaFile(args.truth_fasta)
    except TypeError:
        truth_fasta = None

    logging.basicConfig(
        level = getattr(logging, args.loglevel.upper()),
        format = '%(message)s',
        filename = 'compareGffs.log',
        filemode = 'w'
    )

    logging.warning('TEST')

    if (not args.truth) or (not args.gffs):
        parser.error('You must specify at least two gffs')
    elif (not args.labels):
        parser.error('You must provide labels for your gffs')
    # elif (not args.alignments):
    #     parser.error('You must provide alignments for your gffs')

    transcripts_list = set()
    with open(args.transcripts_list, 'r') as f:
        for line in f:
            transcripts_list.add(line.strip())

    # Store data
    truth = None
    data = {}

    if args.truth:
        fpath = './.truth.pickle'
        if os.path.isfile(fpath):
            logging.debug('Loading {} ...'.format(fpath))
            truth = pickle.load(open(fpath, 'rb'))
        else:
            logging.debug('Parsing truth ...')
            truth = groupGff(
                args.truth
            )
            logging.debug('DONE')
            pickle.dump(truth, open(fpath, 'wb'))

    for i,gff in enumerate(args.gffs):
        label = args.labels[i]
        is_gmap = False
        if label.lower() == 'gmap':
            is_gmap = True
        fpath = './.{}.pickle'.format(label)
        if os.path.isfile(fpath):
            logging.debug('Loading {} ...'.format(fpath))
            data[label] = pickle.load(open(fpath, 'rb'))
        else:
            logging.debug('Parsing {} ...'.format(label))
            data[label] = {
                'gff': groupGff(
                    gff,
                    is_gmap
                )
            }
            logging.debug('DONE')
            pickle.dump(data[label], open(fpath, 'wb'))
        try:
            data[label]['fasta'] = pysam.FastaFile(args.fastas[i])
        except TypeError:
            logging.warning('No fasta specified for {}'.format(label))
            data[label]['fasta'] = None

    if truth:
        for gene in truth:
            # logging.debug('gene: "{}"'.format(gene))
            for transcript in truth[gene]['mRNAs']:
                truth_features = truth[gene]['mRNAs'][transcript]['features']
                used = False
                if transcript not in transcripts_list:
                    continue
                # print('transcript: {}'.format(transcript))
                truth_cdss = []
                truth_exons = []
                truth_nblocks = 0
                truth_ncds = 0
                truth_blockSizes = []
                truth_is_reverse = False
                for feature in truth_features:
                    if feature.strand == '-':
                        truth_is_reverse = True
                    if feature.feature == 'exon':
                        truth_nblocks += 1
                        truth_blockSizes.append((feature.end - feature.start))
                        truth_exons.append(feature)
                    elif feature.feature == 'CDS':
                        truth_ncds += 1
                        # start -1 since this is fed into pysam fetch
                        truth_cdss.append(feature)
                truth_blockSizes = sorted(truth_blockSizes)
                truth_cdss = sorted(truth_cdss, key = lambda x: x.start)
                truth_cds_sizes = [x.end - x.start for x in truth_cdss]
                truth_parts, truth_protein = getProtein(truth_cdss, truth_fasta, truth_is_reverse)
                truth_exon_seqs = getExonSeqs(truth_exons, truth_fasta)
                # subset = False
                for key in data:
                    if used:
                        continue
                    try:
                        data[key]['gff'][gene]
                    except KeyError:
                        logging.debug('{} not in {}'.format(gene, key))
                        logging.debug(data[key]['gff'].keys())
                        sys.exit()
                        continue
                    try:
                        features = data[key]['gff'][gene]['mRNAs'][transcript]['features']
                    except KeyError:
                        logging.debug('"{}" not in "{}"'.format(transcript, key))
                        logging.debug('truth[{}]: "{}"'.format(gene, truth[gene]))
                        logging.debug('{}[{}]: "{}"'.format(key, gene, data[key]['gff'][gene]))
                        data[key]['gff'][gene]['mRNAs'][transcript] = {
                            'nblocks': -1
                        }
                        continue
                    cdss = []
                    nblocks = 0
                    ncds = 0
                    blockSizes = []
                    exons = []
                    is_reverse = False
                    for feature in features:
                        if feature.strand == '-':
                            is_reverse = True
                        if feature.feature == 'exon':
                            nblocks += 1
                            blockSizes.append((feature.end - feature.start))
                            exons.append(feature)
                        elif feature.feature == 'CDS':
                            ncds += 1
                            # start -1 since this is fed into pysam fetch
                            cdss.append(feature)
                    blockSizes = sorted(blockSizes)
                    cdss = sorted(cdss, key = lambda x: x.start)
                    v3scaffold = cdss[0].seqname
                    cds_sizes = [x.end - x.start for x in cdss]
                    can_fix = False
                    diff_blocks = truth_nblocks - nblocks
                    new_cdss = cdss[:]
                    seq = data[key]['fasta'].fetch(v3scaffold)
                    seqlen = len(seq)
                    if is_reverse != truth_is_reverse:
                        logging.warning('Transcript "{}" is not of the same strand in {} compared to truth'.format(transcript, key))
                    try:
                        parts, protein = getProtein(cdss, data[key]['fasta'], is_reverse)
                        exon_seqs = getExonSeqs(exons, data[key]['fasta'])
                    except TypeError:
                        logging.warning('{} does not have a fasta file specified'.format(key))
                    # print('v3 scaffold name: {}'.format(v3scaffold))
                    # print('Truth ncds: {}'.format(truth_ncds))
                    # print('{} ncds: {}'.format(key, ncds))
                    # print('Truth cds sizes: {}'.format(truth_cds_sizes))
                    # print('{} cds sizes: {}'.format(key, cds_sizes))
                    # print('Truth cdss: {}'.format(truth_cdss))
                    # print('{} cdss: {}'.format(key, cdss))
                    # new_cdss = mapBlocks(v3scaffold, seq, truth_cdss, truth_parts)
                    new_cdss = mapBlocks(v3scaffold, seq, truth_cdss, truth_parts)
                    new_exons = mapBlocks(v3scaffold, seq, truth_exons, truth_exon_seqs)
                    if not new_cdss:
                        continue
                    for cds in new_cdss:
                        print(cds)
                    if not new_exons:
                        continue
                    for exon in new_exons:
                        print(exon)
                    used = True
                    # mapBlocks(query_scaffold, scaffold_sequence, target_cdss, target_seqs)

                ##### Commented as of 2018-03-28 #####
                #     for i, truth_block in enumerate(truth_cdss):
                #         print('i: {}'.format(i))
                #         print('truth block: {}'.format(truth_block))
                #         truth_seq = truth_fasta.fetch(truth_cdss[i][0], truth_cdss[i][1], truth_cdss[i][2])
                #         truth_block_size = truth_cds_sizes[i]
                #         pop = True
                #         if i >= ncds:
                #             search_seq = seq[cdss[-1][2]:]
                #             similar = False
                #             pop = False
                #         else:
                #             similar = areBlocksSimilar(cdss[i], parts[i], truth_cdss[i], truth_parts[i])
                #         print('Are blocks similar?: {}'.format(similar))
                #         if not similar:
                #             if i >= ncds - 1:
                #                 start = cdss[-1][2]
                #                 end = seqlen
                #                 search_seq = seq[start:]
                #             else:
                #                 start = cdss[i][2]
                #                 end = cdss[i + 1][1]
                #                 search_seq = seq[start:end]
                #             match = search_seq.find(truth_seq)
                #             if match:
                #                 new_block = (v3scaffold, start + match, end)
                #                 new_cdss.insert(i, new_block)
                #                 if pop:
                #                     new_cdss.pop(i + 1)
                #     print('new_cdss: {}'.format(new_cdss))
                # print('~'*50)

                #     print('v3 scaffold name: {}'.format(v3scaffold))
                #     print('Truth nblocks: {}'.format(truth_nblocks))
                #     print('{} nblocks: {}'.format(key, nblocks))
                #     print('Truth cds sizes: {}'.format(truth_cds_sizes))
                #     print('{} cds sizes: {}'.format(key, cds_sizes))
                #     # Need to fetch sequence and translate
                #     try:
                #         parts, protein = getProtein(cdss, data[key]['fasta'], is_reverse)
                #     except TypeError:
                #         logging.warning('{} does not have a fasta file specified'.format(key))
                #     # 3 cases to consider:
                #     # Case 1: number of blocks are equal
                #     if nblocks == truth_nblocks:
                #         print('Case 1')
                #         pass
                #     # Case 2: alternate annotation has one or more additional blocks
                #     elif nblocks > truth_nblocks:
                #         print('Case 2')
                #         pass
                #     # Case 3: maker annotation has one or more additional blocks
                #     else:
                #         # Here we consider cases where the number of alternate
                #         # blocks is less than the number of maker blocks.
                #         #! It is important to note that blocks may still
                #         # vary in size!
                #         # So we must first determine if each block is different
                #         # and if so, should it be considered a new block or
                #         # should it be modified to match
                #         # For case 3 we have several sub-cases to consider
                #         # Subcase A: Extra block is the left-most block
                #         # Subcase B: Extra block is the right-most block
                #         # Subcase C: Extra block is in between two other blocks
                #         # First need to figure out which subcase this falls under
                #         # Subcase A
                #         print('Case 3')
                #         first_block_similar = areBlocksSimilar(cdss[0], parts[0], truth_cdss[0], truth_parts[0])
                #         last_block_similar = areBlocksSimilar(cdss[-1], parts[-1], truth_cdss[-1], truth_parts[-1])
                #         print('{} first block: {}'.format(key, cdss[0]))
                #         print('Truth first block: {}'.format(truth_cdss[0]))
                #         print('{} last block: {}'.format(key, cdss[-1]))
                #         print('Truth last block: {}'.format(truth_cdss[-1]))
                #         print('Are the first blocks similar?: {}'.format(first_block_similar))
                #         print('Are the last blocks similar?: {}'.format(last_block_similar))
                #         # Subcase A
                #         if last_block_similar and not first_block_similar:
                #             print('Subcase A')
                #         # Subcase B
                #         elif first_block_similar and not last_block_similar:
                #             print('Subcase B')
                #         # Subcase C
                #         else:
                #             print('Subcase C')
                # print('~'*50)
