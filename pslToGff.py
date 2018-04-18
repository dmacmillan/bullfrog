import argparse
import os
import sys
import logging
from itertools import groupby
from customFunctions import *
from gff import *

def parse_psl(psl, feature, exons = {}):
  # Keep track of Qnames so as to avoid looking at
  # secondary alignments
  to_analyze = {}
  with open(psl, 'r') as f:
    # Read 5 lines to skip the header
    for i in range(5):
      f.readline()
    for line in f:
      cols = line.strip().split('\t')
      match = int(cols[0])
      mismatches = int(cols[1])
      cols[0] = match
      # This mRNA name is derived from
      # the maker
      mrna_name = cols[9]
      scaffold_name = cols[13]
      blocksizes = [int(x) for x in cols[18].strip(',').split(',')]
      nblocks = len(blocksizes)
      if (args.mismatch_tolerance) and (mismatches > args.mismatch_tolerance):
        skip[mrna_name] = True
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
      if mrna_name in skip:
        continue
      cols = to_analyze[mrna_name]['cols']
      scaffold_name = to_analyze[mrna_name]['scaffold_name']
      # psl uses a 0-based coordinate system
      # while gff uses a 1-based
      tstarts = [int(x)+1 for x in cols[20].strip(',').split(',')]
      blocksizes = [int(x) for x in cols[18].strip(',').split(',')]
      strand = cols[8]
      # Why is this condition here? Removing...
      # if exons or exons == {}:
        # Get the gene name
      gene_name = mrna_name.rsplit('-',2)[0]
      if not exons:
        exons = {
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
      elif gene_name not in exons:
        exons[gene_name] = {
          'scaffold_name': scaffold_name,
          'mrna_name': {
            mrna_name: {
              'start': tstarts[0],
              'end': tstarts[-1] + blocksizes[-1],
              'strand': strand
            }
          }
        }
      elif mrna_name not in exons[gene_name]['mrna_name']:
        exons[gene_name]['mrna_name'][mrna_name] = {
          'start': tstarts[0],
          'end': tstarts[-1] + blocksizes[-1],
          'strand': strand
        }
      else:
        exons[gene_name]['mrna_name'][mrna_name]['start'] = min(
          exons[gene_name]['mrna_name'][mrna_name]['start'],
          tstarts[0]
        )
        exons[gene_name]['mrna_name'][mrna_name]['end'] = max(
          exons[gene_name]['mrna_name'][mrna_name]['end'],
          tstarts[-1] + blocksizes[-1]
        )
      for i,start in enumerate(tstarts):
        out = Gff()
        out.source = 'pslToGff'
        out.feature = feature
        out.start = tstarts[i]
        out.end = tstarts[i] + blocksizes[i] - 1
        out.score = '0'
        # These boundaries can be fixed in this case
        out.seqname = scaffold_name
        out.strand = strand
        out.frame = '0'
        attributes = {
          'ID': '{}:{}:{}'.format(mrna_name, feature, i),
          'Parent': mrna_name
        }
        out.attribute = attributes
        print(out)
  return exons

parser = argparse.ArgumentParser(
  description='Given psl files, return a gff file',
  formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

# Had to use a file of files as input to avoid too many arguments erryr
parser.add_argument('-epl', '--exon_psls_list', help='A file containing a list of exon psl file paths')
parser.add_argument('-cpl', '--cds_psls_list', help='A file containing a list of cds psl file paths')
parser.add_argument(
  '-mmt',
  '--mismatch_tolerance',
  type=int,
  help='The maximum number of mismatches to tolerate per transcript. If any feature of the' \
  ' transcript has mismatches > this value the entire transcript will not be output.'
)
parser.add_argument("-l", "--log", default='warning', choices=('debug', 'info', 'warning', 'error', 'critical'), help='set the logging level. default = "warning"')

args = parser.parse_args()

if (not args.exon_psls_list) or (not args.cds_psls_list):
  raise parser.error('You must provide at least one set of features')

# Logging
# logging.basicConfig(
#   filename = os.path.join(os.getcwd(), '{}.log'.format(os.path.basename(__file__))),
#   level = getattr(logging, args.log.upper()),
#   filemode = 'w'
# )

# Keep track of gene sizes
genes = {}

# Keep track of exon boundaries
exons = {}

# Keep track of transcripts to skip
skip = {}

with open(args.exon_psls_list, 'r') as f:
  for psl in f:
    psl = psl.strip()
    exons = parse_psl(
      psl,
      'exon',
      exons = exons
    )

with open(args.cds_psls_list, 'r') as f:
  for psl in f:
    psl = psl.strip()
    parse_psl(
      psl,
      'CDS'
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
    'ID': gene
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
      'Parent': gene
    }
    out.attribute = attributes
    print(out)
