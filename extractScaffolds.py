import argparse
from gff import Gff

parser = argparse.ArgumentParser(
  description='Given a gff file, return a gff file where' \
  ' the ID attribute contains the seqname',
  formatter_class=argparse.ArgumentDefaultsHelpFormatter
)

parser.add_argument('gff', help='Input gff file')

args = parser.parse_args()

for line in Gff.read(args.gff):
  if line.seqname in line.attribute['ID']:
    print(line)
