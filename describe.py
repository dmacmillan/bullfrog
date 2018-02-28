import sys
from scipy import stats
from customFunctions import dict_to_string
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description = 'Return statistics based on a list of numbers. Can accept from stdin',
        formatter_class = argparse.ArgumentDefaultsHelpFormatter
    )
    
    parser.add_argument(
        '-i',
        '--infile',
        help = 'Read from file input'
    )
    
    args = parser.parse_args()

    data = []
    
    stream = None
    if args.infile:
        stream = open(args.infile, 'r')
    else:
        stream = sys.stdin
    for line in stream:
        data.append(float(line.strip()))
    
    print(dict_to_string(compute_stats(data)))
