
import argparse

def parseArgs():
    parser = argparse.ArgumentParser(prog='envFinder')

    parser.add_argument('filter_file', nargs='+',
                        help='DTAFilter file to read.')

    args = parser.parse_args()

    return args
