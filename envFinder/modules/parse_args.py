
import argparse

def parseArgs():
    parser = argparse.ArgumentParser(prog='envFinder')

    parser.add_argument('--env_co', default = 0.8, type = float,
                        help = 'Envelope corelation score cuttoff.')

    parser.add_argument('--mz_step_margin', default = 2, type= int,
                        help = 'Margin above and below envelope in plot.')

    parser.add_argument('filter_file', nargs='+',
                        help='DTAFilter file to read.')


    args = parser.parse_args()

    return args
