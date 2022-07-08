
import sys
import os
import argparse
import subprocess
import os.path
import re
from math import ceil

from .modules import parent_parser

ENV_FINDER_EXE = 'envoMatch'

def makePBS(mem, ppn, walltime, wd, args, shell='tcsh'):
    pbsName = '{}/{}.pbs'.format(wd, ENV_FINDER_EXE)
    _flags = ' '.join(['--{} {}'.format(k,v) for k, v in args.items() if v is not None and k != 'input_file'])

    sys.stdout.write('Writing {}...'.format(pbsName))
    with open(pbsName, 'w') as outF:
        outF.write("#!/usr/bin/env {}\n".format(shell))
        outF.write('#PBS -l mem={}gb,nodes=1:ppn={},walltime={}\n\n'.format(mem, ppn, walltime))
        outF.write('cd {}\n'.format(wd))
        outF.write('{} {} {} > stdout.txt\n'.format(ENV_FINDER_EXE, _flags, args['input_file']))

    sys.stdout.write('Done!\n')
    return pbsName


def makeLSF(mem, ppn, walltime, wd, args, shell='bash'):
    # parse walltime for lsf to hh:mm format
    match=re.search(r'^(\d+):(\d{1,2}):(\d{1,2})$', walltime)
    if match:
        _walltime='{}:{}'.format(match.group(1), match.group(2))
    else:
        raise RuntimeError('Could not parse walltime: "{}"'.format(walltime))    

    lsfName = '{}/{}.lsf'.format(wd, ENV_FINDER_EXE)
    _flags = ' '.join(['--{} {}'.format(k,v) for k, v in args.items() if v is not None and k != 'input_file'])

    sys.stdout.write('Writing {}...'.format(lsfName))
    with open(lsfName, 'w') as outF:
        outF.write("#!/usr/bin/env {}\n".format(shell))
        outF.write('#BSUB -W {} -n {} -R "rusage[mem={}]" -R span[hosts=1] -q short -J {}\n'.format(_walltime, ppn, mem * 1024, lsfName))
        outF.write('#BSUB -o stdout.%J.%I -e stderr.%J.%I\n\n')
        outF.write('cd {}\n'.format(wd))
        outF.write('{} {} {} > stdout.txt\n'.format(ENV_FINDER_EXE, _flags, args['input_file']))
    sys.stdout.write('Done!\n')

    return lsfName


JOB_TYPES = {'pbs': {'makeJob': makePBS, 'qsub': 'qsub'},
             'lsf': {'makeJob': makeLSF, 'qsub': 'bsub <'}}


def main():
    parser = argparse.ArgumentParser(prog='qsub_envoMatch', parents=[parent_parser.PARENT_PARSER],
                                     description='Submit {} job to the queue.'.format(ENV_FINDER_EXE))

    parser.add_argument('-m', '--mem', default=8, type=int,
                        help='Amount of memory to allocate per job in gb. Default is 4.')

    parser.add_argument('-p', '--ppn', default=8, type=int,
                        help='Number of processors to allocate per job. Default is 8.')

    parser.add_argument('--jobType', choices=['pbs', 'lsf'], default='pbs',
                        help='Job type. Default is "pbs"')

    parser.add_argument('--shell', default='tcsh', help='The shell to use in job files')

    parser.add_argument('-w', '--walltime', default='12:00:00',
                        help='Walltime per job in the format hh:mm:ss. Default is 12:00:00.')

    parser.add_argument('-g', '--go', action='store_true', default=False,
                        help='Should job be submitted? If this flag is not supplied, program will be a dry run. '
                             'job file will written but job will not be submitted.')

    args = parser.parse_args()
    parent_args = parent_parser.PARENT_PARSER.parse_known_args()[0]

    # calc nThread
    if args.parallel and args.nThread is None:
        _nThread = args.ppn * 2
    elif not args.parallel and args.nThread is None:
        _nThread=1
    else:
        _nThread = args.nThread
        args.ppn = ceil(args.nThread / 2) * 2

    # Manually check args
    n_arg_errors = 0
    error_message = '\nThere were errors in the options you specified...'
    if not os.path.isfile(args.input_file):
        error_message += '\n\tSpecified input file: {} does not exist on path:\n\t\t{}\n'.format(args.input_file,
                                                                                                 os.path.abspath(args.input_file))
        n_arg_errors += 1
    if n_arg_errors > 0:
        sys.stderr.write(error_message)
        return -1

    sys.stdout.write('\nRequested job with {} processor(s) and {}gb of memory...\n'.format(args.ppn, args.mem))
    #get wd
    wd = os.path.dirname(os.path.abspath(args.input_file))

    envoMatch_args = {arg: getattr(args, arg) for arg in vars(parent_args)}
    _nThread = args.nThread
    envoMatch_args['nThread'] = _nThread
    envoMatch_args['plotEnv'] = '' if args.plotEnv else None
    envoMatch_args['splitPlots'] = '' if args.splitPlots else None
    envoMatch_args['verbose'] = '' if args.verbose else None

    jobName = JOB_TYPES[args.jobType]['makeJob'](args.mem, args.ppn, args.walltime, wd, envoMatch_args, shell=args.shell)
    command = '{} {}'.format(JOB_TYPES[args.jobType]['qsub'], jobName)
    if args.verbose and args.go:
        sys.stdout.write('{}\n'.format(command))
    if args.go:
        proc = subprocess.Popen([command], cwd=wd, shell=True)
        proc.wait()

if __name__ == '__main__':
    main()

