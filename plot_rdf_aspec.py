#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 15:55:38 2016

@author: Petr Voytsik
"""

import argparse
import os.path
from subprocess import Popen, PIPE
import sys

import matplotlib.pyplot as plt
import numpy as np


def plot_spec(file_name, chann, int_num, max_amp=None, show=False):
    """
    """
    freq = []
    ch1 = []
    ch2 = []
    ch3 = []
    ch4 = []

    with open(file_name, 'r') as file:
        title = file.readline()[1:].strip()

    freq, ch1, ch2, ch3, ch4 = np.loadtxt(file_name, unpack=True)

    fig, ax = plt.subplots()

    ax.set_ymargin(0.05)
    if not show:
        ax.xaxis.set_ticks(np.arange(-16, 16+1, 4))

    ax.plot(freq, ch1, label='USB1')
    ax.plot(-freq, ch2, label='LSB1')
    ax.plot(freq, ch3, label='USB2')
    ax.plot(-freq, ch4, label='LSB2')
    ax.set_title(title, size=16)
    ax.set_xlabel('Frequency (MHz)', size=16)
    ax.set_ylabel('Amplitude', size=16)
    ax.set_ylim(ymin=0)
    ax.set_xlim(-16, 16)
    if max_amp:
        ax.set_ylim(ymax=max_amp)

    int_time = 2 * 31.25e-9 * chann * int_num
    ax.text(0.02, 0.95, 'Integration time: {:.4g} s'.format(int_time),
            transform=ax.transAxes)

    # ax.legend(loc='lower center', ncol=4)
    ax.legend(loc='best', ncol=2)
    ax.grid()

    fig.tight_layout()
    if show:
        plt.show()
    else:
        file_name = os.path.basename(file_name) + '.pdf'
        plt.savefig(file_name)

    plt.close()


def mk_spec(file_name, chann, int_num, offset=0):
    """
    Generete autospectrum using ``rdf_aspec``.

    """
    cmd = ['rdf_aspec', file_name, str(chann), str(int_num), str(offset)]

    print('Info: start processing file ', file_name)
    with Popen(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True) \
            as aspec:
        outs, errs = aspec.communicate()
        if aspec.returncode:
            print('Error: rdf_aspec failed at file {}'.format(file_name),
                  file=sys.stderr)
            print(errs, file=sys.stderr)
            return None

    return outs


def main(args):
    """Main"""
    # Check input parameters
    if args.integration_time <= 0:
        print('Error: Intergration time must be positive', file=sys.stderr)
        return 1
    if args.c <= 0:
        print('Error: The number of spectral channels must be positive',
              file=sys.stderr)
        return 1
    if args.offset < 0:
        print('Error: The offset must not be negative', file=sys.stderr)
        return 1

    int_time = args.integration_time
    spec_count = int(int_time / (2 * 31.25e-9 * args.c))
    print('Debug: number of FFT frames = {}'.format(spec_count))

    start_off_list = np.arange(args.offset, 1200-int_time, int_time)[:args.n]

    for file_name in args.file_list:
        if not os.path.isfile(file_name):
            print('Error: file {} does not exist'.format(file_name),
                  file=sys.stderr)
            continue

        # offset = args.offset
        for offset in start_off_list:
            print('Info: offset = {:.1f} s'.format(offset))
            basename = os.path.basename(file_name).split('.')[0]
            out_file_name = '{}_{:04d}.aspec'.format(basename, int(offset))

            if not os.path.isfile(out_file_name) or args.force:
                spec_data = mk_spec(file_name, args.c, spec_count,
                                    offset=offset)
                if spec_data:
                    with open(out_file_name, 'w') as out_file:
                        out_file.write(spec_data)
                else:
                    print('Error: autospec data is empty', file=sys.stderr)
                    continue
            else:
                print('Note: file {} already exists'.format(out_file_name))

            plot_spec(out_file_name, args.c, spec_count, max_amp=args.max_ampl,
                      show=args.show)

    return 0


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('file_list', metavar='FILE', nargs='+',
                        help='input RDF file(s)')

    parser.add_argument('-c', type=int, default=1024,
                        help='number of spectral channels')
#    parser.add_argument('-i', type=int, default=10000,
#                        help='number of spectra to integrate')
    parser.add_argument('-i', '--integration-time', type=float, default=10,
                        help='integration time (s)')
    parser.add_argument('-o', '--offset', type=float, default=0,
                        help='offset is seconds')
    parser.add_argument('-n', type=int, default=1,
                        help='generate N sequence autospectrum')
    parser.add_argument('--show', action='store_true',
                        help='show plot on display')
    parser.add_argument('--max-ampl', type=float,
                        help='amplitude upper limit on plot')
    parser.add_argument('-f', '--force', action='store_true',
                        help='force update .aspec file')

    sys.exit(main(parser.parse_args()))
