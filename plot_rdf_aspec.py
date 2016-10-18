#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 10 15:55:38 2016

@author: Petr Voytsik
"""

import argparse
import matplotlib.pyplot as plt
import numpy as np
import os.path
from subprocess import Popen, PIPE
import sys


def plot_spec(data, chann, int_num, file_name, show=False):
    """
    """
    freq = []
    ch1 = []
    ch2 = []
    ch3 = []
    ch4 = []

    for line in data.split('\n'):
        if not line:
            continue

        if line[0] == '#':
            title = line[1:].strip()
        else:
            cols = line.split()

            try:
                freq.append(float(cols[0]))
                ch1.append(float(cols[1]))
                ch2.append(float(cols[2]))
                ch3.append(float(cols[3]))
                ch4.append(float(cols[4]))
            except ValueError:
                print(line, file=sys.stderr)
                return

    fig, ax = plt.subplots()

    ax.set_ymargin(0.05)
    ax.xaxis.set_ticks(np.arange(-16, 16+1, 4))

    ax.plot(freq, ch1, label='USB')
    ax.plot(-np.asarray(freq), ch2, label='LSB')
    ax.plot(freq, ch3, label='USB')
    ax.plot(-np.asarray(freq), ch4, label='LSB')
    ax.set_title(title, size=16)
    ax.set_xlabel('Frequency (MHz)', size=16)
    ax.set_ylabel('Amplitude', size=16)
    ax.set_ylim(ymin=0)
    ax.set_xlim(-16, 16)

    int_time = 2 * 31.25e-9 * chann * int_num
    ax.text(0.02, 0.95, 'Integration time: {:.3g} s'.format(int_time),
            transform=ax.transAxes)

    ax.legend(loc='lower center', ncol=4)
    ax.grid()

    fig.tight_layout()
    if show:
        plt.show()
    else:
        file_name = os.path.basename(file_name) + '.pdf'
        plt.savefig(file_name)

    plt.close()


def mk_spec(file_name, chann, int_num, show=False):
    """
    """
    cmd = ['rdf_aspec', file_name, str(chann), str(int_num)]

    print('Info: start processing file ', file_name)
    with Popen(cmd, stdout=PIPE, stderr=PIPE, universal_newlines=True) as aspec:
        outs, errs = aspec.communicate()
        if aspec.returncode:
            print('Error: rdf_aspec failed at file {}'.format(file_name),
                  file=sys.stderr)
            print(errs, file=sys.stderr)
            return

    plot_spec(outs, chann, int_num, file_name, show=show)


def main(args):
    """Main"""
    int_time = 2 * 31.25e-9 * args.c * args.i
    print('Info: integration time = {:.3g}'.format(int_time))

    for file_name in args.file_list:
        if not os.path.isfile(file_name):
            print('Error: file {} does not exist'.format(file_name),
                  file=sys.stderr)
            continue

        mk_spec(file_name, args.c, args.i, show=args.show)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('file_list', metavar='FILE', nargs='+',
                        help='input RDF file(s)')

    parser.add_argument('-c', type=int, default=1024,
                        help='number of spectral channels')
    parser.add_argument('-i', type=int, default=10000,
                        help='number of spectra to integrate')
    parser.add_argument('--show', action='store_true',
                        help='show plot on display')

    main(parser.parse_args())
