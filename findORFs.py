#!/usr/bin/env python3
# Name: Rob Lodes (rlodes)
# Group Members: got help from most of the tutors for this lab
'''
Analyzes a FASTA-formatted file containing a sequence of DNA, and finds the ORFs.
Input: a FASTA file
Output: formatted output listing the seq name, and for each ORFs the frame, start, stop and legngth. The frame
        is noted with a +/- to illustrate top/bottom reading of the strand.
Example: python findORFs.py < data.fa -mG=100 -lG > output.txt
    -lG, --longestGene, action='store', nargs='?', const=True, default=False,
                                 help='longest Gene in an ORF'
    -mG, --minGene, type=int, choices=(0, 100, 200, 300, 500, 1000), default=100,
                                 action='store', help='minimum Gene length'
    -s, --start, action='append', default=['ATG'], nargs='?',
                                 help='start Codon'
    -t, --stop, action='append', default=['TAG', 'TGA', 'TAA'], nargs='?',
                                 help='stop Codon'
    -v, --version, action='version', version='%(prog)s 0.1'
'''
from sequenceAnalysis import OrfFinder, FastAreader

########################################################################
# CommandLine
########################################################################
class CommandLine():
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.
    '''

    def __init__(self, inOpts=None):
        '''
        Implement a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(
            description='Program prolog - a brief description of what this thing does',
            epilog='Program epilog - some other stuff you feel compelled to say',
            add_help=True,  # default is True
            prefix_chars='-',
            usage='%(prog)s [options] -option1[default] <input >output'
            )
        self.parser.add_argument('-lG', '--longestGene', action='store', nargs='?', const=True, default=False,
                                 help='longest Gene in an ORF')
        self.parser.add_argument('-mG', '--minGene', type=int, choices=(0, 100, 200, 300, 500, 1000), default=100,
                                 action='store', help='minimum Gene length')
        self.parser.add_argument('-s', '--start', action='append', default=None, nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-t', '--stop', action='append', default=None, nargs='?',
                                 help='stop Codon')  # allows multiple list optionsca
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
########################################################################

def main(inFile=None, options=None):
    '''
    Find some genes.
    '''
    thisCommandLine = CommandLine(options)

    if (thisCommandLine.args.start is None):
        thisCommandLine.args.start = ['ATG']
    if (thisCommandLine.args.stop is None):
        thisCommandLine.args.stop = ['TAG', 'TGA', 'TAA']

    # myFasta = FastAreader(inFile)                     # use this for debugging.
    myFasta = FastAreader()                             # use this one for the command line (also to turn in)

    for head, seq in myFasta.readFasta():
                                                        # preprocess sequence
        seq = seq.upper().replace(" ", "")              # upper case, and remove white space
        print(head)                                     # print the header
        myOrf = OrfFinder(seq,                          # instantiate, call, and send commands to OrfFinder
                          thisCommandLine.args.start,
                          thisCommandLine.args.stop,
                          thisCommandLine.args.longestGene,
                          thisCommandLine.args.minGene)
                                                        # sort the ORFs by reverse length, not rev frame
        sortedOrf = sorted(myOrf.findStrandedOrfs(), key = lambda x: (x[3], -x[1]), reverse = True)
                                                        # x3 is length, -x1 is start, neg reverses to sort accecnding
        for orf in sortedOrf:                           # print out the sorted ORFs
            # print('{:+d} {:>5d}..{:>5d} {:>5d} {}'.format(orf[0], orf[1], orf[2], orf[3], orf[4])) # debug version
            print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(orf[0], orf[1], orf[2], orf[3]))

    # print(thisCommandLine.args)                       # for debugging

if __name__ == "__main__":
    # main(inFile='tass2.fa', options=['-mG=300', '-lG'])   # for debugging
    # main(inFile='lab5test.fa', options=['-mG=0', '-lG'])  # for debugging
    main()  # use this if running from commandline (and for turning in)
