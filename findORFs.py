#!/usr/bin/env python3
# Name: Rob Lodes (rlodes)
# Group Members: None

from sequenceAnalysis import OrfFinder, FastAreader

# def main():
#     # myReader = FastAreader('testGenome.fa')     # use this one for testing
#     myReader = FastAreader()  # Use this for final program to direct in from stdin
#
#     print(myReader)
#
# # Press the green button in the gutter to run the script.
# if __name__ == '__main__':
#     print_hi('PyCharm')


# !/usr/bin/env python3
# Name: Your full name (CATS account username)
# Group Members: List full names (CATS usernames) or “None”


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
        self.parser.add_argument('-s', '--start', action='append', default=['ATG'], nargs='?',
                                 help='start Codon')  # allows multiple list options
        self.parser.add_argument('-t', '--stop', action='append', default=['TAG', 'TGA', 'TAA'], nargs='?',
                                 help='stop Codon')  # allows multiple list options
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 0.1')
        if inOpts is None:
            self.args = self.parser.parse_args()
        else:
            self.args = self.parser.parse_args(inOpts)

########################################################################
# Main
# Here is the main program
#
#
########################################################################

def main(inFile=None, options=None):
    '''
    Find some genes.
    '''
    thisCommandLine = CommandLine(options)

    myFasta = FastAreader(inFile) # use this for debugging.
    # myFasta = FastAreader() # use this one for the command line (also to turn in)

    for head, seq in myFasta.readFasta():
        # preprocess sequence
        # upper case, and remove white space
        seq = seq.upper().replace(" ", "")
        # print(head, seq)
        print(head)
        myOrf = OrfFinder(seq,
                          thisCommandLine.args.start,
                          thisCommandLine.args.stop,
                          thisCommandLine.args.longestGene,
                          thisCommandLine.args.minGene)

        # myOrf.findStrandedOrfs()     # temp function to get things going
        print('myOrf seq ', myOrf.seq)  # test call to myOrf
        sortedOrf = sorted(myOrf.findStrandedOrfs(), key = lambda x: (x[3], -x[1]), reverse = True)
        print('sorted orfs ', sortedOrf)
        for orf in sortedOrf:
            print('{:+d} {:>5d}..{:>5d} {:>5d}'.format(orf[0], orf[1], orf[2], orf[3]))

    print(thisCommandLine.args)

    # print(thisCommandLine.args.longestGene)
    # print(thisCommandLine.args.start)
    # print(thisCommandLine.args.stop)
    # print(thisCommandLine.args.minGene)


if __name__ == "__main__":
    # main(inFile='tass2.fa', options=['-mG=300', '-lG'])  # delete this stuff if running from commandline
    main(inFile='lab5test.fa', options=['-mG=0', '-lG'])  # delete this stuff if running from commandline
    # main()  # use this stuff if running from commandline (and for turning in)
