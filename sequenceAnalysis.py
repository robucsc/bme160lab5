#!/usr/bin/env python3
# Name: Rob Lodes (rlodes)
# Group Members: Got lots of help from Stephen (the tutor) on this assignment
"""
    This file contains three classes NucParams, FastAreader, and ProteinParam

    NucParams converts FastA sequences into codon, amino acid, and nucleotide dictionaries.

    FastAreader defines objects to read FastA files.

    ProteinParam calculates the physical-chemical properties of a protein sequence
"""

class OrfFinder:
    ''' doc string '''
    def __init__(self, seq, start, stop, longestGene, minGene):
        self.seq = seq
        self.start = start
        self.stop = stop
        self.longestGene = longestGene
        self.minGene = minGene

    def makeCodon(self):
        '''
        gerantion of codons
        '''
        for frame in range(3):
            for offsetIndex in range(frame, len(self.seq), 3):
                codon = self.seq[offsetIndex:offsetIndex + 3]
                # print('frame ', frame + 1, 'codon ', codon)     # debug code, the frame +1 is for ease in reading
                yield frame, offsetIndex, codon  # yield for generator

    def findORFs(self, top):
        orfList = []
        startPos = [0]                          # start at zero to accommodate dangling starts
        myCodons = self.makeCodon()             # to allow for the generator to be nexted
        # next(myCodons)                        # skips the first zero codon to account for the generator starting at zero
                                                # use next in the previous - this step is now in the following line
        previousFrame = next(myCodons)[0]       # unpacking version , _, _ = next(myCodons)
        for frame, i, codon in myCodons:
            # check if codon is start
            if frame != previousFrame:                      # handle danging starts for frame 1 and 2
                for start in startPos:
                    # orf is frame, start, stop, length, seq-code section -- in this case stop is the same as length
                    if top:                                 # top strand - section 1
                        orfLength = len(self.seq) - start   # precalc orf length
                        orf = [frame + 0, start + 1, len(self.seq), orfLength,
                               self.seq[start:len(self.seq)] + ' 1']  # +1 offsets for alignment
                    else:                                   # bottom strand - section 2
                        orfLength = len(self.seq) - start   # precalc orf length
                        orf = [(frame + 0) * -1, 1, len(self.seq) - start, orfLength,
                               self.seq[start:len(self.seq) - start] + ' 2']
                    if (orfLength > self.minGene):          # minGene comparison, if smaller don't append
                        orfList.append(orf)                 # add orf to the list
                startPos = [0]                              # end of orf so init start position
            previousFrame = frame                           # set previousFrame for next loop

            if codon in self.start:                         # frame 0
                startPos.append(i)                          # save position
                # print("starts ", codon)                   # debug code
            elif codon in self.stop:
                for start in startPos:
                    # compare to minGene, is the length greater than minGene - not implemented yet
                    # if at start of seq.. always edge case if first stop found
                    if top:                                 # top strand - section 3
                        orfLength = i - start + 3           # precalc orf length
                        orf = [frame + 1, start + 1, i + 2 + 1, orfLength,
                               self.seq[start:i + 2 + 1] + ' 3']  # +1 offsets for alignment
                    else:       # bottom strand non-dangling, and dangling stops - section 4
                        orfLength = i - start + 3           # precalc orf length
                        orf = [(frame + 1) * -1, (len(self.seq) + 1) - (i + 3),
                               (len(self.seq) + 1) - (start + 1), orfLength,
                               self.seq[start:(len(self.seq) + 1) - (start + 1)] + ' 4']
                    if (orfLength > self.minGene):          # minGene comparison, if smaller don't append
                        orfList.append(orf)                 # add orf to the list
                startPos = []                               # end of orf so init start position

        for start in startPos:                              # frame 3 to end, only happens with no stop
            if top:                                         # top strand - section 5
                orfLength = len(self.seq) - start           # precalc orf length
                orf = [previousFrame + 1, start + 1, len(self.seq),
                       orfLength, self.seq[start:len(self.seq)] + ' 5']  # +1 offsets for alignment, +2 for codon size
            else:
                orfLength = len(self.seq) - start           # precalc orf length - section 6
                orf = [(previousFrame + 1) * -1, 1, len(self.seq) - start,
                       orfLength, self.seq[start:len(self.seq) - start]  + ' 6']
            if (orfLength > self.minGene):                  # minGene comparison, if smaller don't append
                orfList.append(orf)                         # add orf to the list
        return orfList

    def revComp(self):                                      # reverse the seq for use in bottom strand
        self.seq = self.seq.replace('T', 'a').replace('A', 't').replace('G', 'c').replace('C', 'g').upper()[::-1]

    def processBottomOrfs(self):
        self.revComp()                                      # reverse the sequence
        bottomOrfs = self.findORFs(0)                       # call findORFS and return orf list, 0 for bottom strand
        return bottomOrfs                                   # return the list

    def findStrandedOrfs(self):
        topStrandOrf = self.findORFs(1)                     # 1 for top strand
        # print('top strand ', topStrandOrf)                # debug code
        bottomStrandOrf = self.processBottomOrfs()          # call code for bottom strand
        # print('bottom strand ', bottomStrandOrf)          # debug code
        return topStrandOrf + bottomStrandOrf               # add strands and return

# it's always hard to find the one you want
class NucParams:
    '''
    Convert FastA sequences into codon, amino acid, and nucleotide dictionaries.
    Input is a FastA string
    Output depends on the method.
    methods: addSequence, aaComposition, nucComposition, codonComposition, and nucCount.
    '''
    rnaCodonTable = {
        # RNA codon table
        # U
        'UUU': 'F', 'UCU': 'S', 'UAU': 'Y', 'UGU': 'C',  # UxU
        'UUC': 'F', 'UCC': 'S', 'UAC': 'Y', 'UGC': 'C',  # UxC
        'UUA': 'L', 'UCA': 'S', 'UAA': '-', 'UGA': '-',  # UxA
        'UUG': 'L', 'UCG': 'S', 'UAG': '-', 'UGG': 'W',  # UxG
        # C
        'CUU': 'L', 'CCU': 'P', 'CAU': 'H', 'CGU': 'R',  # CxU
        'CUC': 'L', 'CCC': 'P', 'CAC': 'H', 'CGC': 'R',  # CxC
        'CUA': 'L', 'CCA': 'P', 'CAA': 'Q', 'CGA': 'R',  # CxA
        'CUG': 'L', 'CCG': 'P', 'CAG': 'Q', 'CGG': 'R',  # CxG
        # A
        'AUU': 'I', 'ACU': 'T', 'AAU': 'N', 'AGU': 'S',  # AxU
        'AUC': 'I', 'ACC': 'T', 'AAC': 'N', 'AGC': 'S',  # AxC
        'AUA': 'I', 'ACA': 'T', 'AAA': 'K', 'AGA': 'R',  # AxA
        'AUG': 'M', 'ACG': 'T', 'AAG': 'K', 'AGG': 'R',  # AxG
        # G
        'GUU': 'V', 'GCU': 'A', 'GAU': 'D', 'GGU': 'G',  # GxU
        'GUC': 'V', 'GCC': 'A', 'GAC': 'D', 'GGC': 'G',  # GxC
        'GUA': 'V', 'GCA': 'A', 'GAA': 'E', 'GGA': 'G',  # GxA
        'GUG': 'V', 'GCG': 'A', 'GAG': 'E', 'GGG': 'G'  # GxG
    }

    def __init__(self, inString=''):
        '''
        class init
        input: FastA string, takes a string and defaults to an empty string
        output: make dictionaries for amino actid, codons, and nucleotide, and then set them to zero
        '''
        self.aaComp = {aminoAcid: 0 for aminoAcid in NucParams.rnaCodonTable.values()}
        self.codonComp = {codon: 0 for codon in NucParams.rnaCodonTable}
        self.nucComp = {nucTide: 0 for nucTide in 'ACGTUN'}

    def addSequence(self, inSeq):
        '''
        addSequence
        Add a sequence to the dictionaries one codon at a time from a FastA string
        input: FastA sequence
        output: creates dictionaries for amino acid compostion, nucleotide composition,
        codon compostion
        '''
        # print(inSeq)
        mySeq = inSeq.replace(' ', '')              # remove spaces from input, and assign to a new name
        for index in range(0, len(mySeq), 3):       # that master for loop
            codon = mySeq[index:index+3].upper()    # set to upper case
                                                    # add the nucs
            for i in codon:                         # step through the codon
                if i in 'ACGTUN':                   # check if a valid nuc
                    self.nucComp[i] += 1            # add the nuc to the count
                                                    # add the codons and amino acids
            codon = codon.replace('T', 'U')         # convert to RNA for dictionary
            if codon in self.rnaCodonTable:         # check to see that the codon is in the dictionary
                self.aaComp[self.rnaCodonTable.get(codon)] += 1
                self.codonComp[codon] += 1

    def aaComposition(self):
        ''' aaComposition dictionary access, keys are amino acids'''
        return self.aaComp

    def nucComposition(self):
        '''nucComposition dictionary access, keys are nucleotides '''
        return self.nucComp

    def codonComposition(self):
        '''codonComposition dictionary access, keys are codons '''
        return self.codonComp

    def nucCount(self):
        '''nucCount - get the count of the nucleotides '''
        return sum(self.nucComp.values())

import sys

class FastAreader:
    '''
    Define objects to read FastA files.

    instantiation:
    thisReader = FastAreader ('testTiny.fa')
    usage:
    for head, seq in thisReader.readFasta():
        print (head,seq)
    '''

    def __init__(self, fname=''):
        '''contructor: saves attribute fname '''
        self.fname = fname

    def doOpen(self):
        ''' Handle file opens, allowing STDIN.'''
        if self.fname is '':
            return sys.stdin
        else:
            return open(self.fname)

    def readFasta(self):
        ''' Read an entire FastA record and return the sequence header/sequence'''
        header = ''
        sequence = ''
        with self.doOpen() as fileH:
            header = ''
            sequence = ''
            # skip to first fasta header
            line = fileH.readline()
            while not line.startswith('>'):
                line = fileH.readline()
            header = line[1:].rstrip()

            for line in fileH:
                if line.startswith('>'):
                    yield header, sequence
                    header = line[1:].rstrip()
                    sequence = ''
                else:
                    sequence += ''.join(line.rstrip().split()).upper()
        yield header, sequence

class ProteinParam:
    '''
    Calculates the physical-chemical properties of a protein sequence
    input: a protein string
         Example sequence: VLSPADKTNVKAAW
    output: Number of amino acids found in the string
            The calculated molecular weight of the protein
            The calculated molar extinction coefficient of the protein
            The calculated mass extinction coefficient of the protein
            The calculated theoretical pI of the protein
            And the amino acid percentage composition of the protein
    Public Methods: aaCount(), pI(), aaComposition(), molarExtinction(), massExtinction(),
            and molecularWeight(self)
    private methods: __init__(protein), _charge_(pH)
    '''

    # These tables are for calculating:
    #     molecular weight (aa2mw), along with the mol. weight of H2O (mwH2O)
    #     absorbance at 280 nm (aa2abs280)
    #     pKa of positively charged Amino Acids (aa2chargePos)
    #     pKa of negatively charged Amino acids (aa2chargeNeg)
    #     and the constants aaNterm and aaCterm for pKa of the respective termini
    #  Feel free to move these to appropriate methods as you like

    # As written, these are accessed as class attributes, for example:
    # ProteinParam.aa2mw['A'] or ProteinParam.mwH2O

    aa2mw = {
        'A': 89.093, 'G': 75.067, 'M': 149.211, 'S': 105.093, 'C': 121.158,
        'H': 155.155, 'N': 132.118, 'T': 119.119, 'D': 133.103, 'I': 131.173,
        'P': 115.131, 'V': 117.146, 'E': 147.129, 'K': 146.188, 'Q': 146.145,
        'W': 204.225, 'F': 165.189, 'L': 131.173, 'R': 174.201, 'Y': 181.189
    }

    mwH2O = 18.015
    aa2abs280 = {'Y': 1490, 'W': 5500, 'C': 125}

    aa2chargePos = {'K': 10.5, 'R': 12.4, 'H': 6}
    aa2chargeNeg = {'D': 3.86, 'E': 4.25, 'C': 8.33, 'Y': 10}
    aaNterm = 9.69
    aaCterm = 2.34

    def __init__(self, protein):                        # this is basically the constructor for ProteinParam
        '''
        Creates, and inits the amino acid count dictionary, and then analyzes the input sequence
        for the count of each amino acid, which is then added to the dictionary's amino acid count
        '''
        self.aminoAcidCount = {aa: 0 for aa in 'ACDEFGHIKLMNPQRSTVWY'}  # a dictionary for the count of amino acids
        self.myProtein = protein.upper().replace(' ', '')                        # make a public-ish object

        for aminoAcid in self.myProtein:
            if aminoAcid in self.aminoAcidCount:        # check to see that the amino acid is in the dictionary
                self.aminoAcidCount[aminoAcid] += 1     # for each amino in the protein increase its count

    def aaCount(self):
        '''
        Count the number of amino acids in the protein string, and record them in the dictionary
        '''
        totalCount = 0                                  # object to hold the total count
        for count in self.aminoAcidCount.values():      # traverse through the aminoAcidCount list
                totalCount += count                     # sum the count
        return totalCount                               # return the total count of amino acids

    def pI(self, precision = 2):                        # binary search extra credit version
        '''
        Calculate the theoretical isoelectric point, estimated by finding the particular pH that
        yields a neutral net charge that is close to 0.
        Input: optional precision value, defaulted to two.
        '''
        lopH = 0.0
        hipH = 14.0
        while (hipH - lopH) > (10**(-1 * precision)):   # setting precision of conditional check
            midpH = lopH + (hipH - lopH) / 2            # making a new mid point based on adjusted hi/low values
            if self._charge_(midpH) > 0:                # is the charge greater than zero at this pH?
                lopH = midpH                            # set the new low at the current mid
            else:
                hipH = midpH                            # if not then set the new hi at the current mid
        return lopH + (hipH - lopH) / 2                 # derive the pH to return

    def pIregular(self):                                # regular version
        '''
        Calculate the theoretical isoelectric point, estimated by finding the particular pH that
        yields a neutral net charge that is close to 0.
        '''
        bestCharge = 100000000                          # This is David's version from class
        for pH100 in range(0, 1400 + 1):                # offset to avoid fp errors
            pH = pH100 / 100                            # adjust offset
            thisCharge = abs(self._charge_(pH))         # call charge on current pH
            if thisCharge < bestCharge:                 # check this charge if higher then pH was found
                bestCharge = thisCharge                 # adjust charge
                bestpH = pH                             # adjust pH
        return bestpH

    def aaComposition(self):
        '''
        Return a keyed single letter amino acid dictionary of all 20 amino acids. This dictionary holds the count
        of each amino acid in the protein sequence.
        '''
        return self.aminoAcidCount                      # return the count of amino acids

    def _charge_(self, pH): # private method
        '''
        Calculates the net charge on the protein at a specific pH.
        input: pH value
        output: calculated total charge of the protein
        '''
        positiveCharge = 0.0                            # init vars
        negitiveCharge = 0.0
        for aminoAcid in self.aa2chargePos.keys():      # loop through
            positiveCharge += self.aminoAcidCount[aminoAcid] * (10 ** (self.aa2chargePos[aminoAcid]) /
                                                                (10**(self.aa2chargePos[aminoAcid]) + 10**pH))
        positiveCharge += (10 ** self.aaNterm) / (10 ** self.aaNterm + 10 ** pH)    # add in the Nterminus

        for aminoAcid in self.aa2chargeNeg.keys():
            negitiveCharge += self.aminoAcidCount[aminoAcid] * (10**pH / (10**(self.aa2chargeNeg[aminoAcid]) + 10**pH))
        negitiveCharge += 10**pH / (10 ** self.aaCterm + 10 ** pH)                  # add in the Cterminus

        return positiveCharge - negitiveCharge          # return the net charge

    def molarExtinction(self, Cystine = True):          # extra credit version
        '''
        Calculates the molar extinction of the protein string
        𝐸 = 𝑁y𝐸y + 𝑁w𝐸w + 𝑁c𝐸c
        Ny is the number of tyrosines, Nw tryptophans, Nc cysteines
        Ey, Ew, and Ec are the extinction coefficients for tyrosine, tryptophan, and cysteine
        '''
        moExt = 0.0                                     # init var
        for aminoAcid in self.aa2abs280:                # traverse through the aa2abs280 dictionary
            moExt += self.aminoAcidCount[aminoAcid] * self.aa2abs280[aminoAcid] # calculate molar extinction
        if not Cystine:                                 # if reducing remove the cystines
            moExt -= self.aminoAcidCount['C'] * self.aa2abs280['C']
        return moExt

    def molarExtinctionRegular(self):                   # Regular version
        '''
        Calculates the molar extinction of the protein string
        𝐸 = 𝑁y𝐸y + 𝑁w𝐸w + 𝑁c𝐸c
        Ny is the number of tyrosines, Nw tryptophans, Nc cysteines
        Ey, Ew, and Ec are the extinction coefficients for tyrosine, tryptophan, and cysteine
        '''
        moExt = 0.0                                     # init var
        for aminoAcid in self.aa2abs280:                # traverse through the aa2abs280 dictionary
            moExt += self.aminoAcidCount[aminoAcid] * self.aa2abs280[aminoAcid] # calculate molar extinction
        return moExt

    def massExtinction(self, Cystine = True):           # extra credit version
        '''
        Calculates the mass extinction
        By dividing the the molar extinction by molecular weight the mass extinction can be derived
        If the molecular weight is "false" this will return zero.
        '''
        myMW = self.molecularWeight()
        myMWsansCystine = self.molecularWeight() - (self.aminoAcidCount.get('C') * (self.aa2mw.get('C') - self.mwH2O))
        if Cystine:                                     # Use this return if oxidizing
            return self.molarExtinction() / myMW if myMW else 0.0
        else:                                           # Use this if reducing (remove the cystines)
            return self.molarExtinction() / myMWsansCystine if myMWsansCystine else 0.0

    def massExtinctionRegular(self):                     # Regular version
        '''
        Calculates the mass extinction
        By dividing the the molar extinction by molecular weight the mass extinction can be derived
        If the molecular weight is "false" this will return zero.
        '''
        myMW = self.molecularWeight()
        return self.molarExtinction() / myMW if myMW else 0.0

    def molecularWeight(self):
        '''
        Calculates the molecular weight of the protein string.
        This is almost just a sum of the molecular weights, but each time an amino acid is added
        to the chain a water molecule is released, sa that must be accounted for in the calculation
        '''
        totalWeight = 0.0                                   # init object to hold the total weight
        if self.aaCount() > 0:                              # check to see if there are amino acids
            totalWeight = self.mwH2O                        # add the water molecule if there is an amino string
        for aminoAcid in self.aminoAcidCount.keys():        # traverse through the aminoAcidCount dictionary
            totalWeight += self.aminoAcidCount.get(aminoAcid) * (self.aa2mw.get(aminoAcid) - self.mwH2O) # sum the weight
        return totalWeight
