#!/usr/bin/env python3
"""
Author : Erick Samera
Date   : 2022-01-08 -> 2022-08-18
Purpose: Calculates the Tm of primers and estimates an appropriate annealing temp. for different polymerases.
"""

from argparse import (
    Namespace,
    ArgumentParser,
    ArgumentDefaultsHelpFormatter,
    RawDescriptionHelpFormatter)
from pathlib import Path
from sys import argv
from math import log
try:
    from pandas import DataFrame
except ModuleNotFoundError:
    print("WARNING!!! Install pandas or Erick will have to write the results manually :(\n")

# --------------------------------------------------
def get_args() -> Namespace:
    """ Get command-line arguments """

    parser = ArgumentParser(
        description='Program parses primers and returns an estimated annealing temperature (°C).',
        #usage='%(prog)s',
        epilog='v1.0.0 : We should put the versioning in the program.',
        formatter_class=ArgumentDefaultsHelpFormatter)

    subparsers = parser.add_subparsers(
        title='program run options',
        required=True,
        description="choose to either load a file or enter a single primer")

    # options for parsing ThermoFisher-formatted files
    # --------------------------------------------------
    file_parser = subparsers.add_parser(
        'file',
        description= \
            "Online ThermoFisher calculator does batch processing with a special file format. " \
            "See bottom for details.",
        help='load a formatted file containing primer information',
        epilog= \
            'acceptable file format:\n'
            'ID#1 Primer#1 ; ID#2 Primer#2 \\n\t\tThermoFisher website format (.txt)\n' \
            'ID#1   Primer#1    ID#2    Primer#2 \\n\t\tTab delimited format (.tab/tsv)\n' \
            'ID#1,Primer#1,ID#2,Primer#2\\n\t\t\tComma-separated values format (.csv)',
        formatter_class=RawDescriptionHelpFormatter)
    file_parser.add_argument(
        'input_path',
        type=Path,
        help='path of Thermo-formatted file (see bottom)')
    file_parser.add_argument(
        '-o',
        '--out',
        dest='out',
        metavar='PATH',
        type=Path,
        default=None,
        help='path of [.csv] to output results, otherwise only prints to console')
    file_parser.add_argument(
        '-p',
        '--pol',
        dest='pol',
        metavar='POL',
        type=str,
        choices=['SuperFi', 'Phusion', 'DreamTaq'],
        default='Phusion',
        help="Specify polymerase to use {SuperFi, Phusion, DreamTaq} (default: Phusion)")
    file_parser.add_argument(
        '-c',
        '--conc',
        dest='conc',
        metavar='CONC',
        type=float,
        default=0.5,
        help="Primer concentration (uM) (default: 0.5)")

    # options for parsing a single primer set in the command line
    # --------------------------------------------------
    primer_parser = subparsers.add_parser(
        'primer',
        description="Input details for single primer processing.",
        help='load a single pair of primers',
        formatter_class=ArgumentDefaultsHelpFormatter)
    primer_parser.add_argument(
        '-o',
        '--out',
        dest='out',
        metavar='PATH',
        type=Path,
        default=None,
        help='path of [.csv] to output results, otherwise only prints to console')
    primer_parser.add_argument(
        '-p',
        '--pol',
        dest='pol',
        metavar='POL',
        type=str,
        choices=['SuperFi', 'Phusion', 'DreamTaq'],
        default='Phusion',
        help="Specify polymerase to use {SuperFi, Phusion, DreamTaq}")
    primer_parser.add_argument(
        '-c',
        '--conc',
        dest='conc',
        metavar='CONC',
        type=float,
        default=0.5,
        help="Primer concentration (uM)")

    required_group = primer_parser.add_argument_group(
        title='primer arguments (required)',
        description="Note: Not case-sensitive, no ambiguous allowed.")
    required_group.add_argument(
        '-f',
        '--fwd',
        dest='f_seq',
        metavar='SEQ',
        type=str,
        default=None,
        required=True,
        help="(5'-3') forward primer sequence")
    required_group.add_argument(
        '-r',
        '--rev',
        dest='r_seq',
        metavar='SEQ',
        type=str,
        default=None,
        required=True,
        help="(5'-3') reverse primer sequence")

    # if no arguments are initially given, output the help and quit
    if len(argv) == 1:
        parser.print_help()
        parser.exit()

    args = parser.parse_args()

    # parser errors for file input
    # --------------------------------------------------

    # parser errors for primer input
    # --------------------------------------------------
    try:
        # primer sequences are sometimes separated into groups of 3, convert this
        args.f_seq = ''.join(args.f_seq.upper().split(' '))
        args.r_seq = ''.join(args.r_seq.upper().split(' '))

        # validate nucleotides
        allowed_nucs = 'ATCG'
        if not all([True if nuc in allowed_nucs else False for nuc in args.f_seq.upper()]) \
        or not all([True if nuc in allowed_nucs else False for nuc in args.r_seq.upper()]):
            parser.error('Only the standard nucleotides are allowed: {A, T, C, G}')
    except AttributeError:
        # The attributes, f_seq and r_seq, are only defined in args if the 'primer' subparser is used.
        # In the case that the file subparser is used, f_seq and r_seq will not be defined as part of args
        # and will throw an AttributeError which is exc
        pass

    return args
# --------------------------------------------------
class Primer:
    """
    A class to represent a single primer.

    Attributes (that you care about):
        name (str): the name of the primer, required by ThermoFisher format
        seq (str): the primer sequenece (5'-3')
        primer_conc (float): concentration (uM) of the primer
        salt_conc (float): concentration (uM) of the salts
    """

    def __init__(self, name: str, seq: str, primer_conc: float) -> None:
        self.name = name
        self.seq = seq.upper()
        self.primer_conc = primer_conc * 1e-6
        self.salt_conc = 50 * 1e-3

    def _calculate_thermodynamics(self, values_arg: dict, inits_arg: dict) -> dict:
        """
        Function performs thermodynamic array calculations for the dS and dH of a given primer set.

        Paramters:
            values_arg (dict): source values to process
            inits_arg (dict): values to intialize prior to calculation

        Returns:
            (dict):
                dS (float): entropy
                dH (float): enthalpy
        """

        total_dH = inits_arg['dH_ini']
        total_dS = inits_arg['dS_ini']

        for i, s in enumerate(self.seq):
            if i < len(self.seq)-1:
                total_dH += values_arg[s][self.seq[i+1]]['dH']
                total_dS += values_arg[s][self.seq[i+1]]['dS']

        return {'dS': total_dS, 'dH': total_dH}

    def _Tm_All97(self) -> float:
        """
        Function to get the Tm using the Allawi-SantaLucia (1997) method.

        Parameters:
            None

        Returns:
            (float): primer Tm (°C)
        """

        Tm_All97_data = {
            'A': {
                'A':{'dS':-22.2,'dH':-7900.0},
                'C':{'dS':-22.4,'dH':-8400.0},
                'T':{'dS':-20.4,'dH':-7200.0},
                'G':{'dS':-21.0,'dH':-7800.0}},
            'C': {
                'A':{'dS':-22.7,'dH':-8500.0},
                'C':{'dS':-19.9,'dH':-8000.0},
                'T':{'dS':-21.0,'dH':-7800.0},
                'G':{'dS':-27.2,'dH':-10600.0}},
            'T': {
                'A':{'dS':-21.3,'dH':-7200.0},
                'C':{'dS':-22.2,'dH':-8200.0},
                'T':{'dS':-22.2,'dH':-7900.0},
                'G':{'dS':-22.7,'dH':-8500.0}},
            'G': {
                'A':{'dS':-22.2,'dH':-8200.0},
                'C':{'dS':-24.4,'dH':-9800.0},
                'T':{'dS':-22.4,'dH':-8400.0},
                'G':{'dS':-19.9,'dH':-8000.0}}}

        dS_ini = 0
        dH_ini = 0

        if self.seq.endswith(('A', 'T')):
            dS_ini += 4.1
            dH_ini += 2300
        if self.seq.endswith(('C', 'G')):
            dS_ini -= 2.8
            dH_ini += 100
        if self.seq.startswith(('A', 'T')):
            dS_ini += 4.1
            dH_ini += 2300
        if self.seq.startswith(('C', 'G')):
            dS_ini -= 2.8
            dH_ini += 100

        thermodynamic_data = self._calculate_thermodynamics(
            values_arg=Tm_All97_data,
            inits_arg={'dH_ini': dH_ini, 'dS_ini': dS_ini})
        res = (thermodynamic_data['dH'] / (1.9872 * log(self.primer_conc / 4.0) + thermodynamic_data['dS'])) \
            + (16.6 * log(0.215273974689348) / log(10)) \
            - 273.15

        res_adj = (res+3)*0.9376798568+4.5185404499

        if res_adj < 0:
            res_adj = 0
        elif res_adj > 95:
            res_adj =  95

        return res_adj

    def _Tm_taq(self) -> float:
        """
        Function to get the Tm using the SantaLucia (1996) method.

        Parameters:
            None

        Returns:
            (float): primer Tm (°C)
        """

        Tm_San96_data = {
            'A': {
                'A':{'dS':-23.6,'dH':-8400.0},
                'C':{'dS':-23.0,'dH':-8600.0},
                'T':{'dS':-18.8,'dH':-6500.0},
                'G':{'dS':-16.1,'dH':-6100.0}},
            'C': {
                'A':{'dS':-19.3,'dH':-7400.0},
                'C':{'dS':-15.6,'dH':-6700.0},
                'T':{'dS':-16.1,'dH':-6100.0},
                'G':{'dS':-25.5,'dH':-10100.0}},
            'T': {
                'A':{'dS':-18.5,'dH':-6300.0},
                'C':{'dS':-20.3,'dH':-7700.0},
                'T':{'dS':-23.6,'dH':-8400.0},
                'G':{'dS':-19.3,'dH':-7400.0}},
        'G': {
                'A':{'dS':-20.3,'dH':-7700.0},
                'C':{'dS':-28.4,'dH':-11100.0},
                'T':{'dS':-23.0,'dH':-8600.0},
                'G':{'dS':-15.6,'dH':-6700.0}}}

        thermodynamic_data = self._calculate_thermodynamics(
            values_arg=Tm_San96_data,
            inits_arg={'dH_ini': 0.0, 'dS_ini': -0.0})
        NaEquiv = 0.15527397
        nucleotide_F_term  = -15.894952
        entropy_correction =  0.368 * (len(self.seq) - 1.0) * log(NaEquiv)
        thermodynamic_data['dS'] += + entropy_correction

        res = thermodynamic_data['dH'] \
        / (thermodynamic_data['dS'] + 1.9872 * nucleotide_F_term) \
        - 273.15

        if res < 0:
            res=0
        elif res > 95:
            res = 95

        if len(self.seq) < 21:
            res_ajusted = res * 0.9085395477132917 - 3.707388372789194
        else:
            res_ajusted = (res + 3) * 0.9085395477132917 - 3.707388372789194

        return res_ajusted

    def return_Tm(self, pol_arg: str) -> float:
        """
        Function calculats melting temperature (Tm) based on polymerase.

        Parameters:
            pol_arg (str): polymerase name from defined choices

        Returns:
            (float): primer melting temperature (°C)
        """

        if pol_arg in ['SuperFi', 'Phusion']:
            Tm = self._Tm_All97()
        elif pol_arg in ['DreamTaq']:
            Tm = self._Tm_taq()
        else:
            Tm = self._Tm_taq()

        return Tm
# --------------------------------------------------
def main() -> None:
    """ Insert docstring here """

    args = get_args()

    # create a list of primers to process
    primers_list = []
    try:
        if args.f_seq and args.r_seq:
            fwd_primer = Primer(name='Forward', seq=args.f_seq, primer_conc=args.conc)
            rev_primer = Primer(name='Reverse', seq=args.r_seq, primer_conc=args.conc)
            primers_list.append({'fwd_primer': fwd_primer, 'rev_primer': rev_primer})
    except AttributeError:
        # The attributes, f_seq and r_seq, are only defined in args if the 'primer' subparser is used.
        # In the case that the file subparser is used, f_seq and r_seq will not be defined as part of args
        # and will throw an AttributeError which is excepted here.
        pass
    try:
        if args.input_path:
            with open(args.input_path.resolve(), 'r', encoding='UTF8') as input_file:
                for line in input_file.readlines():

                    # process as Thermo-formatted file
                    # --------------------------------------------------
                    if len([info.strip() for info in line.split(';')]) == 2:
                        line_info = [info.split(' ') for info in [info.strip() for info in line.split(';')]]
                        fwd_primer = Primer(name=line_info[0][0], seq=line_info[0][1], primer_conc=args.conc)
                        rev_primer = Primer(name=line_info[1][0], seq=line_info[1][1], primer_conc=args.conc)

                    # process as tab-delimited file/tab-separated values file
                    # --------------------------------------------------
                    elif len([info.strip() for info in line.split('\t')]) == 4:
                        line_info = [info.strip() for info in line.split('\t')]
                        fwd_primer = Primer(name=line_info[0], seq=line_info[1], primer_conc=args.conc)
                        rev_primer = Primer(name=line_info[2], seq=line_info[3], primer_conc=args.conc)

                    # process as comma-seperated values file
                    # --------------------------------------------------
                    elif len([info.strip() for info in line.split(',')]) == 4 or len([info.strip() for info in line.split(',')]) == 5:
                        line_info = [info.strip() for info in line.split(',')]
                        fwd_primer = Primer(name=line_info[0], seq=line_info[1], primer_conc=args.conc)
                        rev_primer = Primer(name=line_info[2], seq=line_info[3], primer_conc=args.conc)

                    # process as comma-separated values output-file
                    # --------------------------------------------------
                    elif len([info.strip() for info in line.split(',')]) == 9:
                        line_info = [info.strip() for info in line.split(',')]
                        if line_info[0]:
                            fwd_primer = Primer(name=line_info[1], seq=line_info[2], primer_conc=args.conc)
                            rev_primer = Primer(name=line_info[4], seq=line_info[5], primer_conc=args.conc)
                        else:
                            # In the output-file, the header line contains an empty string in the first position.
                            # This 'else' conditional will pass over the header file and only process the lines thereafter.
                            pass 

                    # add the processed primers to the processing list
                    # --------------------------------------------------
                    try:
                        primers_list.append({'fwd_primer': fwd_primer, 'rev_primer': rev_primer})
                    except UnboundLocalError:
                        # The forward and reverse primer will be defined if they fall into the formats above.
                        # If the line does not follow any formats as written above, the variables, fwd_primer and rev_primer
                        # won't be defined and this will throw an UnboundLocalError.
                        pass
    except AttributeError:
        # The attributes, f_seq and r_seq, are only defined in args if the 'primer' subparser is used.
        # In the case that the file subparser is used, f_seq and r_seq will not be defined as part of args
        # and will throw an AttributeError which is excepted here.
        pass

    # create a list of output to print out to terminal
    output_list = []
    for primer_pair in primers_list:
        primer_pair_Ta = calculate_Ta(
            primer_1_arg=primer_pair['fwd_primer'],
            primer_2_arg=primer_pair['rev_primer'],
            pol_arg=args.pol)
        output_list.append({
            'FWD Primer': primer_pair['fwd_primer'].name,
            'FWD Primer Seq': primer_pair['fwd_primer'].seq,
            'Tm 1 (*C)': primer_pair['fwd_primer'].return_Tm(args.pol),
            'REV Primer': primer_pair['rev_primer'].name,
            'REV Primer Seq': primer_pair['rev_primer'].seq,
            'Tm 2 (*C)': primer_pair['rev_primer'].return_Tm(args.pol),
            f'{args.pol} Ta (*C)': primer_pair_Ta['Ta'],
            'Notes': primer_pair_Ta['note'],})

    if output_list:
        output_DataFrame = DataFrame(output_list)
        print('\n', output_DataFrame)
        if args.out:
            print(f"Output to: {args.out.resolve()}")
            output_DataFrame.to_csv(args.out.resolve())
    else:
        print('Nothing!')

    return None

def calculate_Ta(primer_1_arg: Primer, primer_2_arg: Primer, pol_arg: str) -> dict:
    """
    Function will calculate annealing temperature (Ta) based on polymerase.

    Parameters:
        primer_1_arg (Primer): object representing forward primer
        primer_2_arg (Primer): object representing reverse primer
        pol_arg (str): polymerase name from defined choices

    Returns:
        (dict)
            Ta (float): calculated Ta using specified primers and polymerase
            note (str): warnings and notes from the Ta calculation
    """

    if pol_arg in ['SuperFi', 'Phusion']:
        if min(len(primer_1_arg.seq), len(primer_2_arg.seq)) < 21 \
            and min(primer_1_arg.return_Tm(pol_arg), primer_2_arg.return_Tm(pol_arg)) < 72:
            Ta = min(primer_1_arg.return_Tm(pol_arg), primer_2_arg.return_Tm(pol_arg))
        else:
            Ta = min(min(primer_1_arg.return_Tm(pol_arg), primer_2_arg.return_Tm(pol_arg)), 72)
    elif pol_arg in ['DreamTaq']:
        if min(primer_1_arg.return_Tm(pol_arg), primer_2_arg.return_Tm(pol_arg)) < 72:
            Ta = min((max(primer_1_arg.return_Tm(pol_arg), primer_2_arg.return_Tm(pol_arg))), 72)
        else:
            Ta = 72
    else:
        if min(primer_1_arg.return_Tm(pol_arg), primer_2_arg.return_Tm(pol_arg)) < 72:
            Ta = min((max(primer_1_arg.return_Tm(pol_arg), primer_2_arg.return_Tm(pol_arg))), 72)
        else:
            Ta = 72

    # notes to output in case of errors
    # --------------------------------------------------
    note = ''

    if abs(primer_1_arg.return_Tm(pol_arg) - primer_2_arg.return_Tm(pol_arg)) >= 5:
        note += "Tm difference of more than 5C or greater is not recommended. "
    if Ta < 45:
        note += "Annealing temperature lower than 45C is not recommended. "
    if (len(primer_1_arg.seq) < 7) or (len(primer_2_arg.seq) < 7):
        note += "Both primers need to be longer than 7 nt. "
    if (69 < primer_1_arg.return_Tm(pol_arg) < 72) and (69 < primer_2_arg.return_Tm(pol_arg) < 72):
        note += "A 2-step protocol (combined annealing/extension) is recommended when primer Tm values are higher than 69C, using 72C for annealing step. "
    if Ta >= 72:
        note += "Annealing temperature should not exceed 72C. "

    return {'Ta': Ta, 'note': note}
# --------------------------------------------------
if __name__ == '__main__':
    main()
