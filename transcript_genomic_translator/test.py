"""
Created on Feb 15, 2021

@author: ywkim
"""

import os
import argparse
import numpy as np
import pandas as pd
import translator


def parse_arguments():
    """Parses arguments

    :return: arguments
    """

    parser = argparse.ArgumentParser(description="arguments for testing")
    parser.add_argument('-d', action='store_true', help='output files will be deleted')

    return parser.parse_args()


class TestTranslator:
    def set_up(self):
        folder = os.path.dirname(os.path.abspath(__file__)) + '/'
        # create files
        df_1 = pd.DataFrame({'tr': ['TR1', 'TR2'],
                             'chr': ['CHR1', 'CHR2'],
                             'pos': [3, 10],
                             'cigar': ['8M7D6M2I2M11D7M', '20M']})
        df_1.to_csv(folder + 'file1.tsv', sep='\t', header=False, index=False)

        df_2 = pd.DataFrame({'tr': ['TR1', 'TR2', 'TR1', 'TR2'],
                             'pos': ['4', 0, 13, 10]})  # testing different dtypes
        df_2.to_csv(folder + 'file2.tsv', sep='\t', header=False, index=False)

        # create bad input mapping file
        df_bad_1 = pd.DataFrame({'tr': ['TR1', 'TR2'],
                                 'chr': ['CHR1', 'CHR2'],
                                 'pos': [3, 10],
                                 'cigar': ['8M7D6M2I2M11D7M', np.nan]})  # missing a CIGAR string
        df_bad_1.to_csv(folder + 'bad_file1.tsv', sep='\t', header=False, index=False)

    def test_translator(self):
        answer = {0: 3, 1: 4, 2: 5, 3: 6, 4: 7, 5: 8, 6: 9, 7: 10, 8: 18, 9: 19, 10: 20, 11: 21, 12: 22, 13: 23, 14: 24,
                  15: 24, 16: 24, 17: 25, 18: 37, 19: 38, 20: 39, 21: 40, 22: 41, 23: 42, 24: 43}
        assert answer == translator.transcript_to_genome('8M7D6M2I2M11D7M', 3), "The translator returns wrong values"

        print("test 1 done")

    def test_main(self):
        folder = os.path.dirname(os.path.abspath(__file__)) + '/'
        translator.main(folder + 'file1.tsv', folder + 'file2.tsv', 'output_test')  # outputs output_test.tsv

        df_test = pd.read_csv('output_test.tsv', sep='\t', header=None,
                              names=['transcript', 'query_position', 'chromosome', 'mapped_position'])

        chromosomes = df_test['chromosome'].values.tolist()
        positions = df_test['mapped_position'].values.tolist()

        chr_answer = ['CHR1', 'CHR2', 'CHR1', 'CHR2']
        pos_answer = [7, 10, 23, 20]

        assert chr_answer == chromosomes, "The translator returns wrong chromosomes"
        assert pos_answer == positions, "The translator returns wrong positions"

        print("test 2 done")

    def test_main_bad_input(self):
        folder = os.path.dirname(os.path.abspath(__file__)) + '/'
        translator.main(folder + 'bad_file1.tsv', folder + 'file2.tsv', 'output_bad_file')

        df_test = pd.read_csv('output_bad_file.tsv', sep='\t', header=None,
                              names=['transcript', 'query_position', 'chromosome', 'mapped_position'])

        chromosomes = df_test['chromosome'].values.tolist()
        positions = df_test['mapped_position'].values.tolist()

        chr_answer = ['CHR1', 'CHR1']
        pos_answer = [7, 23]

        assert chr_answer == chromosomes, "The translator returns wrong chromosomes"
        assert pos_answer == positions, "The translator returns wrong positions"

        print("test 3 for bad input done")

    def tear_down(self):
        os.remove('file1.tsv')
        os.remove('file2.tsv')
        os.remove('bad_file1.tsv')
        os.remove('output_test.tsv')
        os.remove('output_bad_file.tsv')


if __name__ == '__main__':

    TestTranslator().set_up()
    TestTranslator().test_translator()
    TestTranslator().test_main()
    TestTranslator().test_main_bad_input()

    args = parse_arguments()

    if args.d:
        TestTranslator().tear_down()    # this line is to delete test output files
