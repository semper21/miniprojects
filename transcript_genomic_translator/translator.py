"""
Created on Feb 14, 2021

@author: ywkim
"""

import sys
import re
import numpy as np
import pandas as pd


def transcript_to_genome(cigar_: str, start_pos: int) -> dict:
    """Uses CIGAR string to map transcript positions to genomic positions, which is then stored in a dictionary

    :param cigar_: CIGAR string
    :param start_pos: 0-based starting position
    :return: dictionary {transcript_position:genomic_position}
    """

    dict_ = {}
    l = re.findall(r'\d+|\D+', cigar_)

    diff = start_pos
    t_pos = 0
    for i, x in zip(l[::2], l[1::2]):   # looping through two elements at a time (ex. (8,M), (7,D)...)
        # Here, we're assuming we're only using the original CIGAR operators (D, I, M)
        if x == 'D':
            diff += int(i)
        elif x == 'I':
            diff -= int(i)
            for p in range(int(i)):
                dict_[t_pos] = t_pos + diff + (int(i)-p)
                t_pos += 1
        elif x == 'M':
            for p in range(int(i)):
                dict_[t_pos] = t_pos + diff
                t_pos += 1
        else:
            raise ValueError("Unexpected CIGAR operator")

    return dict_


def check_for_nan(df: pd.DataFrame) -> list:
    """Checks dataframe for NaNs and returns transcripts that have missing values

    :param df: dataframe to check for NaNs
    :return: list of transcripts that have missing values
    """
    df_bool = df.isna().any(axis=1)
    idx = list(np.where(df_bool == True))[0]

    return df.iloc[idx]['transcript'].values.tolist()


def main(input_file_1: str, input_file_2: str, file_name: str):
    """Stores transcript to genome positions in a nested dictionary,
       from which genomic positions of query transcript positions are retrieved.

    :param input_file_1: A four column (tab-separated) file containing the transcripts
    :param input_file_2: A two column (tab-separated) file indicating a set of queries
    :param file_name: Name of the output file
    :return: -
    """

    # below lines assumes input file has no headers
    df_mapping = pd.read_csv(input_file_1, sep='\t', header=None,
                             names=['transcript', 'chromosome', 'pos', 'cigar'])
    # check for NaN so that rows with NaN values would be automatically deleted
    nan_tsc = list(set(check_for_nan(df_mapping)))
    for tsc in nan_tsc:
        print("missing transcript: " + tsc + " cannot be mapped")

    df_mapping = df_mapping.dropna()

    # if the query file contains NaN, ValueError will be raised to inform the user
    try:
        df_query = pd.read_csv(input_file_2, sep='\t', header=None,
                               names=['transcript', 'query_position'], dtype={'query_position': np.int64})
    except ValueError:
        raise ValueError("Query file contains NaN, please check your input")

    # turn values in df_mapping into lists
    transcripts = df_mapping['transcript'].values.tolist()
    positions = df_mapping['pos'].astype('int64').values.tolist()
    cigars = df_mapping['cigar'].values.tolist()

    # dictionary to map transcript to chromosome
    chr_dict = dict(zip(df_mapping.transcript, df_mapping.chromosome))

    translator = {}  # nested dictionary for each transcript
    # Making dictionaries of {[transcript coordinate]:genomic coordinates} for each transcript
    for idx, transcript in enumerate(transcripts):
        cigar = cigars[idx]
        start = positions[idx]
        translator[transcript] = transcript_to_genome(cigar, start)

    # Retrieve genomic coordinates for each queries
    chromosomes = []
    mapped_pos = []
    queries = df_query.values
    for query in queries:
        q_tsc = query[0]    # query transcript
        q_pos = query[1]    # query position
        try:
            chromosomes.append(chr_dict[q_tsc])
            try:
                mapped_pos.append(translator[q_tsc][q_pos])
            except KeyError:    # query position that has not been mapped - will be interpreted as being out of range
                print(q_pos + " : query position out of range and will be ignored")
                mapped_pos.append(np.nan)
        except KeyError:                # query transcript that has not been mapped either because it doesn't exist
            chromosomes.append(np.nan)  # in the input file with CIGAR strings or because it contained NaN values
            mapped_pos.append(np.nan)
            if q_tsc in nan_tsc:
                pass
            else:
                print(q_tsc + " cannot be mapped and will be ignored")

    df_query['chromosome'] = chromosomes
    df_query['mapped_position'] = mapped_pos

    df_query = df_query.dropna().astype({'mapped_position': np.int64})
    df_query.to_csv(file_name + '.tsv', sep='\t', header=False, index=False)


if __name__ == '__main__':
    mapping_file = sys.argv[1]  # A four column (tab-separated) file containing the transcripts
    query_file = sys.argv[2]  # A two column (tab-separated) file indicating a set of queries

    main(mapping_file, query_file, 'output')
