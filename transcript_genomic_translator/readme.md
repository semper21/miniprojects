
## Description
This script is for translating transcript coordinates to genomic coordinates.

### Environment
* Python 3.6
* Pytest 6.2.2
* Pandas 1.0.1+
* Numpy 1.18.1+

### Usage
Required arguments for the translator.py script are as follows:

|**argument**|**description**|
|------------|---------------|
|input_file_1|A four column tab-separated file containing the transcripts|
|input_file_2|A two column tab-separated file indicating a set of queries|

Usage is as follows:
```shell script
python translation.py input_file_1 input_file_2
```

Output: output.tsv in the current working directory

#### Testing
I have tested this script to make sure it maps transcript coordinates to genomic coordinates correctly and can handle missing values in input files.

Optional arguments are as follows:

|**argument**|**description**|
|------------|---------------|
|-d |pass this flag if deleting output files after testing is desired|

Usage is as follows:
```shell script
python test.py -d
```
or
```shell script
pytest test.py 
```
Note: -d argument should not be used for pytest

### Any Other Important Information
Assumptions:
1. Input files are tab-separated
2. Input files have no headers
3. Query file will not contain NaNs, and if it does, they will ignored
4. CIGAR string used here is the original version ([Original Exonerate CIGAR from 2001](https://github.com/NBISweden/GAAS/blob/master/annotation/knowledge/cigar.md#original-exonerate-cigar-from-2001))

Strengths:
- Mapping transcript coordinates to genomic coordinates happen once in the beginning, not every time for each query position
- Mapped coordinates are stored in a dictionary - use of dictionary allows quick look up 

Weaknesses:
- As the number of transcripts to map increases, utilizing dictionary may take up a lot of space (though this shouldn't be a huge problem  since memory isn't very expensive)
- Where there's missing values in the input files, they are handled as NaN, which means int will be converted to float

