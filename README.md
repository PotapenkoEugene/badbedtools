# badbedtools
Bad bedtools copy

# Tests
Tests are performed on large files to exclude errors, so the test run takes 2 minutes. All the results of the program are compared with the results of the corresponding bedtools functions.

# Functionality
This tool disgustingly copies bedtools functions such as:
    
    -sort
    
    -merge
    
    -subtract
    
    -subtract with A options
    
    -intersect
    
    -getfasta

# Bad description of functions:

import badtools as bad


## sort() - Function for sorting <bed/gff/vcf>


Basic usage (for bed files):

**bad.sort(file)**

  
    Order the intervals in a file.
  
        :param bed_vcf_gff_file: File in <bed/gff/vcf> format  
    
        :param extension: <bed/gff/vcf>
  
        :return: string
    
## merge() - Function to merge overlapped intervals

Basic usage (for bed files):

**bad.merge(file)**

or if you want to merge almost overlapping intervals

**bad.merge(file,max_distance=100)**

    Merges overlapping sorted BED/GFF/VCF entries into a single interval.
    
    :param bed_vcf_gff_file: File in <bed/gff/vcf> format.
    
    :param extension: <bed/gff/vcf>.
    
    :param max_distance: Maximum distance between features allowed for features to be merged. Default is 0
    
    :param sorting: Add sorting step. default = True
    
    :return: string


## subtract_A - Function for delete intervals in file1 that overlaps with intervals in file2


Basic usage (for bed files):

**bad.subtract_A(file1, file2)**

    Delete intervals that is overlapped by another feature(s)
    
    :param bed_vcf_gff_file: File in <bed/gff/vcf> format.
    
    :param bed_vcf_gff_subtracted_file: Subracted file in <bed/gff/vcf> format.
    
    :param extension_1: <bed/gff/vcf>
    
    :param extension_2: <bed/gff/vcf>
    
    :param sorting: Add sorting step. Default = True
    
    :param merging: Add merging step. Default = True
    
    :return: string
        
## subtract() - Function for cutting some intervals that overlapped by another intervals


Basic usage (for bed files):

**bad.subtract(file1, file2)**


    Removes the portion(s) of an interval that is overlapped by another feature(s)
    
    :param bed_vcf_gff_file: File in <bed/gff/vcf> format.
    
    :param bed_vcf_gff_subtracted_file: Subracted file in <bed/gff/vcf> format.
    
    :param extension_1: <bed/gff/vcf>
    
    :param extension_2: <bed/gff/vcf>
    
    :param sorting: Add sorting step. Default = True
    
    :param merging: Add merging step. Default = True
    
    :return: string

## intersect() - Function to return overlaps beetween two intervals files


Basic usage (for bed files):

**bad.intersect(file1, file2)**


    Report overlaps between two feature files.
    
    :param bed_vcf_gff_file: File in <bed/gff/vcf> format
    
    :param bed_vcf_gff_intersect_file: File in <bed/gff/vcf> format
    
    :param extension: <bed/gff/vcf>
    
    :param extension_int: <bed/gff/vcf>
    
    :param subtract: Return non intersected regions
    
    :param sorting: Add sorting step. default = True
    
    :param merging: Add merging step. default = True
    
    :return: string

## getfasta() - extract DNA seq from fasta

    Extract DNA sequences from a fasta file based on feature coordinates.
    
    :param bed_vcf_gff: File in bed_vcf_gff formats
    
    :param fastafile: Fasta file !!! without spaces in sequence !!!
    
    :param outname: Name for output file.
    
    :param extension: Extension for interval file. Default = 'bed'
    
    :return: None (write file in the specified directory)

## getnospacefasta() - If you haven't fasta without space, use it

    If you haven't fasta without space, use it
    
    :param fastafile: File in fasta format
    
    :param outname: Output dir and name
    
    :return: None (write file in the specified directory)


## Speed up your work
By default, the functions perform all preprocessing steps with your files.
But as you can see in functions intersect(), subtract() and subtract_A() you can disable sorting and merging steps if you already have
sorted and merged files. And, accordingly, in merge() function you can disable sorting step.
