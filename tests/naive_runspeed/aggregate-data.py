#!/usr/bin/env python3

## Reads and aggregates the data from the speed tests.

import sys
import os
import re
import statistics

knlDir = 'knl-tests'
sbDir = 'sb-tests'

def extract_files(dirName):
    """Extracts the relevant log files from the given directory path.

    Return format: a dict of filename-to-metadata mappings. Filenames are
    strings, metadata are tuples of (arch, Npts, dim, M)"""

    output = {}

    if dirName == knlDir:
        arch = 'KN'
    else:
        arch = 'SB'

    for (dirpath, dirnames, filenames) in os.walk(dirName):
        for fname in filenames:
            if fname[-4:] == ".log":
                parentdir = dirpath.split("/")[-1]
                components = parentdir.split(".")
                nPtString  = components[-2]
                dimString  = components[-1]
                MString    = fname.split("-")[0]

                nPts = int(nPtString[1:])
                dim  = int(dimString[1:])
                M    = int(MString)

                FQFname = os.path.join(dirpath,fname) # fully qualified fname

                output[FQFname] = (arch, nPts, dim, M)

    return output

def get_md_dict():
    """Get metadata dictionary for all existing log files."""

    knlEntries = extract_files(knlDir)
    sbEntries = extract_files(sbDir)
    return {**knlEntries, **sbEntries}

def gen_times(fname):
    """Look in file to get a set of timings.

    Returns a tuple of timings: (mutation, fitness), where either entry
    is a list of all mutation/fitness times discovered in fname."""

    output = ([],[])

    mutationRegex = re.compile("mutations: [0-9]+\.[0-9]+")
    fitnessRegex = re.compile("fitness: [0-9]+\.[0-9]+")
    numericRegex = re.compile("[0-9]+.[0-9]+")

    with open(fname,'r') as infile:
        for line in infile.readlines():
            fitnessMatch = fitnessRegex.search(line)
            mutationMatch = mutationRegex.search(line)
            if mutationMatch:
                timeString = numericRegex.search(mutationMatch.group(0))
                mutationTime = float(timeString.group(0))
                output[0].append(mutationTime)
            if fitnessMatch:
                timeString = numericRegex.search(fitnessMatch.group(0))
                fitnessTime = float(timeString.group(0))
                output[1].append(fitnessTime)

    return output

def stats(info):
    """Gets statistical info for a tuple of ([mutationtime],[fitnesstime])"""

    (mutation,fitness) = info

    mutationStats = {}
    fitnessStats = {}

    mutationStats["mean"] =  statistics.mean(mutation)
    mutationStats["median"]  = statistics.median(mutation)
    mutationStats["stdev"]  = statistics.stdev(mutation)

    fitnessStats["mean"] =  statistics.mean(fitness)
    fitnessStats["median"]  = statistics.median(fitness)
    fitnessStats["stdev"]  = statistics.stdev(fitness)

    return(mutationStats,fitnessStats)

def pretty_print_info(tup):
    """Pretty print the 4-tuple of experiment identifiers"""
    (arch, npts, dim, M) = tup

    print(arch + " with " + str(npts) + " points in " + str(dim) + " dimensions and M=" + str(M))

def pretty_print_stats(pair):
    """Pretty print the statistical results from stats()"""
    (mut,fit) = pair

    print("        Mutation \t Fitness")
    print("Mean:   " + str(round(mut['mean'],3)) + " \t\t " + str(round(fit['mean'],3)))
    print("Median: " + str(round(mut['median'],3)) + " \t\t " + str(round(fit['median'],3)))
    print("Stdev:  " + str(round(mut['stdev'],3)) + " \t\t " + str(round(fit['stdev'],3)))


def main():
    outputs = {}
    results = get_md_dict()
    for fn in results:
        outputs[results[fn]] = stats(gen_times(fn))

    # Print SB info
    for info in outputs:
        if info[0] == "SB":
            pretty_print_info(info)
            pretty_print_stats(outputs[info])
            print("\n\n")

    # Print KNL info
    for info in outputs:
        if info[0] == "KN":
            pretty_print_info(info)
            pretty_print_stats(outputs[info])
            print("\n\n")



if __name__ == "__main__":
    main()
