#!/usr/bin/env python2
from __future__ import with_statement

import sys
import os
import itertools

# Returns true if significant differences (> tolerance) are found,
# returns false otherwise.
def files_differ(file1, file2, tolerance, verbose=True, keepGoing=False):
    """Find numerical differences of a star-discrepancy output file."""

    foundDifference = False
    assert tolerance >= 0, "Tolerance cannot be negative!"

    linect = 1
    with open(file1,'r') as inf1, open(file2, 'r') as inf2:

        for rawLine1,rawLine2 in itertools.izip(inf1,inf2):
        # Read each file line-by-line breaking up on spaces
            line1 = rawLine1.rstrip()
            line2 = rawLine2.rstrip()

            array1 = line1.split(" ")
            array2 = line2.split(" ")

        # Test that the lines have the same # of elements
            if len(array1) != len(array2):
                if verbose:
                    print("Test failed: on line " + str(linect) + ", different number of entries.")
                    print("in files " + file1 + " and " + file2)
                if not keepGoing:
                    return True
                else:
                    foundDifference = True

            itemct = 1
            # Test each element on the line for differences
            for a,b in itertools.izip(array1, array2):
                delta = abs(float(a) - float(b))
                if delta > tolerance:
                    if verbose:
                        print("Test failed on line " + str(linect) + ", element " + str(itemct))
                        print("Allowable error is "+str(tolerance)+ ", actual diff is " +str(delta))
                        print("in files " + file1 + " and " + file2)
                    if not keepGoing:
                        return True
                    else:
                        foundDifference = False
                itemct += 1

        linect += 1
    return foundDifference

def main(args):
    if len(args) < 2:
        print("Usage: %s <outdir> <refdir> [-s] [-k]"%args[0])
        print("")
        print("testdir: path of directory where the output to be tested lives")
        print("refdir: path of directory where reference output lives")
        print("-s : activates silent mode (do not print test results)")
        print("-k : keep going (try remaining tests even if one fails)")
        sys.exit(1)

    verbose = True
    keepGoing = False

    outdir = args[1]
    refdir = args[2]

    # Directory exists and is not empty
    if not os.path.isdir(outdir) or not os.listdir(outdir):
        print("You gave me " + outdir + " as outdir, but it is either empty or")
        print("not a directory. Please double-check your arguments and try again.")
        sys.exit(1)
    if not os.path.isdir(refdir) or not os.listdir(refdir):
        print("You gave me " + refdir + " as outdir, but it is either empty or")
        print("not a directory. Please double-check your arguments and try again.")
        sys.exit(1)

    testFileList = os.listdir(outdir)
    refFileList = os.listdir(refdir)

    # All files to be tested should have corresponding entries in the reference

    missingFiles = set(testFileList).difference(set(refFileList))
    if missingFiles:
        print("Some of the files in " + outdir + " are not in " + refdir)
        print(missingFiles)
        sys.exit(1)

    if "-s" in args or "--silent" in args:
        verbose = False
    if "-k" in args or "--keep-going" in args:
        keepGoing = True

    anyDifferent = False

    for testFile in testFileList:
        fname1 = os.path.join(outdir, testFile)
        fname2 = os.path.join(refdir, testFile)
#        print("Testing " + fname1 + " against " + fname2)
        differ = files_differ(fname1,fname2,1e-3,verbose,keepGoing)

        if differ and not keepGoing:
            sys.exit(1)
        elif differ:
            anyDifferent = True

    if anyDifferent:
        sys.exit(1)
    else:
        sys.exit(0)

if __name__ == "__main__":
    main(sys.argv)
