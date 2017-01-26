#!/usr/bin/env python2
from __future__ import with_statement

import os
import sys
import re

def files_differ(file1, file2, tolerance, verbose=True, keepGoing = False):
    """Determine if stdout files differ"""

    f1lines = []
    f2lines = []

    with open(file1,'r') as inf1:
        for line in inf1.readlines():
            if "Star discrepancy after" in line:
                f1lines.append(line)

    with open(file2,'r') as inf2:
        for line in inf2.readlines():
            if "Star discrepancy after" in line:
                f2lines.append(line)

    f1Final = f1lines[-1].split()
    f2Final = f2lines[-1].split()

    try:
        f1Gen = int(f1Final[3])
        f2Gen = int(f2Final[3])
        f1Discrepancy = float(f1Final[-1])
        f2Discrepancy = float(f2Final[-1])

        if f1Gen == f2Gen and abs(f1Discrepancy - f2Discrepancy) < tolerance:
            return False
        else:
            print("Failed stdout: values differed by more than tolerance.")
            print("For file " + file1)
            return True

    except ValueError:
        if verbose:
            print("Failed stdout: Stdout does not have right format!")
            print("For file " + file1)
        return True

    return False

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
        print("[WARN]: Some of the files in " + outdir + " are not in " + refdir)
        print(missingFiles)

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
