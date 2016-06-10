import fileinput
import sys
import os

def replaceSearchKeys(file, searchExp, replaceExp):
    for line in fileinput.input(file, inplace=1):
        if searchExp in line:
            line = line.replace(searchExp,replaceExp)
        sys.stdout.write(line)

def printFile(cmakeFile):
    with open(cmakeFile) as f:
        for line in f:
            print(line)

if __name__ == '__main__':

    # First argument is the executable name
    # Second argument is the CPP source name

    ## Modifying the new CMakeLists.txtfor the new project
    cmakeFile = "CMakeLists.txt"

    # Executable filename
    execName = sys.argv[1]

    # Source cpp filename
    sourceName = sys.argv[2]

    # Search keys
    searchKeyExecutable = "<EXEC_NAME>"
    searchKeySource = "<EXEC_SOURCES>"

    # Modify the text file
    replaceSearchKeys(cmakeFile,searchKeyExecutable,execName)
    replaceSearchKeys(cmakeFile,searchKeySource,sourceName)
