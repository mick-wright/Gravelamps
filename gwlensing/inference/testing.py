import sys
import os
from configparser import ConfigParser

def main():
    #Instantiate the Configuration Parser
    config = ConfigParser()

    #If user hasn't given a useable ini file, raise exception
    if not os.path.isfile(sys.argv[1]):
        raise IOError("Input ini file not given")

    #Check that the Configuration Parser can read the ini file
    try:
        config.read("sys.argv[1]")
    except IOError:
        print("Input file cannot be read") 
