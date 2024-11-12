#!/usr/bin/python

# Macro to create lists for new upcDst production 
# from embedding MuDst proccesing

import os
from subprocess import Popen
from os.path import exists
import re
import sys
import time


def list_files_in_directory(directory_path):
    try:
        # Get a list of all files and directories in the specified directory
        files_and_directories = os.listdir(directory_path)
        
        # Filter out directories, keeping only files
        files = [f for f in files_and_directories if os.path.isfile(os.path.join(directory_path, f))]
        
        return files
    except FileNotFoundError:
        return "The directory '{}' does not exist.".format(directory_path)
    except Exception as e:
        return "An error occurred: {}".format(e)

#_____________________________________________________________________________
if __name__ == "__main__":

    lists_path = os.getcwd() + "/lists_proton"
    with open("config_template.in", 'r') as source_file, open("config_ProtonProtonbar.in", 'w') as dest_file:
        for line in source_file:
            if 'add_input' in line:
                for run in list_files_in_directory(lists_path):
                    run_number = run[0:run.rfind(".")]
                    input_line = "add_input " + run_number + " " + os.path.abspath(run) + "\n"
                    dest_file.write(input_line)
            else:
                dest_file.write(line)  # Write the original line















