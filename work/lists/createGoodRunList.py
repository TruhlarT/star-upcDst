#!/usr/bin/python

def extract_run_number(line):
    # Split the line by '/' and get the last element
    parts = line.split('/')
    # The run number is the last part without the file extension
    run_number = parts[-1].split('.')[0]
    #print run_number
    return run_number

def main():
    # Open the first file for reading
    with open('extendedUpc.list', 'r') as file1:
        lines_file1 = file1.readlines()

    # Open the second file for reading
    with open('mcTest.list', 'r') as file2:
        run_numbers_file2 = set(file2.read().splitlines())

    #print run_numbers_file2

    # Open a new file for writing
    with open('forwardRPTest.list', 'w') as output_file:
        for line in lines_file1:
            # Extract the run number from the line
            run_number = extract_run_number(line)
            # Check if the run number is in the second file
            substring_found = False
            for pathToFile in run_numbers_file2:
                if run_number in pathToFile:
                    substring_found = True
                    break

            # If so, write the line to the output file
            if substring_found:
                #print "yes"
                output_file.write(line)

if __name__ == "__main__":
    main()



