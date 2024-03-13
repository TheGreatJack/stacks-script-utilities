import argparse

# Function to parse arguments, take input file, take dictionary file

def parse_arguments():
    """Parses command line arguments using argparse.

    Returns:
        A namespace object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Get population proportions from VCF file and population file")
    parser.add_argument("ustacks_log_file", type=str, help="Path to the catalog tags file")
    parser.add_argument("dictionary_file", type=str, help="Path to the dictionary file")
    return parser.parse_args()


# Function to load dictionary data from a dictionary file, 

def load_dictionary(dictionary_file):
    """Loads dictionary data from a dictionary file.

    Args:
        dictionary_file: Path to the dictionary file.

    Returns:
        A dictionary containing the dictionary data.
    """
    dictionary = {}
    with open(dictionary_file, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            else:
                line = line.strip()
                key, value = line.split('\t')
                dictionary[key] = value
    return dictionary




# "ustacks_log_parser" is a a function to read data from ustacks log data. I want to output the usual data which is: 
# "sample	loci assembled	depth of cov	max cov	number reads incorporated	% reads incorporated"
# I want to add to that: Population (from pop file), a "merge rate" for primary stacks including gaps, secondary reads incorporation
# stdev of cov


def ustacks_log_parser(ustacks_log_file, population_dictionary):
    """Parses the ustacks log file and prints the population proportions for each locus.

    Args:
        ustacks_log_file: Path to the ustacks log file.
        population_dictionary: Population dictionary
    """
    # Open ustacks log file
    with open(ustacks_log_file, 'r') as f:
        # Initalize variables


        # Read file line by line
        while True:
            line = f.readline()
            #print(line)
            if not line: # If for some reason the vcf file is empty
                print("ustacks log file empty???",)
                break
            if line.startswith("Sample"):
                # Print header
                print("sample",
                    "population",
                    "loci_assembled",
                    "depth_of_cov",
                    "stdev_of_cov",
                    "max_cov",
                    "number_reads_incorporated",
                    "percent_reads_incorporated",
                    "primary_reads",
                    "secondary_reads",
                    "stacks_before_primary",
                    "stacks_after_primary",
                    "stacks_before_primary_gaps",
                    "stacks_after_primary_gaps",
                    "secondary_reads_incorp",
                    "secondary_reads_incorp_gaps", sep = "\t")

                sample = line.split(" ")[-1].strip().replace("'","")
                population = population_dictionary[sample]
                while sample !="-" and population!="-":
                    sample, population = gather_sample_data(f, population_dictionary, sample, population)
                break


# Function to get all sample related data, save in variables and then start again with next sample

def gather_sample_data(file_handle, dictionary_file, sample, population):
    """Gets all sample related data, saves in variables and then starts again with next sample.

    Args:
        sample: The sample name.
        population: The population name.
    """
    # Initialize variables
    loci_assembled = "-"
    depth_of_cov = "-"
    stdev_of_cov = "-"
    max_cov = "-"
    number_reads_incorporated = "-"
    percent_reads_incorporated = "-"
    primary_reads = "-"
    secondary_reads = "-"
    stacks_before_primary = "-"
    stacks_after_primary = "-"
    stacks_before_primary_gaps = "-"
    stacks_after_primary_gaps = "-"
    secondary_reads_incorp = "-"
    secondary_reads_incorp_gaps = "-"

    while True:
        line = file_handle.readline()
        if not line: # If end of file is reached
            new_sample = "-"
            new_population = "-"

            print(sample, 
                population, 
                loci_assembled, 
                depth_of_cov, 
                stdev_of_cov, 
                max_cov, 
                number_reads_incorporated, 
                percent_reads_incorporated, 
                primary_reads, 
                secondary_reads, 
                stacks_before_primary, 
                stacks_after_primary, 
                stacks_before_primary_gaps, 
                stacks_after_primary_gaps, 
                secondary_reads_incorp, 
                secondary_reads_incorp_gaps, sep = "\t")

            break
        if line.startswith("Sample"):
            new_sample = line.split(" ")[-1].strip().replace("'","")
            new_population = population_dictionary[new_sample]

            print(sample, 
                population, 
                loci_assembled, 
                depth_of_cov, 
                stdev_of_cov, 
                max_cov, 
                number_reads_incorporated, 
                percent_reads_incorporated, 
                primary_reads, 
                secondary_reads, 
                stacks_before_primary, 
                stacks_after_primary, 
                stacks_before_primary_gaps, 
                stacks_after_primary_gaps, 
                secondary_reads_incorp, 
                secondary_reads_incorp_gaps, sep = "\t")
            break
        elif line.startswith("Final number of"):
            loci_assembled = line.split(" ")[-1].strip()
        elif line.startswith("Final coverage"):
            line_split = line.replace(";","").split(" ")
            depth_of_cov = line_split[2].split("=")[1]
            stdev_of_cov = line_split[3].split("=")[1]
            max_cov = line_split[4].split("=")[1]
            read_data = line_split[5].split("=")[1] 
            number_reads_incorporated = read_data.split("(")[0]
            percent_reads_incorporated = read_data.split("(")[1].split(")")[0].replace("%","")
        elif "stacks representing" in line:
            if "primary" in line:
                primary_reads = line.split(" ")[5].strip()
            elif "secondary" in line:
                secondary_reads = line.split(" ")[6].strip()
        elif line.startswith("Assembling stacks (max."):
            line = file_handle.readline()
            line_split = line.split(" ")
            stacks_before_primary = line_split[3].strip()
            stacks_after_primary = line_split[6].strip().replace(";","")
        elif line.startswith("Assembling stacks, allowing"):
            line = file_handle.readline()
            line_split = line.split(" ")
            stacks_before_primary_gaps = line_split[3].strip()
            stacks_after_primary_gaps = line_split[6].strip()
        elif line.startswith("Merging secondary stacks (max."):
            line = file_handle.readline()
            line_split = line.split(" ")
            secondary_reads_incorp = line_split[3].strip()
        elif line.startswith("Merging secondary stacks, allowing"):
            line = file_handle.readline()
            line_split = line.split(" ")
            secondary_reads_incorp_gaps = line_split[3].strip()
    return new_sample, new_population



if __name__ == "__main__":
    args = parse_arguments()
    population_dictionary = load_dictionary(args.dictionary_file)
    ustacks_log_parser(args.ustacks_log_file, population_dictionary)
