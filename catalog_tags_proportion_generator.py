import argparse
import gzip
import numpy as np


# Function to parse arguments, take input file, take dictionary file

def parse_arguments():
    """Parses command line arguments using argparse.

    Returns:
        A namespace object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Get population proportions from VCF file and population file")
    parser.add_argument("catalog_tags_file", type=str, help="Path to the catalog tags file")
    parser.add_argument("dictionary_file", type=str, help="Path to the dictionary file")
    return parser.parse_args()


# Function to load dictionary data from a dictionary file, 
# adds a 1-based index to each item in the file dictionary

def load_dictionary(dictionary_file):
    """Loads dictionary data from a dictionary file.

    Args:
        dictionary_file: Path to the dictionary file.

    Returns:
        A dictionary containing the dictionary data.
    """
    dictionary = {}
    index = 1 
    with open(dictionary_file, 'r') as f:
        while True:
            line = f.readline()
            if not line:
                break
            else:
                line = line.strip()
                key, value = line.split('\t')
                #dictionary[key] = value
                dictionary[index] = value
                index += 1
    return dictionary

# Function to read gzip catalog.tags.tsv line by line

def read_catalog_tags_file(catalog_tags_file, population_dictionary):
    """Reads catalog.tags.tsv file line by line.

    Args:
        catalog_tags_file: Path to the catalog.tags.tsv file.

    Returns:
        A list of lines from the catalog
    """

    unique_populations = list(set(population_dictionary.values()))
    unique_populations.sort()

    print("Locus", "Locus_diversity", "Sample_number", "Locus_length", *unique_populations, sep="\t")
    with gzip.open(catalog_tags_file, 'rt') as f:
        line = f.readline() # Skip first line
        if not line.startswith("#"):
            print("File doesnt start with comment")
        while True:
            line = f.readline()
            if not line: # If for some reason the vcf file is empty
                print("Catalog file empty???",)
                break
            line = line.strip().split("\t")
            locus = line[1]
            sample_provenance = line[4]
            sample_provenance = sample_provenance_extractor(sample_provenance, population_dictionary)
            locus_diversity = simpsons_diversity(sample_provenance)
            locus_length = len(line[5])
            samples_number = len(sample_provenance)
            pop_proportions = get_population_proportions(sample_provenance, unique_populations)

            print(locus, locus_diversity, samples_number, locus_length, *pop_proportions.values(), sep="\t")

# Function to handle sample provenance of a particular locus
            
def sample_provenance_extractor(sample_provenance_raw, population_dictionary):
    population_list = []

    sample_provenance = sample_provenance_raw.split(",")

    for sample in sample_provenance:
        sample = sample.split("_")[0]
        population = population_dictionary[int(sample)]
        population_list.append(population)

    return population_list

# Function to estimate the population diversity in the locus

def simpsons_diversity(data):
    """
    Calculates Simpson's Diversity Index for a given data set.

    Args:
        data: A list or NumPy array containing data points.

    Returns:
        The Simpson's Diversity Index as a float.
    """
    data_counts = np.unique(data, return_counts=True)[1]
    data_counts_1 = data_counts - 1

    numerator = np.sum(data_counts * data_counts_1)
    denom_1 = np.sum(data_counts)
    denom = denom_1 * (denom_1 -1)
    return 1 - numerator/denom


# Function to get population proportions from a list with population classifications

def get_population_proportions(population_list, unique_populations):
    """
    Calculates the population proportions for a given list of population classifications.

    Args:
        population_list: A list of population classifications.

    Returns:
        A dictionary containing the population proportions.
    """
    # Initiatilize dictionary
    population_proportions = dict.fromkeys(unique_populations, 0)

    # Add counts
    for population in population_list:
        population_proportions[population] += 1

    # Normalize proportions
    locus_total_samples = sum(population_proportions.values())
    for population, proportion in population_proportions.items():
        population_proportions[population] = proportion / locus_total_samples
    
    return population_proportions




if __name__ == "__main__":
    args = parse_arguments()
    population_dictionary = load_dictionary(args.dictionary_file)
    read_catalog_tags_file(args.catalog_tags_file, population_dictionary)
