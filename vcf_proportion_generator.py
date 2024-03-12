import argparse

# Function to parse arguments, take input file, take dictionary file

def parse_arguments():
    """Parses command line arguments using argparse.

    Returns:
        A namespace object containing parsed arguments.
    """
    parser = argparse.ArgumentParser(description="Get population proportions from VCF file and population file")
    parser.add_argument("vcf_file", type=str, help="Path to the dictionary file")
    parser.add_argument("dictionary_file", type=str, help="Path to the dictionary file")
    return parser.parse_args()


# Function to load dictionary data from a dictionary file

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

# "vcf_reader" is a a function to read data from VCF file. Line by line, when all lines that belong to a locus are detected another function estimates the 
# proportion of each population in that locus, using the dictionary data. The special function prints the Locus and proportion of each population. 

def vcf_reader(vcf_file, population_dictionary):
    """Loads VCF data from a VCF format file.

    Args:
        input_file: Path to the VCF format file.

    """
    with open(vcf_file, 'r') as f:
        # Read file line by line

        # Header reading
        while True:
            line = f.readline()
            if not line: # If for some reason the vcf file is empty
                print("VCF file empty???",)
                break
            if line.startswith("#CHROM"):
                vcf_sample_names = line.strip().split("\t")
                vcf_sample_names = vcf_sample_names[9:] #Just get the sample names
                sample_number = len(vcf_sample_names)
                print("#Sample_number:", sample_number)
                unique_populations = list(set(population_dictionary.values()))
                unique_populations.sort()
                break
        
        print("Locus", "sample_proportion","sample_max" , "sample_average",*unique_populations, sep="\t")
        
        # Initialize data with the first reading jaja
        locus_data = []
        ns_data = []

        line = f.readline()
        # Split line by tab
        line = line.strip().split("\t")
        previous_locus = line[0]
        info_field = info_field_extractor(line[7]) 
        data_of_locus = line[9:]

        locus_data.append(data_of_locus)
        ns_data.append(info_field)   

        while True:
            line = f.readline()
            if not line:
                # Print the last locus
                ns_average,locus_total_samples,population_proportions = locus_data_processor(locus_data, ns_data, vcf_sample_names , unique_populations,population_dictionary)
                # Print result of locus 
                print(previous_locus,locus_total_samples/sample_number, ns_average,*population_proportions.values(), sep="\t")
                break        
            # Split line by tab
            line = line.strip().split("\t")
            current_locus = line[0]
            info_field = info_field_extractor(line[7]) 
            data_of_locus = line[9:]
            if current_locus != previous_locus:
                # Proccess locus_data and ns_data
                ns_average,locus_total_samples,population_proportions = locus_data_processor(locus_data, ns_data, vcf_sample_names , unique_populations,population_dictionary)
                # Print result of locus 
                print(previous_locus,locus_total_samples/sample_number, locus_total_samples,ns_average,*population_proportions.values(), sep="\t")
                # Re define locus_data and ns_data with the current data_of_locus and info field
                previous_locus = current_locus
                locus_data = []
                ns_data = []
                locus_data.append(data_of_locus)
                ns_data.append(info_field)
            else:
                locus_data.append(data_of_locus)
                ns_data.append(info_field)      
                

# Small info field handler to output the NS field easily

def info_field_extractor(info_field):
    """Extracts the NS field from the INFO field of a VCF record.

    Args:
        info_field: The INFO field of a VCF record.

    Returns:
        The NS field of the VCF record.
    """
    ns_field = None
    for field in info_field.split(";"):
        if field.startswith("NS="):
            ns_field = field.split("=")[1]
            break
    return ns_field



# Locus data processer

def locus_data_processor(locus_data, ns_data, vcf_sample_names , unique_populations,population_dictionary):
    """Calculates the proportion of each population in a locus.

    Args:
        locus_data: A list of lists containing the data for each sample at a locus.
        ns_data: A list of strings containing the NS field for each sample at a locus.
        vcf_sample_names: A list of strings containing the names ofunique_populations = list(set(population_dictionary.values())) the samples.
        population_dictionary: A dictionary containing the population of each sample.

    Returns:
        A dictionary containing the proportion of each population in the locus.
    """
    # Get highest number in ns data
    ns_data = [int(x) for x in ns_data]
    ns_average = sum(ns_data) / len(ns_data)

    population_proportions = {}

    # Create dictionary to store population proportions
    population_proportions = dict.fromkeys(unique_populations, 0)

    # Iterate over samples and calculate proportions
    for sample_info in zip(vcf_sample_names, *locus_data):
        #print(sample_info)
        sample_name = sample_info[0]
        data = sample_info[1:]

        # Check for samples with non-empty data
        if any(x != "./." for x in data):
            #print("samples_with_data:",sample_name, data) 
            population = population_dictionary[sample_name]
            population_proportions[population] = population_proportions[population] + 1 
            #print("population", population,sample_name,population_proportions[population])
    # Normalize proportions
    locus_total_samples = sum(population_proportions.values())
    for population, proportion in population_proportions.items():
        print("norm",proportion,    locus_total_samples)
        population_proportions[population] = proportion / locus_total_samples
    
    # Return population proportions<ctrl63> 
    return ns_average,locus_total_samples,population_proportions



if __name__ == "__main__":
    args = parse_arguments()
    population_dictionary = load_dictionary(args.dictionary_file)
    vcf_reader(args.vcf_file, population_dictionary)
