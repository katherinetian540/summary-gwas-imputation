import pandas as pd

def add_chr_prefix(file_path, output_path):
    # Read the file
    data = pd.read_csv(file_path, sep='\t', nrows=1000)

    data['chr'] = 'chr' + data['chr'].astype(str)
    # Save the modified data to a new file
    data.to_csv(output_path, sep='\t', index=False)

if __name__ == "__main__":
    input_file = '/Users/katherinetian/Downloads/oncoarray_bcac_public_release_oct17.txt'  # Path to your input GWAS file
    output_file = '/Users/katherinetian/Downloads/gwas_subset_with_chr.txt'  # Output file path

    # Apply function to add 'chr' to chromosome and save the subset
    add_chr_prefix(input_file, output_file)
    
    # Load and print the head of the new subset to verify changes
    data = pd.read_csv(output_file, sep='\t')
    print(data.head())
