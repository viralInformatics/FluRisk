import os
from collections import defaultdict

import matplotlib
import numpy as np
import pandas as pd

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import argparse
from matplotlib.ticker import FixedLocator


def plot_and_save_combined(protein_data, phenotype_data, filename):
    fig, (ax2, ax1) = plt.subplots(1, 2, figsize = (14, 7))

    ax2.bar(phenotype_data.keys(), phenotype_data.values(),
            color = plt.cm.Spectral(np.linspace(0.15, 1, len(phenotype_data))),
            edgecolor = 'grey', alpha = 0.8, linewidth = 0.5)
    ax2.set_title(f'Marker counts for different phenotypes')
    keys = ["Human-adaptation", "Receptor binding alteration", "Drug resistance",
            "Transmissibility alteration", "Mammalian virulence"]
    # ax2.set_xticklabels(keys, rotation = 45, ha = 'right')

    ax2.set_xticks(range(len(keys)))
    ax2.set_xticklabels(keys, rotation = 45, ha = 'right')
    ax2.xaxis.set_major_locator(FixedLocator(range(len(keys))))


    ax1.bar(protein_data.keys(), protein_data.values(),
            color = plt.cm.viridis(np.linspace(0, 1, len(protein_data))),
            edgecolor = 'grey', alpha = 0.65, linewidth = 0.5)
    ax1.set_title(f'Marker counts for different protein types')
    # ax1.set_xticklabels(protein_data.keys(), rotation = 45, ha = 'right')

    ax1.set_xticks(range(len(protein_data.keys())))
    ax1.set_xticklabels(list(protein_data.keys()), rotation = 45, ha = 'right')
    ax1.xaxis.set_major_locator(FixedLocator(range(len(protein_data.keys()))))

    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


def main():
    parser = argparse.ArgumentParser(description = 'Generate bar plots for protein and phenotype data.')
    parser.add_argument('-i', '--input_dir', required = True, help = 'Input directory containing CSV files.')
    parser.add_argument('-o', '--output_dir', required = True, help = 'Output directory where plots will be saved.')

    args = parser.parse_args()
    input_dir = args.input_dir
    output_dir = args.output_dir

    directories = ['adaptation', 'binding', 'resistance', 'transmissibility', 'virulence']
    result = defaultdict(lambda: defaultdict(set))
    phenotype_dict = defaultdict(dict)

    for directory in directories:
        path = os.path.join(input_dir, directory)
        for filename in os.listdir(path):
            file_path = os.path.join(path, filename)

            if filename.endswith('.csv'):
                df = pd.read_csv(file_path, sep = ',')
                deleted_columns = [f"{directory.title()} Markers", "Protein Type"] if directory != "resistance" else [
                    "Resistance Markers", "Protein Type", "Resistance_level"]
                df.drop_duplicates(subset = deleted_columns, inplace = True)
                file_key = filename.rsplit("_", 1)[0]
                phenotype_dict[file_key][directory] = len(df)

                for _, row in df.iterrows():
                    protein = row[2]
                    if "H3 numbering" in protein:
                        protein = "HA"
                    elif "N2 numbering" in protein:
                        protein = "NA"
                    marker = row[1]
                    result[file_key][protein].add(marker)

    for file_key in result:
        for protein in result[file_key]:
            result[file_key][protein] = len(result[file_key][protein])
    protein_dict = dict(result)
    phenotype_dict = dict(phenotype_dict)

    os.makedirs(output_dir, exist_ok = True)
    image_filenames_combined = []
    for sample, protein_counts in protein_dict.items():
        phenotype_counts = phenotype_dict[sample]
        filename = f'{sample.replace(" ", "_").replace("(", "").replace(")", "").replace(".", "")}-combined.png'
        plot_and_save_combined(protein_counts, phenotype_counts, os.path.join(output_dir, filename))
        image_filenames_combined.append(filename)


if __name__ == '__main__':
    main()
    # main("./","total_barplot")
