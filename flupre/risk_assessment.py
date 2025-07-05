# !/usr/bin/python
# -*- coding: utf-8 -*-
# Created on 2024/3/25 13:09
import argparse
import os
import pickle
from collections import defaultdict

import numpy as np
import pandas as pd

base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MARKER1 = os.path.join(base_dir, 'data', 'markers_dict_four.pickle')
MARKER2 = os.path.join(base_dir, 'data', 'markers_dict_resistance.pickle')
phenotype_descriptions = {
    'transmissibility': 'Mammalian Transmissibility',
    'binding': 'Human Receptor Binding Ability',
    'adaptation': 'Mammalian Adaptability',
    'virulence': 'Mammalian Virulence'
}

drugs = ["oseltamivir", "zanamivir", "peramivir", "laninamivir", "baloxavir", "adamantane"]  # "laninamivir",
level_order = {'unknown': 0, 'low': 1, 'middle': 2, 'high': 3}
resistance_weights = {
    'high': 4,
    'middle': 3,
    'low': 2,
    'unknown': 1
}

def calculate_relative_position(sample_value, possible_values):
    if sample_value == 0:
        return 0.0
    values, counts = np.unique(possible_values, return_counts = True)
    non_zero_values = values[values > 0]
    non_zero_counts = counts[values > 0]
    adjusted_non_zero_cdf = np.cumsum(non_zero_counts) / np.sum(non_zero_counts)
    non_zero_values = np.insert(non_zero_values, 0, 0)
    adjusted_non_zero_cdf = np.insert(adjusted_non_zero_cdf, 0, 0)
    return np.interp(sample_value, non_zero_values, adjusted_non_zero_cdf)


def replace_prob_host(row):
    if row['Host'] == 'Avian':
        row['Probability'] = 1 - float(row['Probability'])
    return row

def replace_prob_virulence(row):
    if row['Virulence Level'] == 'Avirulent':
        row['Probability'] = 1 - float(row['Probability'])
    return row

def replace_prob_binding(row):
    if row['Binding Type'] == 'α2-3':
        row['Probability'] = 1 - float(row['Probability'])
    return row

def generate_phenotype_stats(input_dir, annotation_dir, host_dir, virulence_dir,binding_dir):
    phenotype_dict = defaultdict(dict)
    phenotype_categories = ['transmissibility', 'binding', 'adaptation', 'virulence', 'resistance']


    host_df = pd.read_csv(f"{host_dir}/prediction.csv")
    host_df = host_df.apply(replace_prob_host,axis = 1)
    host_df.drop_duplicates(inplace = True)

    binding_df = pd.read_csv(f"{binding_dir}/prediction.csv")
    binding_df = binding_df.apply(replace_prob_binding,axis = 1)
    binding_df.drop_duplicates(inplace = True)

    vir_df = pd.read_csv(f"{virulence_dir}/prediction.csv")
    vir_df = vir_df.apply(replace_prob_virulence,axis = 1)
    for category in phenotype_categories:
        path = os.path.join(input_dir, category)

        for filename in os.listdir(path):
            file_key = filename.split('_markers')[0]
            file_path = os.path.join(path, filename)
            drugs = ["oseltamivir", "zanamivir", "peramivir", "laninamivir", "baloxavir", "adamantane"]
            level_order = {'unknown': 0, 'low': 1, 'middle': 2, 'high': 3}
            resistance_weights = {
                'high': 4,
                'middle': 3,
                'low': 2,
                'unknown': 1
            }

            if filename.endswith('_markers.csv'):
                df = pd.read_csv(file_path)
                if category == "resistance":
                    df['level_order'] = df['Resistance_level'].map(level_order)
                    for drug in drugs:
                        drug_df = df[df['Drug'] == drug].copy()
                        if not drug_df.empty:
                            drug_df_filtered = drug_df.sort_values('level_order', ascending = False).drop_duplicates(
                                subset = ["Resistance Markers", "Protein Type"], keep = 'first'
                            )
                            weighted_count = drug_df_filtered['Resistance_level'].map(resistance_weights).sum()
                            phenotype_dict[file_key][drug] = weighted_count
                        else:
                            phenotype_dict[file_key][drug] = 0

                    df.drop(columns = 'level_order', inplace = True)
                else:
                    df.drop_duplicates(subset = [f"{category.title()} Markers", "Protein Type"], inplace = True)
                    phenotype_dict[file_key][category] = len(df)
    annotation_dict = read_for_anno_files(annotation_dir)
    for df, key in zip([host_df, binding_df, vir_df], ['adaptation', 'binding', 'virulence']):
        for _, row in df.iterrows():
            strain_id = row['Strain ID']
            probability = row['Probability']
            if strain_id in phenotype_dict:
                phenotype_dict[strain_id][key] = probability
    phenotype_dict = update_phenotype_dict(phenotype_dict, annotation_dict)
    # print(phenotype_dict)
    return phenotype_dict

def update_phenotype_dict(phenotype_dict, annotation_dict):
    for strain_id in phenotype_dict:
        if strain_id not in annotation_dict or annotation_dict[strain_id] == ['Unknown']:
            phenotype_dict[strain_id]['virulence'] = 0
            phenotype_dict[strain_id]['adaptation'] = 0

    return phenotype_dict

def read_for_anno_files(directory):
    result_dict = {}

    files = [f for f in os.listdir(directory) if f.endswith('_for_anno.csv')]
    for file in files:
        key = file.replace('_annotation_for_anno.csv', '')
        file_path = os.path.join(directory, file)
        df = pd.read_csv(file_path)

        if 'Protein Abbreviation' not in df.columns:
            unique_proteins = df['ProteinType'].drop_duplicates().tolist()
        else:
            unique_proteins = df['Protein Abbreviation'].drop_duplicates().tolist()
        result_dict[key] = unique_proteins

    return result_dict



def save_risk_values_to_csv(sample_values, four_markers, resistance_markers, output_dir, test_name, all_data, host_dir,
                            virulence_dir):
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    phenotypes = list(phenotype_descriptions.keys())

    relative_positions_four = {
        phenotype: calculate_relative_position(sample_values[phenotype],
                                               four_markers[phenotype_descriptions[phenotype]])
        for phenotype in phenotypes if phenotype not in ['adaptation', 'virulence','binding']
    }

    for phenotype in ['adaptation', 'virulence','binding']:
        if phenotype in phenotypes:
            relative_positions_four[phenotype] = sample_values[phenotype]

    relative_positions_resistance = {
        drug: calculate_relative_position(sample_values[drug], resistance_markers[drug])
        for drug in resistance_markers.keys()
    }

    risk_values = {**relative_positions_four, **relative_positions_resistance}
    # print(risk_values)

    df = pd.DataFrame([risk_values])
    df = df.rename(phenotype_descriptions, axis = 1)

    csv_file_path = os.path.join(output_dir, f"{test_name}_risk_values.csv")
    df.to_csv(csv_file_path, index = False)

    df['Isolate'] = test_name.replace('_risk', '')
    df = df[['Isolate'] + [col for col in df.columns if col != 'Isolate']]
    df.columns = ['Isolate', "Mammalian-transmissibility", "Human-receptor-binding", "Mammalian-adaptation",
                  "Mammalian-virulence", "Oseltamivir", "Zanamivir", "Peramivir", "Laninamivir","Baloxavir", "Adamantane"]
    df = df.loc[:, ['Isolate', "Mammalian-adaptation", "Mammalian-virulence", "Mammalian-transmissibility",
                    "Human-receptor-binding", "Oseltamivir", "Zanamivir", "Peramivir", "Laninamivir", "Baloxavir", "Adamantane"]]
    all_data.append(df)



def main():
    parser = argparse.ArgumentParser(description = 'Generate risk values CSV for markers data.')
    parser.add_argument('-i', '--input_dir', required = True, help = 'Input directory containing Markers CSV files.')
    parser.add_argument('-a', '--annotation_dir', required = True, default = 'result/', help = 'Input directory containing sequence annotation CSV files.')
    parser.add_argument('-hp', '--host_dir', default = 'ada_prediction/', help = 'Input directory containing host prediction files.')
    parser.add_argument('-vp', '--virulence_dir', default = 'vir_prediction/',
                        help = 'Input directory containing virulence prediction files.')
    parser.add_argument('-bp', '--binding_dir', default = 'bin_prediction/',
                        help = 'Input directory containing receptor binding perfernece prediction files.')
    parser.add_argument('-o', '--output_dir', required = True, help = 'Output directory where CSV files will be saved.')
    args = parser.parse_args()

    input_dir = args.input_dir
    annotation_dir = args.annotation_dir
    output_dir = args.output_dir
    host_dir = args.host_dir
    virulence_dir = args.virulence_dir
    binding_dir = args.binding_dir

    all_data = []
    with open(MARKER1, 'rb') as handle:
        loaded_markers_four = pickle.load(handle)
    with open(MARKER2, 'rb') as handle:
        loaded_markers_res = pickle.load(handle)
    phenotype_dict = generate_phenotype_stats(input_dir, annotation_dir, host_dir, virulence_dir, binding_dir)
    # print(phenotype_dict)
    for test_name, sample_values in phenotype_dict.items():
        save_risk_values_to_csv(sample_values, loaded_markers_four, loaded_markers_res, output_dir, test_name, all_data,
                                host_dir, virulence_dir)
    if all_data:
        summary_df = pd.concat(all_data)
        summary_df.to_csv(os.path.join(output_dir, 'summary_risk_values.csv'), index = False)


if __name__ == "__main__":
    main()
