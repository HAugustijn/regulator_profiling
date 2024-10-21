#!/usr/bin/env python3

# Import statements
import argparse
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def get_arguments():
    """Parsing the arguments"""
    parser = argparse.ArgumentParser(description="",
                                     usage='''
-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
## python visualize_results.py
Mandatory arguments:
    -I  Provide the ratio in/out BGC file (pfam_in_out_bgc.tsv)
    -C  Provide the ratio BGC class file (class_bgc_counts.tsv)
    -Cc Provide the ratio file with the BGC subclasses (class_bgc_counts_type2.tsv)
    -B  Provide the MIBiG classification file (class_bgc_counts_type3.tsv)
    -R  Provide the ratio CDSs in BGCs file (ratio_CDSs_in_BGCs.tsv)
    -M  Provide the MIBiG annotation file (BGC_types_mibig_sortedV2.txt)
    -O  Provide the path to the output folder where results should be deposited

-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-_-
''')

    parser.add_argument("-I", "--input", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-C", "--classes", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-Cc", "--classesc", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-B", "--bgc_counts", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-R", "--ratio", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-M", "--mibig", help=argparse.SUPPRESS, required=True)
    parser.add_argument("-O", "--outdir", help=argparse.SUPPRESS, required=True)

    return parser.parse_args()

def main():
    args = get_arguments()

    df = pd.read_csv(args.input, sep='\t')
    df_classes = pd.read_csv(args.classes, sep='\t')
    df_classesc = pd.read_csv(args.classesc, sep='\t', index_col=0)
    df_bgc_counts = pd.read_csv(args.bgc_counts, sep='\t')
    df_bgc_ratio = pd.read_csv(args.ratio, sep='\t')
    df_mibig = pd.read_csv(args.mibig, sep='\t')

    # Make box plot BGC ratio
    sns.boxplot(y=df_bgc_ratio['3'], width=.3)
    plt.ylabel('Ratio hits in BGCs (%)')
    plt.gcf().set_size_inches(2, 4)
    plt.ylim(-5, 100)
    plt.savefig(f"{args.outdir}ratio_bgcs_box_all.svg", format="svg")

    # Identifying outliers in the second dataset
    Q1 = df['%'].quantile(0.25)
    Q3 = df['%'].quantile(0.75)
    IQR = Q3 - Q1

    lower_bound = Q1 - 1.5 * IQR
    upper_bound = Q3 + 1.5 * IQR

    outliers = df[(df['%'] < lower_bound) | (df['%'] > upper_bound)]

    # Merge df and df_classes on the category column (assumed to be the first column in df_classes)
    merged_df = pd.merge(df_bgc_counts, df, left_on=df_bgc_counts.columns[0], right_on='Unnamed: 0')
    merged_df_classes = pd.merge(df_classesc, df['%'], left_on=df_classesc.columns[0], right_index=True)


    # Filter the merged dataframe to only include outliers
    merged_outliers = merged_df[merged_df['Unnamed: 0'].isin(outliers['Unnamed: 0'])]
    merged_outliers_classes = merged_df_classes[(merged_df_classes.index).isin(outliers['Unnamed: 0'])]

    # Sort the merged_outliers dataframe by the % column from high to low
    sorted_outliers = merged_outliers.sort_values(by='%', ascending=True).reset_index(drop=True)
    sorted_outliers = sorted_outliers.loc[:, (sorted_outliers != 0).any(axis=0)]

    sorted_outliers_c = merged_outliers_classes.sort_values(by='%', ascending=True)
    sorted_outliers_c = sorted_outliers_c.loc[:, (sorted_outliers_c != 0).any(axis=0)]

    sort_list = ["NRPS", "NRP-metallophore", "CDPS", "isocyanide-nrp", "NAPAA", "prodigiosin", "thioamide-NRP",
                 "NRPS-like", "T1PKS", "T2PKS", "HR-T2PKS", "T3PKS", "transAT-PKS-like", "transAT-PKS", "hglE-KS",
                 "arylpolyene", "PKS-like", "lassopeptide", "lanthipeptide-class-i", "lanthipeptide-class-ii",
                 "lanthipeptide-class-iii", "lanthipeptide-class-iv", "lanthipeptide-class-v", "crocagin",
                 "cyanobactin", "guanidinotides", "LAP", "linaridin", "lipolanthine", "methanobactin", "ranthipeptide",
                 "redox-cofactor", "thioamitides", "thiopeptide", "triceptide", "RRE-containing", "RiPP-like",
                 "terpene", "oligosaccharide", "2dos", "amglyccycl", "aminocoumarin", "aminopolycarboxylic-acid",
                 "betalactone", "blactam", "butyrolactone", "ectoine", "furan", "hserlactone", "hydrogen-cyanide",
                 "indole", "melanin", "NI-siderophore", "nucleoside", "phenazine", "phosphoglycolipid", "phosphonate",
                 "pyrrolidine", "resorcinol", "other"]

    filtered_df_c = sorted_outliers_c.reindex(sort_list, axis=1)
    filtered_df_c = filtered_df_c.loc[:, filtered_df_c.sum(axis=0) > 0]
    normalized_classes = filtered_df_c.copy()

    def adjust_value(value):
        if value >= 75:
            return 100
        elif value > 0 and value < 25:
            return 25
        elif value >= 25 and value < 50:
            return 50
        elif value >= 50 and value < 75:
            return 75
        else:
            return 0

    # Applying the function to each element in the DataFrame
    for index, row in normalized_classes.iterrows():
        total_count = row.sum()
        avg_row = row / total_count * 100
        normalized_classes.loc[index] = avg_row.apply(adjust_value)

    sorted_normalized_classes = normalized_classes.reindex(list(sorted_outliers['Unnamed: 0']))

    sorted_normalized_classes.to_csv(f"{args.outdir}/class_BGCs_itol.tsv", sep="\t")

    bgc_to_functions = {}
    for _, row in df_mibig.iterrows():
        bgc = row['Accession']
        functions = row[row == True].index.tolist()
        bgc_to_functions[bgc] = functions

    transformed_data = {}
    for index_counts, row_counts in df_bgc_counts.iterrows():
        results_dict = row_counts.to_dict()
        reg_name = results_dict['Unnamed: 0']
        if reg_name not in transformed_data.keys():
            transformed_data[reg_name] = {}

        for key in results_dict.keys():
            if not key == 'Unnamed: 0':
                if 'Unnamed' in key or 'NA' in key:
                    function = 'NA'
                else:
                    function = bgc_to_functions[key]
                if isinstance(function, list):
                    for func in function:
                        if func not in transformed_data[reg_name].keys():
                            transformed_data[reg_name][func] = results_dict[key]
                        else:
                            transformed_data[reg_name][func] += results_dict[key]
                else:
                    if function not in transformed_data[reg_name].keys():
                        transformed_data[reg_name][function] = results_dict[key]
                    else:
                        transformed_data[reg_name][function] += results_dict[key]

    transformed_df = pd.DataFrame(transformed_data).T
    filtered_df_merged = transformed_df[(transformed_df.index).isin(sorted_outliers['Unnamed: 0'])]
    filtered_df = filtered_df_merged.loc[:, filtered_df_merged.sum(axis=0) > 10]
    filtered_df = filtered_df.reindex(sorted(filtered_df.columns), axis=1)
    filtered_df = filtered_df.reindex(list(filtered_df.columns[1:]) + ['NA'], axis=1)

    normalized_classes = filtered_df.copy()
    for index, row in normalized_classes.iterrows():
        total_count = row.sum()
        normalized_classes.loc[index, normalized_classes.columns] = row / total_count * 100

    df_hm = normalized_classes.reindex(list(sorted_outliers['Unnamed: 0']))
    df_hm.to_csv(f"{args.outdir}/class_BGCs_MIBiG.tsv", sep="\t")

    categories = normalized_classes.columns
    fig, ax = plt.subplots(figsize=(5, 6))

    color_ranges = ["#B2D58C", "#115366", "#F9F1A5", "#85D2DC", "#FABD68", "#9D323F", "#6A8E3C", "#FF69B4",
                    "#8A2BE2", "#7FFF00", "#DC143C", "#00CED1", "#FF8C00", "#8B0000", "#00FA9A", "#8B008B",
                    "#00FF7F", "#7B68EE", "#FF4500", "#4682B4", "#D2691E", "#FF1493", "#6495ED", "#cccccc"]

    bottom = [0] * len(df_hm)
    x_labels = df_hm.index
    for i, category in enumerate(categories):
        values = df_hm[category]
        ax.bar(x_labels, values, bottom=bottom, label=category, color=color_ranges[i % len(color_ranges)], width=0.7)
        bottom += values

    ax.set_xticklabels(x_labels, rotation=45, ha='right')
    ax.legend(title='Classes', bbox_to_anchor=(1.05, 1), loc='upper left')

    plt.savefig(f"{args.outdir}mibig_classes.svg", format="svg")

    plt.tight_layout()
    plt.show()

    # Merge df and df_classes on the category column (assumed to be the first column in df_classes)
    merged_df = pd.merge(df_classes, df, left_on=df_classes.columns[0], right_on='Unnamed: 0')

    # Filter the merged dataframe to only include outliers
    merged_outliers = merged_df[merged_df['Unnamed: 0'].isin(outliers['Unnamed: 0'])]

    # Sort the merged_outliers dataframe by the % column from high to low
    sorted_outliers = merged_outliers.sort_values(by='%', ascending=True).reset_index(drop=True)

    # Normalize the class counts to match the percentage scores from df
    normalized_classes = sorted_outliers.copy()
    for index, row in normalized_classes.iterrows():
        percentage = row['%']
        total_count = row[1:7].sum()
        normalized_classes.loc[index, normalized_classes.columns[1:7]] = row[1:7] * (percentage / total_count)

    # Create a stacked bar plot for the normalized outlier data
    categories = df_classes.columns[1:]
    fig, ax = plt.subplots(figsize=(6, 10))

    # Plotting the horizontal stacked bars for normalized outliers
    bottom = [0] * len(normalized_classes)
    y_labels = normalized_classes[normalized_classes.columns[0]]
    for category in categories:
        values = normalized_classes[category]
        ax.barh(y_labels, values, left=bottom, label=category)
        bottom += values

    ax.set_xlim(0, 110)
    ax.set_ylabel('Categories')
    ax.set_xlabel('Percentage')
    ax.set_title('Horizontal Stacked Bar Plot of Outlier Categories')
    ax.legend(title='Classes', bbox_to_anchor=(1.05, 1), loc='upper left')

    # Add labels for the outliers
    def addlabels(y, x, k):
        for i in range(len(y)):
            ax.text(x[i] + 0.5, i, f"n={k[i]}", va='center', fontsize=7, alpha=0.5)

    addlabels(sorted_outliers['Unnamed: 0'], sorted_outliers['%'], sorted_outliers['sum'])

    plt.savefig(f"{args.outdir}ratio_bgc_classes.svg", format="svg")

    plt.tight_layout()
    plt.show()


if __name__ == "__main__":
    main()
