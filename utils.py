import numpy as np
import pandas as pd
import re

species_dict = {"Bd": 'Brevundimonas_diminuta_HAMBI_18',
                "Ct": 'Comamonas_testosteroni_HAMBI_403',
                "Ec": 'Escherichia_coli_ATCC_11303',
                "Jl": 'Janthinobacterium_lividum_HAMBI_1919',
                "Sm": 'Serratia_marcescens_ATCC_13880',
                "Sc": 'Sphingomonas_capsulata_HAMBI_103',
                "Pf": 'Pseudomonas_fluorescens_SBW25'}

exp_dict = {"PS": "bacteria-ciliate",
            "NP": "bacteria_alone"}

def pos_to_char(pos):
    return chr(pos + 97)

def read_replicate_variant_data(replicate, f, sampling_data_file, min_2_consecutive_threshold=0.1, day_cutoff=800):  
    
    sampling_df = pd.read_csv(sampling_data_file, sep="\t")
    org = species_dict[replicate.split("_")[1]]
    exp = exp_dict[replicate.split("_")[0]]

    sampling_df_pred = sampling_df[(sampling_df["bacterial_strain"] == org) & (sampling_df["organism_type"] == 'ciliate')]
    sampling_df = sampling_df[(sampling_df["bacterial_strain"] == org) & (sampling_df["organism_type"] == 'bacteria')]
    sampling_df_rep = sampling_df[(sampling_df["replicate"] == int(replicate.split("_")[-1][1:])) & (sampling_df["experiment"] == exp)]
    sampling_df_pred_rep = sampling_df_pred[(sampling_df_pred["replicate"] == int(replicate.split("_")[-1][1:])) & (sampling_df_pred["experiment"] == exp)]
  
    header = None
    with open(f, "r") as file:
        for line in file:
            if re.search("^#CHROM", line):
                header = line.split("\t")
    header[0] = header[0][1:]
    header[-1] = header[-1][:-1]

    df = pd.read_csv(f, sep="\t", comment="#", header=None)
    df.columns = header
    samples = df.columns.values[9:]

    chosen_days = []
    chosen_samples = []
    chosen_dates = []

    od600 = []
    pred_counts = []

    for x in samples:
        if ("PS" in x or "NP" in x) and "merged" in x:
            date_split = x.split("_")[3:6]
            new_date = date_split[2] + "/" + date_split[1] + "/" + date_split[0]
        elif "Pf_BP" in x or "Pf_B_" in x:
            new_date = x.split("_")[3]
        elif ("PS" in x or "NP" in x) and not "merged" in x:
            date_split = x.split("_")[3:]
            new_date = date_split[2] + "/" + date_split[1] + "/" + date_split[0]

        if "D" in new_date:
            ups = sampling_df_rep[sampling_df_rep["day_in_experiment"] == int(new_date[1:])]["unreliable_pop_size"]
            day = new_date[1:]
            od = sampling_df_rep[sampling_df_rep["day_in_experiment"] == int(new_date[1:])]["OD600"]
            pred_count = sampling_df_pred_rep[sampling_df_pred_rep["day_in_experiment"] == int(new_date[1:])]["pred_corrected"]
        else:
            ups = sampling_df_rep[sampling_df_rep["sampling_date"] == new_date]["unreliable_pop_size"]
            day = sampling_df_rep[sampling_df_rep["sampling_date"] == new_date]["day_in_experiment"]
            od = sampling_df_rep[sampling_df_rep["sampling_date"] == new_date]["OD600"]
            pred_count = sampling_df_pred_rep[sampling_df_pred_rep["sampling_date"] == new_date]["pred_corrected"]

        if len(ups) == 1:
            if ups.values[0] == "no":
                if type(day) == str:
                    chosen_days.append(int(day))
                else:
                    chosen_days.append(day.values[0])
                chosen_samples.append(x)
                chosen_dates.append(new_date)
                od600.append(od.values[0])

                if len(pred_count) != 0:
                    pred_counts.append(pred_count.values[0])
                else:
                    pred_counts.append(0)

    flag_for_deletion = np.array([False if len(x.split(",")) != 1 else True for x in df["ALT"]])
    df = df[flag_for_deletion]
    df = df[list(df.columns.values[:9]) + chosen_samples]

    day_cutoff_idx = np.where(np.array(chosen_days) < day_cutoff)[0]
    chosen_dates = list(np.array(chosen_dates)[day_cutoff_idx])
    chosen_days = list(np.array(chosen_days)[day_cutoff_idx])
    chosen_samples = list(np.array(chosen_samples)[day_cutoff_idx])
    od600 = list(np.array(od600)[day_cutoff_idx])
    pred_counts = list(np.array(pred_counts)[day_cutoff_idx])

    chosen_dates = [x for _,x in sorted(zip(chosen_days, chosen_dates))]
    chosen_samples = [x for _,x in sorted(zip(chosen_days, chosen_samples))]
    od600 = [x for _,x in sorted(zip(chosen_days, od600))]
    pred_counts = [x for _,x in sorted(zip(chosen_days, pred_counts))]
    
    pop_counts_df = pd.DataFrame(index=chosen_samples, columns=["Prey", "Predator"])
    pop_counts_df["Prey"] = od600
    pop_counts_df["Predator"] = pred_counts

    df = df.loc[:, list(df.columns.values[:9]) + chosen_samples]

    afs = df.iloc[:, 9:].to_numpy().astype(str)
    min_2_consecutive = []

    for i in range(afs.shape[0]):
        min_2_consecutive_flag = False
        prev_val = -1
        for j in range(afs.shape[1]):
            afs[i, j] = afs[i, j].split(":")[2]
            if prev_val > min_2_consecutive_threshold and float(afs[i, j]) > min_2_consecutive_threshold:
                min_2_consecutive_flag = True
            prev_val = float(afs[i, j])

        min_2_consecutive.append(min_2_consecutive_flag)    

    df = df.loc[min_2_consecutive]

    return df, chosen_samples, pop_counts_df