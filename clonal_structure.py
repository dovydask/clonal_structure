import glob
import re
import sys
import copy
import ast
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D
from utils import *


class Node(object):
    def __init__(self, label, parent=None, time_of_emergence=None, emerged_history=None):
        self.label = label
        self.parent = parent
        self.time_of_emergence = time_of_emergence
        self.emerged_history = []
        self.children = []
        
    def add_child(self, child):
        self.children.append(child)

    def add_to_emerged_history(self, emerged):
        self.emerged_history.append(emerged)    
        

def process_pyclone_output(replicate, chosen_samples, pyclone_results_dir="./pyclone_output/"):
    
    results = pd.read_csv(pyclone_results_dir + replicate + "_final.tsv", sep="\t")
    n_samples = len(results["sample_id"].unique())
    n_clusters = len(results["cluster_id"].unique())
    
    curves = []
    for c in sorted(results["cluster_id"].unique()):
        cluster_df = results[results["cluster_id"] == c]
        mut_df = cluster_df[cluster_df["mutation_id"] == cluster_df["mutation_id"].unique()[0]]
        mut_df = mut_df.set_index("sample_id").loc[chosen_samples].reset_index()
        y = mut_df["cellular_prevalence"].astype(float)
        curves.append(y.to_numpy())
    curves = np.array(curves)
    
    return results, curves


colorblind_colors_8 = {1: "#800000",
                       2: "#000075",
                       3: "#ffe119",
                       4: "#dcbeff",
                       5: "#f58231",
                       6: "#4363d8",
                       7: "#000000",
                       8: "#a9a9a9",
                       9: "#f032e6"}

species_names_dict = {"Bd": "Brevundimonas diminuta",
                      "Ct": "Comamonas testosteroni",
                      "Pf": "Pseudomonas fluorescens",
                      "Sc": "Sphingomonas capsulata",
                      "Sm": "Serratia marcescens"}


def marginalize(tree, curves):
    adjacency_df = pd.DataFrame(columns=["Parent", "Identity", "Time of emergence"])
    for node in tree:
        if node.label != "0":
            adjacency_df.loc[len(adjacency_df), :] = [int(node.parent), int(node.label), int(node.time_of_emergence)]
    idn_to_letter = {0: "0"}
    idns = [0] + list(adjacency_df["Identity"])
    for i in range(len(idns))[1:]:
        idn_to_letter[idns[i]] = pos_to_char(i-1).upper()

    adjacency_df_letters = adjacency_df.replace({"Identity": idn_to_letter}).replace({"Parent": idn_to_letter})
    clones = []
    clones_idx = []
    for i in range(len(adjacency_df_letters)):
        p, c, t = adjacency_df_letters.iloc[i]
        pi, ci, ti = adjacency_df.iloc[i]
        clone = []
        clone_idx = []
        clone.append(adjacency_df_letters.iloc[i, 1])
        clone_idx.append(adjacency_df.iloc[i, 1])

        while p != "0":
            p = adjacency_df_letters[adjacency_df_letters["Identity"] == c]["Parent"].values[0]
            pi = adjacency_df[adjacency_df["Identity"] == ci]["Parent"].values[0]
            c = p
            ci = pi
            if p != "0":
                clone.append(p)
            if pi != 0:
                clone_idx.append(pi)

        clones.append("".join(clone[::-1]))
        clones_idx.append(clone_idx)

    order = []
    for node in tree:
        if node.label != "0":
            order.append(int(node.label))

    cluster_clone_match_df = pd.DataFrame({"Cluster": np.array(order) - 1, "Clone": clones, "Clone letter": [x[-1] for x in clones], "Time of emergence": adjacency_df_letters["Time of emergence"]})

    tdf = pd.DataFrame({"Clone": clones, "Indices": [x[::-1] for x in clones_idx]})
    tdf["Depth"] = [len(x) for x in tdf["Indices"]]
    tdf["Q"] = [curves[x[-1]-1, -1] for x in tdf["Indices"]]
    tdf["F"] = np.zeros(len(tdf))

    clone_adjacency_df = pd.DataFrame(columns=["Parent", "Identity"])
    for dp in sorted(tdf["Depth"].unique(), reverse=True):
        idx = np.where(tdf["Depth"] == dp)[0]
        for i in idx:
            c = tdf.loc[i, "Clone"]
            p = c[:-1]
            clone_adjacency_df.loc[i] = [p, c]
    clone_adjacency_df = clone_adjacency_df.sort_index()
    clone_adjacency_w_root_df = clone_adjacency_df.copy()

    for i in range(len(clone_adjacency_w_root_df)):
        if clone_adjacency_w_root_df.iloc[i, 0] == "":
            clone_adjacency_w_root_df.iloc[i, 0] = "0"

    for i in range(len(clone_adjacency_df)):
        if clone_adjacency_df.iloc[i, 0] == "":
            clone_adjacency_df.iloc[i, 0] = "0"

    clone_fractions_df = pd.DataFrame(columns=["Generation", "Identity", "Population"])
    clone_fractions_df["Generation"] = list(range(curves.shape[1]))*len(clones)
    clone_fractions_df["Identity"] = np.array([[x]*curves.shape[1] for x in clones]).reshape(1, -1)[0]
    clone_fractions_df["Population"] = np.zeros(len(clone_fractions_df))

    total_error_corrected = 0
    for t in range(curves.shape[1]):
        clones_at_time_t = cluster_clone_match_df[cluster_clone_match_df["Time of emergence"] <= t]["Clone"]
        indices_at_time_t = []
        for c in clones_at_time_t:
            clone_indices = []
            for clone_letter in c:
                clone_indices.append(cluster_clone_match_df[cluster_clone_match_df["Clone letter"] == clone_letter]["Cluster"].values[0] + 1)
            indices_at_time_t.append(clone_indices)
        tdf = pd.DataFrame({"Clone": clones_at_time_t, "Indices": indices_at_time_t})
        tdf["Depth"] = [len(x) for x in tdf["Indices"]]
        tdf["Q"] = [curves[x[-1]-1, t] for x in tdf["Indices"]]
        tdf["F"] = np.zeros(len(tdf))
        tdf = pd.concat([pd.DataFrame({"Clone": "0", "Indices": [0], "Depth": 0, "Q": 1, "F": 0}, index=[0]), tdf]).reset_index(drop = True)

        for i in range(len(tdf))[1:]:
            tdf.loc[i, "Clone"] = "0" + tdf.loc[i, "Clone"]       

        for i in range(len(tdf)):
            c = tdf.loc[i, "Clone"]
            q_sum = 0
            for j in range(i, len(tdf)):
                x = tdf.loc[j, "Clone"]
                if c == x[:len(c)] and len(c) == len(x)-1:
                    q_sum += tdf.loc[j, "Q"]   
            F = tdf.loc[i, "Q"] - q_sum
            if F < 0: F = 0
            tdf.loc[i, "F"] = F

        # MARGINALIZATION
        for i in range(len(tdf)):
            diff = -100000
            while diff < 0:
                c = tdf.loc[i, "Clone"]

                fsum = 0
                for j in range(i, len(tdf)):
                    if c in tdf.loc[j, "Clone"]: 
                        fsum += tdf.loc[j, "F"]

                diff = 1 - fsum
                if diff < 0:
                    for j in range(len(tdf)):
                        if c in tdf.loc[j, "Clone"]:
                            # tdf.loc[j, "F"] -= 0.00001
                            tdf.loc[j, "F"] -= tdf.loc[j, "F"]/1000
                            total_error_corrected += np.abs(tdf.loc[j, "F"]/1000)
                            
        for i in range(len(tdf))[1:]:
            c = tdf.loc[i, "Clone"][1:]
            clone_fractions_df.loc[clone_fractions_df[(clone_fractions_df["Identity"] == c) & (clone_fractions_df["Generation"] == t)].index.values[0], "Population"] = tdf.loc[i, "F"]

    clone_fractions_w_root_df = clone_fractions_df.copy()
    for t in range(curves.shape[1]):
        f = 1 - clone_fractions_df[clone_fractions_df["Generation"] == t]["Population"].sum()
        if f <= 0:
            f = 0
        clone_fractions_w_root_df.loc[len(clone_fractions_w_root_df)] = [t, "0", f]
    clone_fractions_w_root_df = clone_fractions_w_root_df.sort_values(["Identity", "Generation"])

    return clone_fractions_df, clone_fractions_w_root_df, clone_adjacency_df, clone_adjacency_w_root_df, cluster_clone_match_df, tdf, total_error_corrected


def compute_forest(curves, order, emerging_cutoff=0.15, min_consistency_score=0.95, max_negative_difference=-0.01):
    curves_w_root = np.vstack([np.ones(curves.shape[1])-np.max(curves, axis=0), curves])
    data = curves_w_root
    jungle = {}
    jungle_scores = pd.DataFrame(columns=["Children score", "Subclone score", "Total score", "Consistent time point score", "Tree"])
    current_fi = 0
    
    print("Begin compute trees...")

    print("Order: " + str(order[:, 0]))
    emerged_clones = ["0"]
    fi = 1
    forest = {}
    forest[str(fi)] = [Node("0", time_of_emergence=0)]

    for i in range(len(order)):
        
        em = order[i, 0]
        t = int(order[i, 1])
        time_slice = data[:, t].copy()

        if not em in emerged_clones:
            forest_keys = list(forest.keys())
            for tree_key in forest_keys:
                tree = forest[tree_key]
                if len(tree) == len(emerged_clones):
                    available_parents = [x.label for x in tree]
                    candidate_trees = [copy.deepcopy(tree) for i in range(len(available_parents))]
                    new_trees = copy.deepcopy(candidate_trees)
                    assert len(candidate_trees) == len(available_parents)
                    for ci in range(len(candidate_trees)):
                        candidate_tree = candidate_trees[ci]
                        new_tree = copy.deepcopy(new_trees[ci])
                        for ni in range(len(new_tree)):
                            node = new_tree[ni]
                            if node.label == available_parents[ci] and not em in [x.label for x in new_tree]:
                                node.add_child(em)
                                new_tree.append(Node(em, parent=available_parents[ci], time_of_emergence=t))

                        forest[str(fi)] = copy.deepcopy(new_tree)
                        fi += 1

        emerged_clones.append(em)

    print("Done.")
    print(str(len(forest.keys())) + " trees computed, scoring...", end=" ")

    complete_trees = []
    for tree_key in forest.keys():
        if len(forest[tree_key]) == data.shape[0]:
            complete_trees.append(tree_key)

    tree_scores = pd.DataFrame(index=complete_trees, columns=["Children score", "Subclone score", "Total score", "Consistent time point score"])
    tree_scores["Children score"] = 0.0
    tree_scores["Subclone score"] = 0.0
    tree_scores["Total score"] = 0.0
    tree_scores["Consistent time point score"] = 0.0
    tree_scores["Tree"] = ""

    n_all = len(complete_trees)
    print(n_all)
    ti = 1
    for tree_key in complete_trees:
        tree = forest[tree_key]
        node_scores = []
        for node in tree:
            t = node.time_of_emergence
            cf = data[int(node.label), t:]

            sum_children = np.zeros(data.shape[1]-t)
            for c in np.unique(node.children):
                sum_children += data[int(c), t:]

            cdiff = cf - sum_children

            if node.label == "0":
                sdiff = np.ones(data.shape[1]-t) - cf
            else:
                sdiff = data[int(node.parent), t:] - cf

            children_con_sum = np.sum(np.abs(cdiff[cdiff < 0]))
            subclone_con_sum = np.sum(np.abs(sdiff[sdiff < 0]))

            children_con_sum_2 = np.sum(np.array(cdiff > max_negative_difference).astype(int))
            subclone_con_sum_2 = np.sum(np.array(sdiff > max_negative_difference).astype(int))

            tree_scores.loc[tree_key, "Children score"] += children_con_sum
            tree_scores.loc[tree_key, "Subclone score"] += subclone_con_sum
            tree_scores.loc[tree_key, "Total score"] += children_con_sum + subclone_con_sum
            tree_scores.loc[tree_key, "Consistent time point score"] += children_con_sum_2 + subclone_con_sum_2
            if node.parent != None:
                tree_scores.loc[tree_key, "Tree"] += str(node.label) + "_" + str(node.parent) + ","
            else:
                tree_scores.loc[tree_key, "Tree"] += str(node.label) + "_" + "NA" + ","

        ti += 1

        if ti % 10000 == 0:
            print(ti)
        
        
    print("Done.")    
    return forest, tree_scores


def compute_and_plot_clonal_structure(replicate, tree, variant_df, pyclone_results, pyclone_curves, chosen_samples, rep_time_df, save_prefix, pyclone_results_dir, show_plot=True, day_cutoff=800, save_results=True, smoothing_std=1, xticksize=10):
    
    curves_w_root = np.vstack([np.ones(pyclone_curves.shape[1])-np.max(pyclone_curves, axis=0), pyclone_curves])
    data = curves_w_root
    clone_fractions_df, clone_fractions_w_root_df, clone_adjacency_df, clone_adjacency_w_root_df, cluster_clone_match_df, tdf, total_error_corrected = marginalize(tree, pyclone_curves)
    
    rep_days = ast.literal_eval(rep_time_df.loc[replicate, "Days"])
    rep_samples = ast.literal_eval(rep_time_df.loc[replicate, "Samples"])
    
    assert len(rep_days) == len(rep_samples)
    rep_days = [rep_days[i] for i in range(len(rep_days)) if rep_samples[i] in chosen_samples]
    
    gen_to_day_dict = {i: rep_days[i] for i in range(len(rep_days))}
    clone_fractions_w_root_df_2 = clone_fractions_w_root_df.copy()
    clone_fractions_w_root_df_2["Generation"] = [gen_to_day_dict[x] for x in clone_fractions_w_root_df["Generation"]]

    order = [int(x.label) for x in tree][1:]
    curves_post = pyclone_curves[list(np.array(order)-1)]
    marginalized_df = pd.DataFrame(index=range(curves_post.shape[1]))
    for clone in clone_fractions_df["Identity"].unique():
        x = clone[-1]
        fsums = []
        for t in range(curves_post.shape[1]):
            clone_df = clone_fractions_df[clone_fractions_df["Generation"] == t]
            f_sum = 0
            for c in clone_df["Identity"]:
                if x in c:
                    f_sum += clone_df[clone_df["Identity"] == c]["Population"].values[0]

            fsums.append(f_sum)
        marginalized_df[x] = fsums

    x = clone_fractions_w_root_df_2['Generation'].unique()
    y_table = _get_y_values(clone_fractions_w_root_df_2, clone_adjacency_w_root_df, smoothing_std)

    final_order = y_table.columns.values
    Y = y_table.to_numpy().T

    muller_colors = {'0': "#ffffff"}
    ci = 1
    assert len(order) <= 9
    for c in np.array(order):
        clone = tdf.loc[ci, "Clone"][1:]
        muller_colors[clone] = colorblind_colors_8[ci]
        ci += 1
    muller_colors = pd.Series(muller_colors)

    custom_lines = []
    spaces = " "*np.max([len(clone_adjacency_df["Identity"][i]) for i in range(len(clone_adjacency_df))])*3

    x1 = rep_days

    fig, axs = plt.subplots(2, 1, figsize=(17, 6), sharex=True)

    mutation_counts = {}
    for i in range(len(df)):
        mut_id = df.iloc[i, 0] + "_" + str(df.iloc[i, 1])
        mut_df = pyclone_results[pyclone_results["mutation_id"] == mut_id]
        if len(mut_df) != 0:
            assert np.all(mut_df["cluster_id"] == mut_df["cluster_id"].values[0])
            cluster_id = mut_df["cluster_id"].values[0]
            mut_afs = [float(x.split(":")[2]) for x in df.iloc[i].loc[chosen_samples].values]
            axs[0].plot(x1, mut_afs, color=muller_colors[cluster_clone_match_df[cluster_clone_match_df["Cluster"] == cluster_id]["Clone"].values[0]], alpha=0.3, linestyle=(0, (1, 1)))
            mutation_counts[cluster_id] = len(pyclone_results[pyclone_results["cluster_id"] == cluster_id]["mutation_id"].unique())

    # Mutation cohort/marginalized cohort plot
    ci = 0
    for i in range(len(curves_post)):
        c = int(order[i])-1
        axs[0].plot(x1, curves_post[i], color=colorblind_colors_8[ci+1], lw=2, path_effects=[pe.Stroke(linewidth=5, foreground='w'), pe.Normal()])
        custom_lines.append(Line2D([0], [0], color=colorblind_colors_8[ci+1], lw=4, label=clone_adjacency_df["Identity"][i] + spaces[len(clone_adjacency_df["Identity"][i]):] + str(mutation_counts[c])))
        axs[0].plot(x1, marginalized_df.to_numpy().T[i], color=colorblind_colors_8[ci+1], linestyle="--", lw=2)
        ci += 1

    custom_lines.append(Line2D([0], [0], color="grey", lw=2, linestyle=(0, (1, 1)), label="Variant"))
    custom_lines.append(Line2D([0], [0], color="grey", lw=2, linestyle="solid", label="Inferred cohort"))
    custom_lines.append(Line2D([0], [0], color="grey", lw=2, linestyle="dashed", label="Marginalized cohort"))

    axs[0].legend(handles=custom_lines, title_fontsize=15, loc="upper left", title="Mutation cohorts\nCohort / n variants", bbox_to_anchor=(1.05, 1.05), prop={'family': 'monospace', 'size': 15})
    axs[0].tick_params(axis='x', which='both', bottom=False, top=False, labelbottom=False, labelsize=12)
    axs[0].tick_params(axis='y', which='major', labelsize=12)
    axs[0].set_ylim(0, 1.05)
    axs[0].set_xlim(0, day_cutoff)
    axs[0].set_ylabel("Cohort frequency", fontsize=15)
    axs[0].set_title(species_names_dict[replicate.split("_")[1]] + " replicate " + str(int(replicate.split("_")[2][1:])))

    # Muller plot
    axs[1].stackplot(x.astype(np.float64), Y.astype(np.float64), colors=muller_colors.loc[final_order])
    axs[1].tick_params(axis='both', which='major', labelsize=12)
    axs[1].set_xlabel("Day in experiment", fontsize=15)
    axs[1].set_ylabel("Fraction", fontsize=15)
    axs[1].set_ylim(0, 1.05)
        
    axs[0].set_xticks(x1, x1)
    axs[1].set_xticks(x1, x1)

    axs[0].tick_params(axis="y", labelsize=15)
    axs[1].tick_params(axis="y", labelsize=15)
    axs[1].tick_params(axis="x", labelsize=xticksize, rotation=90)

    axs[0].vlines(x1, 0, 1.05, color="black", alpha=0.05)
    axs[1].vlines(x1, 0, 1.05, color="black", alpha=0.05)

    plt.margins(0)
    plt.subplots_adjust(hspace=0.1)

    if save_results:
        cluster_clone_match_df.to_csv(save_prefix + "_cluster_clone_match.tsv", sep="\t", index=False)
        clone_adjacency_w_root_df.to_csv(save_prefix + "_clone_adjacency_w_root.tsv", sep="\t", index=False)
        clone_adjacency_df.to_csv(save_prefix + "_clone_adjacency.tsv", sep="\t", index=False)
        clone_fractions_df.to_csv(save_prefix + "_clone_fractions.tsv", sep="\t", index=False)
        clone_fractions_w_root_df.to_csv(save_prefix + "_clone_fractions_w_root.tsv", sep="\t", index=False)
        plt.savefig(save_prefix + "_mutation_cohort_muller.png", bbox_inches='tight', dpi=300)
        plt.savefig(save_prefix + "_mutation_cohort_muller.svg", bbox_inches='tight', dpi=300, transparent=True)

    if show_plot:
        plt.show()
    else:
        plt.close()
        
    print("Marginalization cost: ", total_error_corrected)
    
    return clone_fractions_df, total_error_corrected


def _get_y_values(populations_df, adjacency_df, smoothing_std):

    ordering = _get_strains_ordering(adjacency_df)
    population_size_max = populations_df.groupby('Generation')['Population'].sum().max()
    generations = populations_df['Generation'].max() - populations_df['Generation'].min()

    pivot = populations_df.pivot(index='Generation', columns='Identity', values='Population').sort_index()
    pivot = pivot.rolling(generations, 1, True, 'gaussian').mean(std=smoothing_std).clip(0, population_size_max)

    Y = pivot[ordering] / 2

    # Avoid middle lines. Double leaf clones.
    keep = [0]

    for i, c in enumerate(Y.columns[1:], 1):
        if c == Y.columns[i - 1]:
            Y.iloc[:, i] *= 2
            keep.pop()

        keep.append(i)

    return Y.iloc[:, keep]


def _get_strains_ordering(adjacency_df):

    children_by_parent = adjacency_df.groupby('Parent')['Identity'].apply(lambda x: list(sorted(x)))

    def get_inner_order(identity):

        children_identities = children_by_parent.get(identity, [])

        if len(children_identities) == 0:
            return [identity, identity]

        inner = [get_inner_order(c) for c in children_identities]
        return [identity] + sum(inner, []) + [identity]

    order = []

    identities = list(set(adjacency_df['Identity'].values) | set(adjacency_df['Parent'].values))

    for strain in sorted(identities):
        if strain not in order:
            order += get_inner_order(strain)

    return np.array(order)