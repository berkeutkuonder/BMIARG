import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from MaCS_treereader import get_true_arg

def median(lst):
    """
    

    Parameters
    ----------
    lst : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return np.percentile(lst, 50, interpolation = "midpoint")

def mean(lst):
    """
    

    Parameters
    ----------
    lst : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    return np.mean(lst)

def mode(lst):
    """
    

    Parameters
    ----------
    lst : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    from scipy.stats import mode
    return float(mode(lst)[0])

def outlier_bound(lst):
    """
    

    Parameters
    ----------
    lst : TYPE
        DESCRIPTION.

    Returns
    -------
    low : TYPE
        DESCRIPTION.
    high : TYPE
        DESCRIPTION.

    """
    median = np.percentile(lst, 50, interpolation = "midpoint")
    Q1 = np.percentile(lst, 25, interpolation = "midpoint")
    Q3 = np.percentile(lst, 75, interpolation = "midpoint")
    IQR = Q3 - Q1
    if min(lst) > (Q1 - 1.96 * IQR):
        low = median - min(lst)
    else:
        low = median - (Q1 - 1.96 * IQR)
    if max(lst) < (Q3 + 1.96 * IQR):
        high = max(lst) - median
    else:
        high = (Q3 + 1.96 * IQR) - median
    return (low, high)

def confidence_interval(mu, lst):
    """
    

    Parameters
    ----------
    mu : TYPE
        DESCRIPTION.
    lst : TYPE
        DESCRIPTION.

    Returns
    -------
    low : TYPE
        DESCRIPTION.
    high : TYPE
        DESCRIPTION.

    """
    sd = np.std(lst)
    low = 1.96 * sd
    if mu - low < 0:
        low = mu
    high = 1.96 * sd
    return (low, high)

def plot_TMRCA_notfull(ARGs, seqlen, base = "median", save_data = False,
                       show_plot = True, group = 50, save_plot = False):
    """
    

    Parameters
    ----------
    ARGs : TYPE
        DESCRIPTION.
    seqlen : TYPE
        DESCRIPTION.
    base : TYPE, optional
        DESCRIPTION. The default is "median".
    save_data : TYPE, optional
        DESCRIPTION. The default is False.
    show_plot : TYPE, optional
        DESCRIPTION. The default is True.
    group : TYPE, optional
        DESCRIPTION. The default is 50.
    save_plot : TYPE, optional
        DESCRIPTION. The default is False.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    true_arg = get_true_arg()
    true_trees = true_arg.localtrees

    sites_TRUE = [true_trees[0].site]
    TMRCAs_TRUE = [true_trees[0].TMRCA]
    append_site = sites_TRUE.append
    append_TMRCA = TMRCAs_TRUE.append

    for i in range(1, len(true_trees)):
        append_site(true_trees[i].site - 1)
        append_TMRCA(true_trees[i - 1].TMRCA)
        append_site(true_trees[i].site)
        append_TMRCA(true_trees[i].TMRCA)
    append_site(seqlen) # To make sure that line goes till the end
    append_TMRCA(TMRCAs_TRUE[-1])

    TMRCA_dict = {}
    for i in range(len(ARGs)):
        for tree in ARGs[i].localtrees:
            site = tree.site
            TMRCA = tree.TMRCA
            if site in TMRCA_dict:
                TMRCA_dict[site].append(TMRCA)
            else:
                TMRCA_dict[site] = [TMRCA]

    TMRCA_dict_new = {}
    for key in TMRCA_dict.keys():
        new_key = (key - 1) // group
        if new_key in TMRCA_dict_new:
            TMRCA_dict_new[new_key] += TMRCA_dict[key]
        else:
            TMRCA_dict_new[new_key] = TMRCA_dict[key]

    sites = []
    TMRCAs = []
    low_errors = []
    high_errors = []
    append_site = sites.append
    append_TMRCA = TMRCAs.append
    append_low = low_errors.append
    append_high = high_errors.append
    for i in range(0, 10000): ## This starts from 0, be careful!! ince you used //
        if i in TMRCA_dict_new:
            if i == 0:
                append_site(1)
            else:
                append_site(i * group - group / 2)
            if base == "median":
                append_TMRCA(median(TMRCA_dict_new[i]))
                low, high = outlier_bound(TMRCA_dict_new[i])
            elif base == "mean":
                mu = mean(TMRCA_dict_new[i])
                append_TMRCA(mu)
                low,high = confidence_interval(mu, TMRCA_dict_new[i])
            elif base == "mode":
                mu = mode(TMRCA_dict_new[i])
                append_TMRCA(mu)
                low, high = confidence_interval(mu, TMRCA_dict_new[i])
            else:
                raise ValueError("Invalid base entered")
            append_low(low)
            append_high(high)

    algorithm = ARGs[0].algorithm
    errors = np.array([low_errors,high_errors])

    if show_plot:
        plt.errorbar(sites, TMRCAs,  yerr = errors, fmt='o', color="black",
                     ecolor="lightgray", elinewidth = 3, capsize = 0)
        plt.plot(sites_TRUE, TMRCAs_TRUE, 'r', zorder = 10)
        plt.title(algorithm + " - TMCRAs")
        plt.xlabel("Site")
        plt.ylabel("TMCRA")
        plt.show()

    if save_plot:
        plt.clf()
        plt.errorbar(sites, TMRCAs,  yerr = errors, fmt='o', color="black",
                 ecolor="lightgray", elinewidth = 3, capsize = 0)
        plt.plot(sites_TRUE, TMRCAs_TRUE, 'r', zorder = 10)
        plt.xlabel("Site")
        plt.ylabel("TMCRA")
        plt.tight_layout()
        plt.savefig(algorithm + "_BL.pdf")
        plt.clf()

    if save_data:
        data = np.array([sites, TMRCAs, low_errors, high_errors])
        data = np.transpose(data)
        df = pd.DataFrame(data = data, columns=["Site", "TMRCA", "low", "high"])
        df.to_csv(algorithm + "_TMRCA.csv", index = False)
        true_data = np.array([sites_TRUE, TMRCAs_TRUE])
        true_data = np.transpose(true_data)
        df = pd.DataFrame(data=true_data, columns=["Site", "TMRCA"])
        df.to_csv("True_TMRCA.csv", index = False)

def plot_TMRCA(ARGs, seqlen, base = "median", save_data = False,
               show_plot = True, group = 50, save_plot = False):
    """
    

    Parameters
    ----------
    ARGs : TYPE
        DESCRIPTION.
    seqlen : TYPE
        DESCRIPTION.
    base : TYPE, optional
        DESCRIPTION. The default is "median".
    save_data : TYPE, optional
        DESCRIPTION. The default is False.
    show_plot : TYPE, optional
        DESCRIPTION. The default is True.
    group : TYPE, optional
        DESCRIPTION. The default is 50.
    save_plot : TYPE, optional
        DESCRIPTION. The default is False.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    true_arg = get_true_arg()
    true_trees = true_arg.localtrees

    sites_TRUE = [true_trees[0].site]
    TMRCAs_TRUE = [true_trees[0].TMRCA]
    append_site = sites_TRUE.append
    append_TMRCA = TMRCAs_TRUE.append

    for i in range(1, len(true_trees)):
        append_site(true_trees[i].site - 1)
        append_TMRCA(true_trees[i - 1].TMRCA)
        append_site(true_trees[i].site)
        append_TMRCA(true_trees[i].TMRCA)
    append_site(seqlen) # To make sure that line goes till the end
    append_TMRCA(TMRCAs_TRUE[-1])

    sites = np.arange(1, seqlen)
    df = pd.DataFrame(data = sites, columns = ["Sites"])

    for i in range(len(ARGs)): #range(2565,2566):
        TMRCA_dict = {}
        for tree in ARGs[i].localtrees:
            site = tree.site
            TMRCA = tree.TMRCA
            TMRCA_dict[site] = TMRCA

        lst = [(k, v) for k, v in TMRCA_dict.items()]
        lst.sort()
        lst.append((seqlen,lst[-1][1])) # To indicate the ending
        col = np.empty([seqlen - 1])
        for j in range(len(lst) - 1):
            col[(lst[j][0] -1 ):lst[j + 1][0]] = lst[j][1]
        df["iter" + str(i)] = col

    sites = np.arange(group // 2, seqlen, group)
    TMRCAs = np.empty([seqlen // group])
    low_errors = np.empty([seqlen // group])
    high_errors = np.empty([seqlen // group])
    for i in range(seqlen // group):
        lst = df.iloc[(i * group):((i + 1) * group), 1:].to_numpy().flatten()
        if base == "median":
            TMRCAs[i] = median(lst)
            low, high = outlier_bound(lst)
        elif base == "mean":
            mu = mean(lst)
            TMRCAs[i] = mu
            low, high = confidence_interval(mu, lst)
        elif base == "mode":
            mu = mode(lst)
            TMRCAs[i] = mu
            low, high = confidence_interval(mu, lst)
        else:
            raise ValueError("Invalid base entered")
        low_errors[i] = low
        high_errors[i] = high

    algorithm = ARGs[0].algorithm
    errors = np.array([low_errors, high_errors])

    if show_plot:
        plt.errorbar(sites, TMRCAs,  yerr = errors, fmt = 'o', color = "black",
                     ecolor=  "lightgray", elinewidth = 3, capsize = 0)
        plt.plot(sites_TRUE, TMRCAs_TRUE, 'r', zorder = 10)
        plt.title(algorithm + " - TMCRAs")
        plt.xlabel("Site")
        plt.ylabel("TMCRA")
        plt.show()

    if save_plot:
        plt.clf()
        plt.errorbar(sites, TMRCAs,  yerr = errors, fmt = 'o', color = "black",
                 ecolor = "lightgray", elinewidth = 3, capsize = 0)
        plt.plot(sites_TRUE, TMRCAs_TRUE, 'r', zorder = 10)
        plt.xlabel("Site")
        plt.ylabel("TMCRA")
        plt.tight_layout()
        plt.savefig(algorithm + "_TMRCA.pdf")
        plt.clf()

    if save_data:
        data = np.array([sites, TMRCAs, low_errors,high_errors])
        data = np.transpose(data)
        df = pd.DataFrame(data=data, columns=["Site", "TMRCA", "low", "high"])
        df.to_csv(algorithm + "_TMRCA.csv", index = False)
        true_data = np.array([sites_TRUE, TMRCAs_TRUE])
        true_data = np.transpose(true_data)
        df = pd.DataFrame(data=true_data, columns=["Site", "TMRCA"])
        df.to_csv("True_TMRCA.csv", index = False)

def plot_branch_length_notfull(ARGs, seqlen, base = "median", show_plot = True,
                               save_data = False, group = 50, save_plot = False):
    """
    

    Parameters
    ----------
    ARGs : TYPE
        DESCRIPTION.
    seqlen : TYPE
        DESCRIPTION.
    base : TYPE, optional
        DESCRIPTION. The default is "median".
    show_plot : TYPE, optional
        DESCRIPTION. The default is True.
    save_data : TYPE, optional
        DESCRIPTION. The default is False.
    group : TYPE, optional
        DESCRIPTION. The default is 50.
    save_plot : TYPE, optional
        DESCRIPTION. The default is False.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """
    true_arg = get_true_arg()
    true_trees = true_arg.localtrees

    sites_TRUE = [true_trees[0].site]
    BLs_TRUE = [true_trees[0].branch_length]
    append_site = sites_TRUE.append
    append_BL = BLs_TRUE.append

    for i in range(1 ,len(true_trees)):
        append_site(true_trees[i].site - 1)
        append_BL(true_trees[i - 1].branch_length)
        append_site(true_trees[i].site)
        append_BL(true_trees[i].branch_length)
    append_site(seqlen) # To make sure that line goes till the end
    append_BL(BLs_TRUE[-1])

    BL_dict = {}
    for i in range(len(ARGs)):
        for tree in ARGs[i].localtrees:
            site = tree.site
            length = tree.branch_length
            if site in BL_dict:
                BL_dict[site].append(length)
            else:
                BL_dict[site] = [length]

    BL_dict_new = {}
    for key in BL_dict.keys():
        new_key = (key - 1) // group
        if new_key in BL_dict_new:
            BL_dict_new[new_key] += BL_dict[key]
        else:
            BL_dict_new[new_key] = BL_dict[key]

    sites = []
    BLs = []
    low_errors = []
    high_errors = []
    append_site = sites.append
    append_BL = BLs.append
    append_low = low_errors.append
    append_high = high_errors.append
    for i in range(0, seqlen): ## This starts from 0, be careful!! ince you used //
        if i in BL_dict_new:
            if i == 0:
                append_site(1)
            else:
                append_site(i * group - group / 2)
            if base == "median":
                append_BL(median(BL_dict_new[i]))
                low,high = outlier_bound(BL_dict_new[i])
            elif base == "mean":
                mu = mean(BL_dict_new[i])
                append_BL(mu)
                low,high = confidence_interval(mu, BL_dict_new[i])
            elif base == "mode":
                mu = mode(BL_dict_new[i])
                append_BL(mu)
                low,high = confidence_interval(mu, BL_dict_new[i])
            else:
                raise ValueError("Invalid base entered")
            append_low(low)
            append_high(high)

    algorithm = ARGs[0].algorithm
    errors = np.array([low_errors, high_errors])

    if show_plot:
        plt.plot(sites_TRUE, BLs_TRUE, 'r', zorder = 10)
        plt.errorbar(sites, BLs,  yerr = errors, fmt='o', color="black",
                     ecolor="lightgray", elinewidth = 3, capsize = 0)
        plt.title(algorithm + " - Branch Lengths")
        plt.xlabel("Site")
        plt.ylabel("Branch Length")
        plt.show()

    if save_plot:
        plt.clf()
        plt.plot(sites_TRUE, BLs_TRUE, 'r', zorder = 10)
        plt.errorbar(sites, BLs,  yerr = errors, fmt = 'o', color = "black",
                     ecolor="lightgray", elinewidth = 3, capsize = 0)
        plt.xlabel("Site")
        plt.ylabel("Branch Length")
        plt.tight_layout()
        plt.savefig(algorithm + "_BL.pdf")
        plt.clf()

    if save_data:
        data = np.array([sites, BLs, low_errors, high_errors])
        data = np.transpose(data)
        df = pd.DataFrame(data=data, columns=["Site", "BL", "low", "high"])
        df.to_csv(algorithm + "_BL.csv", index = False)
        true_data = np.array([sites_TRUE, BLs_TRUE])
        true_data = np.transpose(true_data)
        df = pd.DataFrame(data=true_data, columns=["Site", "BL"])
        df.to_csv("True_BL.csv", index = False)

def plot_branch_length(ARGs, seqlen, base = "median", save_data = False,
                       show_plot = True, group = 50, save_plot = False):
    """
    

    Parameters
    ----------
    ARGs : TYPE
        DESCRIPTION.
    seqlen : TYPE
        DESCRIPTION.
    base : TYPE, optional
        DESCRIPTION. The default is "median".
    save_data : TYPE, optional
        DESCRIPTION. The default is False.
    show_plot : TYPE, optional
        DESCRIPTION. The default is True.
    group : TYPE, optional
        DESCRIPTION. The default is 50.
    save_plot : TYPE, optional
        DESCRIPTION. The default is False.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """

    true_arg = get_true_arg()
    true_trees = true_arg.localtrees

    sites_TRUE = [true_trees[0].site]
    BLs_TRUE = [true_trees[0].branch_length]
    append_site = sites_TRUE.append
    append_BL = BLs_TRUE.append

    for i in range(1, len(true_trees)):
        append_site(true_trees[i].site - 1)
        append_BL(true_trees[i - 1].branch_length)
        append_site(true_trees[i].site)
        append_BL(true_trees[i].branch_length)
    append_site(seqlen) # To make sure that line goes till the end
    append_BL(BLs_TRUE[-1])

    sites = np.arange(1, seqlen)
    df = pd.DataFrame(data = sites, columns=["Sites"])

    for i in range(len(ARGs)): #range(2565,2566):
        BL_dict = {}
        for tree in ARGs[i].localtrees:
            site = tree.site
            branch_length = tree.branch_length
            BL_dict[site] = branch_length

        lst = [(k, v) for k, v in BL_dict.items()]
        lst.sort()
        lst.append((seqlen,lst[-1][1])) # To indicate the ending
        col = np.empty([seqlen - 1])
        for j in range(len(lst) - 1):
            col[(lst[j][0] - 1):lst[j + 1][0]] = lst[j][1]
        df["iter"+str(i)] = col

    sites = np.arange(group // 2, seqlen, group)
    BLs = np.empty([seqlen//group])
    low_errors = np.empty([seqlen // group])
    high_errors = np.empty([seqlen // group])
    for i in range(seqlen // group):
        lst = df.iloc[(i * group):((i + 1) * group), 1:].to_numpy().flatten()
        if base == "median":
            BLs[i] = median(lst)
            low, high = outlier_bound(lst)
        elif base == "mean":
            mu = mean(lst)
            BLs[i] = mu
            low,high = confidence_interval(mu, lst)
        elif base == "mode":
            mu = mode(lst)
            BLs[i] = mu
            low,high = confidence_interval(mu, lst)
        else:
            raise ValueError("Invalid base entered")
        low_errors[i] = low
        high_errors[i] = high


    algorithm = ARGs[0].algorithm
    errors = np.array([low_errors, high_errors])

    if show_plot:
        plt.plot(sites_TRUE, BLs_TRUE, 'r', zorder = 10)
        plt.errorbar(sites, BLs,  yerr = errors, fmt = 'o', color="black",
                     ecolor="lightgray", elinewidth = 3, capsize = 0)
        plt.title(algorithm + ' - Branch Lengths')
        plt.xlabel("Site")
        plt.ylabel("Branch Length")
        plt.show()

    if save_plot:
        plt.clf()
        plt.plot(sites_TRUE, BLs_TRUE, 'r', zorder = 10)
        plt.errorbar(sites, BLs,  yerr = errors, fmt = 'o', color = "black",
                     ecolor = "lightgray", elinewidth = 3, capsize = 0)
        plt.xlabel("Site")
        plt.ylabel("Branch Length")
        plt.tight_layout()
        plt.savefig(algorithm + "_BL.pdf")
        plt.clf()

    if save_data:
        data = np.array([sites, BLs, low_errors, high_errors])
        data = np.transpose(data)
        df = pd.DataFrame(data=data, columns=["Site", "BL", "low", "high"])
        df.to_csv(algorithm + "_BL.csv", index = False)
        true_data = np.array([sites_TRUE, BLs_TRUE])
        true_data = np.transpose(true_data)
        df = pd.DataFrame(data=true_data, columns=["Site", "BL"])
        df.to_csv("True_BL.csv", index = False)

def plot_frequency_sites(ARGs, seqlen, save_data = False, group = 20,
                         show_plot = True, save_plot = False):
    """
    

    Parameters
    ----------
    ARGs : TYPE
        DESCRIPTION.
    seqlen : TYPE
        DESCRIPTION.
    save_data : TYPE, optional
        DESCRIPTION. The default is False.
    group : TYPE, optional
        DESCRIPTION. The default is 20.
    show_plot : TYPE, optional
        DESCRIPTION. The default is True.
    save_plot : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """

    true_arg = get_true_arg()
    TRUE_sites = [tree.site for tree in true_arg.localtrees][1:]

    sites = []
    append_site = sites.append
    for i in range(len(ARGs)):
        for tree in ARGs[i].localtrees[1:]:
            append_site(tree.site)

    sites = (np.array(sites) // group) * group + (group // 2)

    df = pd.DataFrame(data = sites, columns=["Sites"])
    a = df.value_counts(normalize = True)
    a = a.reset_index()
    a = a.append({"Sites": 0, 0:0}, ignore_index = True)
    a = a.append({"Sites": seqlen, 0:0}, ignore_index = True)
    for site in range(group // 2, seqlen, group):
        if site in a["Sites"].values:
            continue
        a = a.append({"Sites": site, 0:0}, ignore_index = True)
    a = a.sort_values(by = ["Sites"])

    algorithm = ARGs[0].algorithm

    if show_plot:
        plt.plot(a["Sites"], a[0], 'k', zorder = 0)
        for site in TRUE_sites:
            plt.axvline(x = site, ymin = 0.92, color = 'r', linewidth = 0.75)
        plt.title(algorithm + " - Recombination Location")
        plt.xlabel("Site")
        plt.ylabel("Frequency")
        plt.show()

    if save_plot:
        plt.clf()
        plt.plot(a["Sites"], a[0], 'k', zorder = 0)
        for site in TRUE_sites:
            plt.axvline(x=site, ymin = 0.92, color = 'r', linewidth = 0.75)
        plt.xlabel("Site")
        plt.ylabel("Frequency")
        plt.tight_layout()
        plt.savefig(algorithm + "_FreqSites.pdf")
        plt.clf()

    if save_data:
        algorithm = ARGs[0].algorithm
        a.to_csv(algorithm + "_FreqSites.csv", index = False)
        true_data = np.array([TRUE_sites])
        true_data = np.transpose(true_data)
        df = pd.DataFrame(data=true_data, columns=["Site"])
        df.to_csv("True_FreqSites.csv", index = False)

def plot_RF_dist(ARGs, seqlen, base = "median", save_data = False, group = 50,
                 show_plot = True, save_plot = False):
    """
    

    Parameters
    ----------
    ARGs : TYPE
        DESCRIPTION.
    seqlen : TYPE
        DESCRIPTION.
    base : TYPE, optional
        DESCRIPTION. The default is "median".
    save_data : TYPE, optional
        DESCRIPTION. The default is False.
    group : TYPE, optional
        DESCRIPTION. The default is 50.
    show_plot : TYPE, optional
        DESCRIPTION. The default is True.
    save_plot : TYPE, optional
        DESCRIPTION. The default is False.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    None.

    """

    true_arg = get_true_arg()
    true_trees = true_arg.localtrees

    sites = np.arange(1, seqlen)
    df = pd.DataFrame(data = sites, columns = ["Sites"])

    for i in range(len(ARGs)): # range(2565,2566):
        RF_dict = {}
        for tree in ARGs[i].localtrees:
            site = tree.site
            for j in range(len(true_trees)):
                if true_trees[j].site > site:
                    dist = tree.robinson_foulds(true_trees[j - 1])
                    break
            else:
                dist = tree.robinson_foulds(true_trees[-1])

            RF_dict[site] = dist
        for tree in true_trees:
            site = tree.site
            trees = ARGs[i].localtrees
            for j in range(len(trees)):
                if trees[j].site > site:
                    dist = tree.robinson_foulds(trees[j - 1])
                    break
            else:
                dist = tree.robinson_foulds(trees[-1])
            RF_dict[site] = dist
        lst = [(k, v) for k, v in RF_dict.items()]
        lst.sort()
        lst.append((seqlen,lst[-1][1])) # To indicate the ending
        col = np.empty([seqlen - 1])
        for j in range(len(lst) - 1):
            col[(lst[j][0] - 1):lst[j + 1][0]] = lst[j][1]
        df["iter"+str(i)] = col

    sites = np.arange(group // 2, seqlen, group)
    RFs = np.empty([seqlen // group])
    low_errors = np.empty([seqlen // group])
    high_errors = np.empty([seqlen // group])
    for i in range(seqlen // group):
        lst = df.iloc[(i * group):((i + 1) * group), 1:].to_numpy().flatten()
        if base == "median":
            RFs[i] = median(lst)
            low, high = outlier_bound(lst)
        elif base == "mean":
            mu = mean(lst)
            RFs[i] = mu
            low,high = confidence_interval(mu, lst)
        elif base == "mode":
            mu = mode(lst)
            RFs[i] = mu
            low,high = confidence_interval(mu, lst)
        else:
            raise ValueError("Invalid base entered")
        low_errors[i] = low
        high_errors[i] = high

    algorithm = ARGs[0].algorithm
    errors = np.array([low_errors, high_errors])

    if show_plot:
        plt.errorbar(sites, RFs,  yerr = errors, fmt = 'o', color = "black",
                     ecolor="lightgray", elinewidth = 3, capsize = 0)
        plt.title(algorithm + " - RF Distance")
        plt.xlabel("Site")
        plt.ylabel("RF Distance")
        plt.show()

    if save_plot:
        plt.clf()
        plt.errorbar(sites, RFs,  yerr = errors, fmt = 'o', color = "black",
                 ecolor="lightgray", elinewidth = 3, capsize = 0)
        plt.xlabel("Site")
        plt.ylabel("RF Distance")
        plt.tight_layout()
        plt.savefig(algorithm + "_RF.pdf")
        plt.clf()

    if save_data:
        data = np.array([sites, RFs, low_errors,high_errors])
        data = np.transpose(data)
        df = pd.DataFrame(data = data, columns = ["Site", "RF", "low", "high"])
        df.to_csv(algorithm + "_RF.csv", index = False)


#%%




