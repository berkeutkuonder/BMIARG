import os
import pandas as pd
from utils import find_matching_par, localtree, find_recombination, ARG
from plotting import plot_TMRCA, plot_branch_length, plot_frequency_sites, plot_RF_dist
from time import time
print("ARGweaver Analysis has started!")
start = time()

def ARGweaver_read_stats(filename):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    df : TYPE
        DESCRIPTION.

    """
    df = pd.read_csv(filename, sep = '\t')
    return df

def ARGweaver_recombs(df, show_plot = True, save_plot = False, save_data = False):
    """
    

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
    show_plot : TYPE, optional
        DESCRIPTION. The default is True.
    save_plot : TYPE, optional
        DESCRIPTION. The default is False.
    save_data : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """

    nrecombs = df["recombs"]

    if show_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(nrecombs)), nrecombs, 'o')
        plt.title("ARGweaver Recombination per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("# recombinations")
        plt.show()

        plt.hist(nrecombs, density=True)
        plt.title("ARGweaver Recombinations")
        plt.xlabel("# recombinations")
        plt.ylabel("Density")
        plt.show()

        plt.plot(range(len(nrecombs)), nrecombs)
        plt.title("ARGweaver Recombination per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("# recombinations")
        plt.show()
        plt.clf()

    if save_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(nrecombs)), nrecombs, 'o')
        plt.title("ARGweaver Recombination per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("# recombinations")
        plt.savefig("ARGweaver_recombs_scatter.pdf")
        plt.clf()

        plt.hist(nrecombs, density = True)
        plt.title("ARGweaver Recombinations")
        plt.xlabel("# recombinations")
        plt.ylabel("Density")
        plt.savefig("ARGweaver_recombs_density.pdf")
        plt.clf()

        plt.plot(range(len(nrecombs)), nrecombs)
        plt.title("ARGweaver Recombination per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("# recombinations")
        plt.savefig("ARGweaver_recombs_line.pdf")
        plt.clf()

    if save_data:
        f = open("ARGweaver_recombs.txt", 'w')
        f.write("iter, nrecombs\n")
        for i in range(len(nrecombs)):
            line = str(i) + ", " + str(nrecombs[i]) + '\n'
            f.write(line)
        f.close()

def ARGweaver_posterior(df, show_plot = True, save_plot = False, save_data = False):
    """
    

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
    show_plot : TYPE, optional
        DESCRIPTION. The default is True.
    save_plot : TYPE, optional
        DESCRIPTION. The default is False.
    save_data : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """

    posterior = df["joint"]

    if show_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(posterior) - 2), posterior[2:])
        plt.title("ARGweaver Log Posterior per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("Log Posterior")
        plt.show()
        plt.clf()

    if save_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(posterior) - 2), posterior[2:])
        plt.title("ARGweaver Log Posterior per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("Log Posterior")
        plt.savefig("ARGweaver_posterior_line.pdf")
        plt.clf()

    if save_data:
        f = open("ARGweaver_posterior.txt", 'w')
        f.write("iter, logposterior\n")
        for i in range(len(posterior)):
            line = str(i) + ", " + str(posterior[i]) + "\n"
            f.write(line)
        f.close()

def ARGweaver_prior(df, show_plot = True, save_plot = False, save_data = False):
    """
    

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
    show_plot : TYPE, optional
        DESCRIPTION. The default is True.
    save_plot : TYPE, optional
        DESCRIPTION. The default is False.
    save_data : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """

    prior = df["prior"]

    if show_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(prior) - 2), prior[2:])
        plt.title("ARGweaver Log prior per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("Log prior")
        plt.show()
        plt.clf()

    if save_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(prior) - 2), prior[2:])
        plt.title("ARGweaver Log prior per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("Log prior")
        plt.savefig("ARGweaver_prior_line.pdf")
        plt.clf()

    if save_data:
        f = open("ARGweaver_prior.txt", 'w')
        f.write("iter, logprior\n")
        for i in range(len(prior)):
            line = str(i) + ", " + str(prior[i]) + "\n"
            f.write(line)
        f.close()

def ARGweaver_likelihood(df, show_plot = True, save_plot = False, save_data = False):
    """
    

    Parameters
    ----------
    df : TYPE
        DESCRIPTION.
    show_plot : TYPE, optional
        DESCRIPTION. The default is True.
    save_plot : TYPE, optional
        DESCRIPTION. The default is False.
    save_data : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    None.

    """

    likelihood = df["likelihood"]

    if show_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(likelihood) - 2), likelihood[2:])
        plt.title("ARGweaver Log likelihood per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("Log likelihood")
        plt.show()
        plt.clf()

    if save_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(likelihood) - 2), likelihood[2:])
        plt.title("ARGweaver Log likelihood per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("Log likelihood")
        plt.savefig("ARGweaver_likelihood_line.pdf")
        plt.clf()

    if save_data:
        f = open("ARGweaver_likelihood.txt", 'w')
        f.write("iter, loglikelihood\n")
        for i in range(len(likelihood)):
            line = str(i) + ", " + str(likelihood[i]) + "\n"
            f.write(line)
        f.close()

def ARGweaver_read_tree(nseq, tree, individuals):
    """
    Recursive function - Reads argweaver tree and returns the list of nodes

    Parameters
    ----------
    nseq: int
        Number of individuals (number of leaves in the tree)
    tree : string
        the string version of tree
    individuals: int list
        ARGweaver presents individuals as n0,n1,n2 etc but n0 does not
        necessarily indicate that n0 is at position 0. This list contains the
        correct positions of individuals (n0,n1,n2)

    Returns
    -------
    int * list list
        the node number of the current tmcra, nodes list

    """
    if tree[0] == '(': # If tree
        close, _ = find_matching_par('(', tree, 0) # Closing index of par + 1
        j = close
        while tree[j] != ']':
            j += 1
        temp = tree[close:j].split("[&&NHX:age=")
        if (len(tree) - 1) > j: # Needed this, so that sorting gives this the last
            # print(j)
            if tree[j + 1] == ';':
                pos = 150
        else:
            pos = int(temp[0].split(':')[0])
        age = float(temp[1])
    else: # If leaf
        temp = tree[:-1].split("[&&NHX:age=")
        pos = int(temp[0].split(':')[0])
        if pos < nseq:
            pos = individuals[pos]
        age = float(temp[1])
        return pos, []

    if tree[1] == '(':
        i, _ = find_matching_par('(', tree, 1)
        while tree[i] != ',': # , represents the seperation of the left and right
            i += 1
        subtree1 = tree[1:i]
        subtree2 = tree[(i + 1):(close - 1)]
        pos1, nodes1 = ARGweaver_read_tree(nseq, subtree1, individuals)
        pos2, nodes2 = ARGweaver_read_tree(nseq, subtree2, individuals)
        nodes = nodes1 + nodes2 + [[age, pos, pos1, pos2]]
        nodes.sort() # Sorting based on age
        for i in range(len(nodes) - 1):
            if nodes[i][0] == nodes[i + 1][0]:
                # Finding the ancestor
                if nodes[i][1] not in nodes[i + 1]:
                    # Swaping them based on ancestry
                    nodes[i], nodes[i + 1] = nodes[i + 1], nodes[i]
        return pos, nodes
    else:
        i = 2
        while tree[i] != ',': # , represents the seperation of the left and right
            i += 1
        subtree1 = tree[1:i]
        subtree2 = tree[(i + 1):(close - 1)]
        pos1, nodes1 = ARGweaver_read_tree(nseq, subtree1, individuals)
        pos2, nodes2 = ARGweaver_read_tree(nseq, subtree2, individuals)
        nodes = nodes1 + nodes2 + [[age, pos, pos1, pos2]]
        nodes.sort() # Sorting based on age
        for i in range(len(nodes) - 1):
            if nodes[i][0] == nodes[i + 1][0]:
                # Finding the ancestor
                if nodes[i][1] not in nodes[i + 1]:
                    # Swaping them based on ancestry
                    nodes[i], nodes[i + 1] = nodes[i + 1], nodes[i]
        return pos, nodes

def ARGweaver_adjust_nodes(nodes, nseq):
    """
    Argweaver node numbering, by default, does not follow a chronological order.
    This function converts node number in a chronological order.

    Parameters
    ----------
    nodes : list list
        List of nodes read from the raw input
    nseq : int
        The number of sequences (number of leaves) in the tree

    Returns
    -------
    None.

    """
    most_recent = nseq
    for node in nodes: # String node indicates the edited node
        n = node[-3]
        for item in nodes:
            if n == item[-1]:
                item[-1] = str(most_recent)
                break
            elif n == item[-2]:
                item[-2] = str(most_recent)
                break
        node[-3] = str(most_recent)
        most_recent += 1
    for node in nodes: # Making all nodes int again
        node[-1] = int(node[-1])
        node[-2] = int(node[-2])
        node[-3] = int(node[-3])

def ARGweaver_create_localtree(site,nodes,nseq,pop):
    """ Creates a local tree object given nodes.

    Parameters
    ----------
    site : int
        Site of the tree
    nodes : list list
        list of nodes in format of [t, c, a, b] where a and b coalesces at c
        at time t
    nseq : int
        The number of sequences (number of leaves) in the tree
    pop: int
        The effective population (this is needed to adjust the coalescence times)

    Returns
    -------
    localtree
        Created local tree object

    """
    times = [0] * nseq
    final_nodes = []
    time_adjust = 4 * pop # Branch lengths are generally shows as 4N
    for node in nodes:
        time = round(node[0] / time_adjust, 6)
        times.append(time)
        if node[-1] > node [-2]:
            new_node = (node[-2], node[-1], node[-3])
        else:
            new_node = (node[-1], node[-2], node[-3])
        final_nodes.append(new_node)
    return localtree(site,times, final_nodes, anc_map = False)

def ARGweaver_read_smc_file(filename):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    f = open(filename)
    inp = f.read()
    f.close()
    lines = inp.split("\n")
    ## This is needed to find the correct positions of the sequences
    individuals = [int(i[1:]) for i in lines[0].split("\t")[1:]]
    nseq = len(individuals)
    trees = []
    for i in range(2, len(lines), 2):
        temp = lines[i].split("\t")
        site = int(temp[1])
        tree = temp[3]
        _ , nodes = ARGweaver_read_tree(nseq, tree, individuals)
        ARGweaver_adjust_nodes(nodes, nseq)
        t = ARGweaver_create_localtree(site, nodes, nseq, pop)
        trees.append(t)
    print("Starting to add recombs!")
    find_recombination(trees)
    # print("Recombinations added successfully!")
    # Checking for invisible recombinations (only applies to ARGweaver)
    for i in range(len(trees) - 1): # Excluding the last tree
        if not trees[i].recomb_exists:
            # print("Adding special recomb at tree", i)
            tree = trees[i]
            site = str(tree.site)
            for j in range(len(lines)):
                if lines[j][5:(5 + len(site))] == site:
                    line = lines[j+1].split("\t")
                    recomb_to_time = round((float(line[-1]) / (4 * pop)), 6)
                    recomb_to = trees[i + 1].times.index(recomb_to_time)
                    recomb_from = line[-4]
                    if int(recomb_from) < nseq:
                        recomb_from_time = (0 + recomb_to_time) / 2
                        trees[i].add_recomb(int(recomb_from), recomb_to, recomb_from_time, recomb_to_time)
                        break
                    k = 0
                    while k < len(lines[j]):
                        index = k + len(recomb_from) + 1
                        if lines[j][k:index] == recomb_from + ':':
                            if lines[j][k - 1] == '(' or lines[j][k - 1] == ')' or \
                            lines[j][k - 1] == ',':
                                time = lines[j][k:].split("[&&NHX:age=")[1]
                                recomb_from_ind_time = round((float(time.split(']')[0]) / (4 * pop)), 6)
                                recomb_from = tree.times.index(recomb_from_ind_time)
                                recomb_from_time = (recomb_from_ind_time + recomb_to_time) / 2
                                trees[i].add_recomb(recomb_from, recomb_to, recomb_from_time, recomb_to_time)
                                break
                        k += 1
                    else:
                        print("Something unexpectd happend at tree", i)
                        print("Filename:", filename)
                        raise ValueError
                    break

    return ARG(None, trees, "ARGweaver")

def ARGweaver_get_parameters(filename):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    pop : TYPE
        DESCRIPTION.
    seqlen : TYPE
        DESCRIPTION.
    max_iter : TYPE
        DESCRIPTION.

    """
    f = open(filename)
    inp = f.read()
    f.close()
    lines = inp.split("\n")
    pop = int(lines[2].split("-N ")[1].split(" -r")[0])
    seqlen = int(lines[4].split("length=")[1].split(',')[0])
    max_iter = int(lines[2].split("-n ")[1].split(" --")[0])
    return pop, seqlen, max_iter
#%%
if not os.getcwd().endswith("ARGweaver"):
    os.chdir("ARGweaver")

pop, seqlen, max_iter = ARGweaver_get_parameters("out.log")

for filename in os.listdir(): # unzipping smc files (if any)
    if filename.endswith(".smc.gz"):
        os.system("gzip -d " + filename)

ARGs = []
append = ARGs.append
for i in range(max_iter + 1) :
    filename = "out." + str(i) + ".smc"
    print("Working on", filename)
    append(ARGweaver_read_smc_file(filename))

if not os.path.exists("plots"):
    print("Creating 'plots' directory")
    os.makedirs("plots")

os.chdir("plots")

plot_TMRCA(ARGs, seqlen, base = "median", save_data = True, show_plot = False, save_plot = True)
print("TMRCA plot is created")
plot_branch_length(ARGs, seqlen, base = "median", save_data = True, show_plot = False, save_plot = True)
print("Branch length plot is created")
plot_frequency_sites(ARGs, seqlen, save_data = True, show_plot = False, save_plot = True)
print("Frequency sites plot is created")
plot_RF_dist(ARGs, seqlen, base = "median", save_data = True, show_plot = False, save_plot = True)
print("RF distance plot is created")

df = ARGweaver_read_stats("../out.stats")
ARGweaver_recombs(df, show_plot = False, save_data = True)
print("Recombination information is saved")
ARGweaver_posterior(df, show_plot = False, save_data = True)
print("Posterior information is saved")
# ARGweaver_prior(df, show_plot = True, save_data = False)
# ARGweaver_likelihood(df, show_plot = True, save_data = False)

index = df["joint"].idxmax()
ARG_map = ARGs[index - 1]
ARG_map.print_arg(dir_name = "map")
print("ARG_map has is created")

ARG_init = ARGs[0]
ARG_init.print_arg(dir_name = "init")
print("ARG_init has is created")

ARG_final = ARGs[-1]
ARG_final.print_arg(dir_name = "final")
print("ARG_final has is created")

end = time()
print("ARGweaver Analysis completed in", (end - start) / 60, "minutes")

#%%
pop, seqlen=10000, 10000
a = ARGweaver_read_smc_file("out.10000.smc")
a.print_info()
#%%
a.localtrees[2].print_tree(keep_latex = True)