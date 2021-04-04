import os
from utils import localtree, ARG, find_recombination
from plotting import plot_TMRCA, plot_branch_length, plot_frequency_sites, plot_RF_dist
from time import time
start = time()
print("Arbores Analysis has started!")

def read_tex_file(filename, remove_zero = False):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    remove_zero : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    f = open(filename)
    inp = f.read()
    f.close()

    tex_trees = inp.split("\input{./phytree.tex}")[:-1] # Last item is not a tree

    localtrees = []
    append = localtrees.append
    for tex in tex_trees:
        tr = read_tex_tree(tex, remove_zero)
        append(tr)

    find_recombination(localtrees)
    return ARG(None, localtrees, "Arbores")

def read_tex_tree(tex, remove_zero = False):
    """
    

    Parameters
    ----------
    tex : TYPE
        DESCRIPTION.
    remove_zero : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    site = int(tex.split("\def\site{")[1].split('}')[0])

    first_nodes = tex.split("\def\pone{{")[1].split("}}")[0]
    first_nodes = [int(i) for i in first_nodes.split(',')]

    second_nodes = tex.split("\def\ptwo{{")[1].split("}}")[0]
    second_nodes = [int(i) for i in second_nodes.split(',')]

    times = tex.split("\def\T{{")[1].split("}}")[0]
    times = [float(i) for i in times.split(',')]

    nseq = int(tex.split("\def\\NL{")[1].split('}')[0])

    third_nodes = list(range(nseq, 2 * nseq - 1))

    nodes = [[first_nodes[i], second_nodes[i], third_nodes[i]] for i in range(len(third_nodes))]

    if not remove_zero:
        return localtree(site, times, nodes)

    #Removing node 0
    for node in nodes:
        if 0 in node: # Lost coalescence when 0 is removed.
            replacer = node[1]
            lost_index = node[2]

    # Replace the lost_index
    for i in range(len(nodes)):
        node = nodes[i]
        if lost_index in node and lost_index != node[2]:
            if node[0] == lost_index:
                nodes[i] = [replacer, node[1], node[2]]
            elif replacer < node[0]:
                nodes[i] = [replacer, node[0], node[2]]
            else:
                nodes[i] = [node[0], replacer, node[2]]
            break

    # Delete lost index
    for i in range(len(nodes)):
        node = nodes[i]
        if lost_index in node:
            del nodes[i]
            break

    # Adjust number of other nodes
    for i in range(len(nodes)):
        for j in range(3):
            if nodes[i][j] > lost_index:
                nodes[i][j] -= 2
            else:
                nodes[i][j] -= 1

    # Deletet times
    del times[lost_index],times[0]

    return localtree(site,times,nodes)

def read_chain(chain, remove_zero = False): # Reads one chains

    items = [float(i) for i in chain.split(' ')[:-1]]

    nseq = int(items[0])

    info_amount = 3 * (nseq - 1) + 5

    if len(items) % info_amount != 0:
        raise NameError("Something is wrong with the file, chain.txt")

    trees_amount = len(items) // info_amount

    localtrees = []
    for i in range(trees_amount):
        tree = items[(i * info_amount):((i + 1) * info_amount)]

        first_nodes = [int(n) for n in tree[1:(nseq)]]
        second_nodes = [int(n) for n in tree[(nseq):(2 * nseq - 1)]]
        third_nodes = list(range(nseq, 2 * nseq - 1))

        nodes = [[first_nodes[i], second_nodes[i], third_nodes[i]] for i in range(len(third_nodes))]
        times = [0] * nseq + tree[(2 * nseq - 1):(3 * nseq - 2)]
        site = int(tree[(3 * nseq - 2)])

        if not remove_zero:
            localtrees.append(localtree(site, times, nodes))
            continue

        #Removing node 0
        for node in nodes:
            if 0 in node: # Lost coalescence when 0 is removed.
                replacer = node[1]
                lost_index = node[2]

        # Replace the lost_index
        for i in range(len(nodes)):
            node = nodes[i]
            if lost_index in node and lost_index != node[2]:
                if node[0] == lost_index:
                    nodes[i] = [replacer, node[1], node[2]]
                elif replacer < node[0]:
                    nodes[i] = [replacer, node[0], node[2]]
                else:
                    nodes[i] = [node[0], replacer, node[2]]
                break

        # Delete lost index
        for i in range(len(nodes)):
            node = nodes[i]
            if lost_index in node:
                del nodes[i]
                break

        # Adjust number of other nodes
        for i in range(len(nodes)):
            for j in range(3):
                if nodes[i][j] > lost_index:
                    nodes[i][j] -= 2
                else:
                    nodes[i][j] -= 1

        # Deletet times
        del times[lost_index],times[0]
        localtrees.append(localtree(site, times, nodes))

    find_recombination(localtrees)
    return ARG(None, localtrees, "Arbores")

def read_chain_file(filename, remove_zero = False):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.
    remove_zero : TYPE, optional
        DESCRIPTION. The default is False.

    Returns
    -------
    ARGs : TYPE
        DESCRIPTION.

    """
    f = open(filename)
    inp = f.read()
    f.close()

    iters = inp.split('\n')[:-1]
    ARGs = []
    append = ARGs.append
    for chain in iters:
        append(read_chain(chain, remove_zero))
    return ARGs

def Arbores_read_diagnostics(filename):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    nrecombs : TYPE
        DESCRIPTION.
    posterior : TYPE
        DESCRIPTION.

    """
    f = open(filename)
    inp = f.read()
    f.close()

    nrecombs = [int(i.split()[10]) for i in inp.split('\n')[:-1]] #Recombination is stored in 11th column
    posterior = [float(i.split()[15]) for i in inp.split('\n')[:-1]] #Posterior is stored in 16th column

    return nrecombs, posterior

def Arbores_recombs(nrecombs, show_plot = True, save_plot = False, save_data = False):
    """
    

    Parameters
    ----------
    nrecombs : TYPE
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

    if show_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(nrecombs)), nrecombs, 'o')
        plt.title("Arbores Recombination per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("# recombinations")
        plt.show()

        plt.hist(nrecombs, density = True)
        plt.title("Arbores Recombinations")
        plt.xlabel("# recombinations")
        plt.ylabel("Density")
        plt.show()

        plt.plot(range(len(nrecombs)), nrecombs)
        plt.title("Arbores Recombination per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("# recombinations")
        plt.show()
        plt.clf()

    if save_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(nrecombs)), nrecombs, 'o')
        plt.title("Arbores Recombination per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("# recombinations")
        plt.savefig("Arbores_recombs_scatter.pdf")
        plt.clf()

        plt.hist(nrecombs, density = True)
        plt.title("Arbores Recombinations")
        plt.xlabel("# recombinations")
        plt.ylabel("Density")
        plt.savefig("Arbores_recombs_density.pdf")
        plt.clf()

        plt.plot(range(len(nrecombs)), nrecombs)
        plt.title("Arbores Recombination per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("# recombinations")
        plt.savefig("Arbores_recombs_line.pdf")
        plt.clf()

    if save_data:
        f = open("Arbores_recombs.txt", 'w')
        f.write("iter, nrecombs\n")
        for i in range(len(nrecombs)):
            line = str(i) + ", " + str(nrecombs[i]) + '\n'
            f.write(line)
        f.close()

def Arbores_posterior(posterior, show_plot = True, save_plot = False, save_data = False):

    if show_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(posterior) - 2), posterior[2:])
        plt.title("Arbores Log Posterior per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("Log Posterior")
        plt.show()
        plt.clf()

    if save_plot:
        import matplotlib.pyplot as plt
        plt.plot(range(len(posterior) - 2), posterior[2:])
        plt.title("Arbores Log Posterior per iteration")
        plt.xlabel("Iteration")
        plt.ylabel("Log Posterior")
        plt.savefig("Arbores_posterior_line.pdf")
        plt.clf()

    if save_data:
        f = open("Arbores_posterior.txt", 'w')
        f.write("iter, logposterior\n")
        for i in range(len(posterior)):
            line = str(i) + ", " + str(posterior[i]) + '\n'
            f.write(line)
        f.close()

if not os.getcwd().endswith("Arbores_1"):
    os.chdir("Arbores_1")

print("Started reading the chain file!")
ARGs = read_chain_file("chain.txt", remove_zero = True)
print("Ended reading the chain file!")

if not os.path.exists("plots"):
    print("Creating 'plots' directory")
    os.makedirs("plots")

os.chdir("plots")
seqlen = 10000

nrecombs, posterior = Arbores_read_diagnostics("../diagnostics.txt")
Arbores_recombs(nrecombs, show_plot = False, save_data = True)
print("Recombination information is saved")
Arbores_posterior(posterior, show_plot = False, save_data = True)
print("Posterior information is saved")

plot_TMRCA(ARGs, seqlen, base = "median", save_data = True, show_plot = False, save_plot = False)
print("TMRCA plot is created")
plot_branch_length(ARGs, seqlen, base = "median", save_data = True, show_plot = False, save_plot = False)
print("Branch length plot is created")
plot_frequency_sites(ARGs, seqlen, save_data = True, show_plot = False, save_plot = False)
print("Frequency sites plot is created")
plot_RF_dist(ARGs, seqlen, base = "median", save_data = True, show_plot = False, save_plot = False)
print("RF distance plot is created")

ARG_map = read_tex_file("../map.tex", remove_zero = True)
ARG_map.print_arg(dir_name = "map")
print("ARG_map has is created")

ARG_init = ARGs[0]
ARG_init.print_arg(dir_name = "init")
print("ARG_init has is created")

ARG_final = ARGs[-1]
ARG_final.print_arg(dir_name = "final")
print("ARG_final has is created")

end = time()
print("Arbores Analysis completed in", (end-start) / 60, "minutes")

