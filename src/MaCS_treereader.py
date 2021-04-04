import os, sys

from utils import find_matching_par, localtree, find_recombination, ARG

def adjust_newick(nseq, tree):
    """
    Adjusts newick tree by labeling all nodes

    Parameters
    ----------
    nseq : int
        number of sequences (leaves) in the tree
    tree : string
        the string version of tree

    Returns
    -------
    string
        edited newick tree
    """
    index = nseq
    final = []
    append = final.append
    i = 0
    while tree[i] != ';':
        if tree[i] == ':':
            if tree[i - 1] == ')':
                append(str(index))
                index += 1
        append(tree[i])
        i += 1
    append(str(index) + ":0")
    return ''.join(final)

def read_tree(tree):
    """
    Recursive function - Reads EDITED newick tree and returns the list of nodes

    Parameters
    ----------
    tree : string
        the string version of EDITED newick tree (by adjust_newick function)

    Returns
    -------
    int * float * list
        the node number of the current tmcra * TMCRA time * nodes list

    """
    if tree[0] == '(': # If tree
        close, _ = find_matching_par('(', tree, 0) # Closing index of par + 1
        j = close
        while tree[j] != ')':
            j += 1
            if j == len(tree):
                break
        temp = tree[close:j].split(":")
        pos = int(temp[0])
        age = float(temp[1])
    else: # If leaf
        temp = tree.split(":")
        pos = int(temp[0])
        age = float(temp[1])
        return pos, age, []

    if tree[1] == '(':
        i, _ = find_matching_par('(', tree, 1)
        while tree[i] != ',': # , represents the seperation of the left and right
            i += 1
        subtree1 = tree[1:i]
        subtree2 = tree[(i + 1):(close - 1)]
        pos1, child_age1, nodes1 = read_tree(subtree1)
        pos2, child_age2, nodes2 = read_tree(subtree2)
        nodes = nodes1 + nodes2 + [[child_age1, pos, pos1, pos2]]
        nodes.sort() # Sorting based on age
        updated_age = child_age1 + age
        return pos, updated_age, nodes
    else:
        i = 2
        while tree[i] != ',': # , represents the seperation of the left and right
            i += 1
        subtree1 = tree[1:i]
        subtree2 = tree[(i + 1):(close - 1)]
        pos1, child_age1, nodes1 = read_tree(subtree1)
        pos2, child_age2, nodes2 = read_tree(subtree2)
        nodes = nodes1 + nodes2 + [[child_age1, pos, pos1, pos2]]
        nodes.sort() # Sorting based on age
        updated_age = child_age1 + age
        return pos, updated_age, nodes

def adjust_nodes(nodes,nseq):
    """
    Node numbering of read_tree's output does not follow a chronological order.
    This function converts node number in a chronological order.

    Parameters
    ----------
    nodes : list list
        List of nodes read from the raw input
    nseq : TYPE
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

def create_localtree(site, nodes, nseq):
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

    Returns
    -------
    localtree
        Created local tree object

    """
    times = [0] * nseq
    final_nodes = []
    for node in nodes:
        times.append(node[0])
        if node[-1] > node [-2]:
            new_node = (node[-2], node[-1], node[-3])
        else:
            new_node = (node[-1], node[-2], node[-3])
        final_nodes.append(new_node)
    times = [round(time, 4) for time in times] # assuring consistency
    return localtree(site, times, final_nodes)

def read_haplotypes_file(filename):
    """
    

    Parameters
    ----------
    filename : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    f = open(filename)
    inp = f.read()
    f.close()
    lines = inp.split('\n')

    nseq = int(lines[0].split(' ')[1])
    seqlen = int(float(lines[0].split(' ')[2]))
    if len(lines) < 5:
        print("The entered %s file does not contain any trees" % filename)
    # Last nseq lines containts genetic info, above 2 lines indicate segsites
    # and position, +1 because the last item in lines list is ''
    trees = lines[4:-(nseq + 3)]
    position_adjust = round(float(lines[-(nseq + 2)].split(' ')[1]) * seqlen) - 1

    current_site = 1
    newick_trees = []
    for tree in trees:
        pos, newick = tree.split(']')
        newick_trees.append((current_site,newick))
        current_site = int(pos[1:]) + current_site

    newick_trees = [ [(i[0] - position_adjust), i[1]] for i in newick_trees]

    i = 0
    while newick_trees[i][0] < 0:
        i+= 1
    newick_trees = newick_trees[(i - 1):]
    newick_trees[0][0] = 1

    localtrees = []
    for item in newick_trees:
        site = item[0]
        tree = item[1]
        tree = adjust_newick(nseq, tree)
        _,_,nodes =read_tree(tree)
        adjust_nodes(nodes, nseq)
        tree = create_localtree(site, nodes, nseq)
        localtrees.append(tree)

    find_recombination(localtrees)

    return ARG(None, localtrees, "MaCS")
#%%
#os.chdir("r10")
for filename in os.listdir():
    if filename.startswith("haplotypes"):
        break
else:
    print("In the working directory, there should be a file whose name starts with haplotypes")
    sys.exit(1)

#%%
true_arg = read_haplotypes_file(filename)

#%%
def get_true_arg():
    """
    

    Returns
    -------
    true_arg : TYPE
        DESCRIPTION.

    """
    return true_arg

if __name__ == "__main__":
    true_arg.print_arg()
