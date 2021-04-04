import sys

def find_matching_par(par, string, par_pos):
    """Finds the corresponding '}' to a given '{'. Only works for this
    type of paranthesis. Does not work with "()" or "[]".

    Parameters
    ----------
    string : par
        Paranthesis to be matched
    string : str
        A string which contains paranthesis
    par_pos : int
        Index value of '{' in the string

    Returns
    -------
    i : int
        The index value of the corresponding '}' + 1
    string[par_pos:i] : str
        substring of the input string starting from par_pos to i
    """
    if par == '{':
        matching_par = '}'
    elif par == '(':
        matching_par = ')'
    elif par == '[':
        matching_par = ']'
    else:
        raise Exception('Input Error - The first argument should be a type of paranthesis')

    if string[par_pos] != par:
        raise Exception('Input Error - string[par_pos] should be %s' % par)

    stack = 1 # Number of open parenthesis unmatched
    for i in range(par_pos + 1, len(string)):
        if stack == 0:
            break
        if string[i] == matching_par:
            stack -= 1
        elif string[i] == par:
            stack += 1
    return (i, string[par_pos:i])

def clean_spaces(inp):
    """Cleans the chars (' ', '\n', '\t') in the strings

    Parameters
    ----------
    inp : str
        A string that will be edited

    Returns
    -------
    str
        Cleaned version without (' ', '\n', '\t')
    """
    clean = [inp[0]] # an accumulator for the final clean version
    append = clean.append
    for j in range(1, len(inp)):
        if inp[j] != ' ' and inp[j] != '\n' and inp[j] != '\t':
            append(inp[j])
        elif inp[j] == ' ' and (inp[j - 1] != ' ' and inp[j - 1] != '\n'):
            append(inp[j])
    return ''.join(clean)

class localtree:
    '''A class used to represent local trees

    Attributes
    ----------
    site : int
        site of the tree
    times : float list
        list of coalescence times, index value represent the coalescence points
        (e. g. index 9 represents event number 9)
    TMRCA : float
        time of the most recent common ancestor (TMRCA)
    nodes : integer tuple list
        list of nodes in format of (a, b, c) where a and b coalesces at c
    recomb_exists: bool
        whether this site has recombination or not
    nseq: int
        number of individuals (n) at age 0
    spot: int
        number of individuals (n) + number of coalescence spots
    '''
    def __init__(self, site, times, nodes, anc_map = True):
        """
        

        Parameters
        ----------
        site : TYPE
            DESCRIPTION.
        times : TYPE
            DESCRIPTION.
        nodes : TYPE
            DESCRIPTION.
        anc_map : TYPE, optional
            DESCRIPTION. The default is True.

        Returns
        -------
        None.

        """
        self.site = site
        self.times = times
        self.TMRCA = self.times[-1]
        self.nodes = nodes
        self.recomb_exists = False
        self.nseq = (len(times) + 1) // 2
        self.spot = len(times)
        length = 0
        for node in nodes:
            length += times[node[2]] - times[node[1]]
            length += times[node[2]] - times[node[0]]
        self.branch_length = length
        if anc_map:
            self.ancestor_map = self.create_ancestor_map(self.nseq, nodes)


    def add_recomb(self, recomb_from_individual,recomb_to_individual,
                   recomb_from_time, recomb_to_time):
        # Adds recombination event to the tree
        self.recomb_exists = True
        self.recomb_from_individual = recomb_from_individual
        self.recomb_to_individual = recomb_to_individual
        self.recomb_from_time = recomb_from_time
        self.recomb_to_time = recomb_to_time

    def print_info(self): #Prints information about the tree
        print("Site: " + str(self.site))
        print("Times: " + str(self.times))
        print("Nodes: " + str(self.nodes))
        if self.recomb_exists:
            print("Recombination from %s to %s at time %s" % (self.recomb_from_individual,self.recomb_to_individual,self.recomb_to_time))
        else:
            print("No recombination")

    def get_latex(self): #Creates the latex code for the tree
        # Empty template of Latex which will be filled
        latex = '\\documentclass{standalone}\n\n\\usepackage{tikz}\n\\usetikzlibrary{shapes,arrows,positioning,shapes.geometric}\n\\usepackage{ifthen}\n\n\\setlength{\\textwidth}{40cm}\n\n\\begin{document}\n\\begin{tikzpicture}\n\n\\def\\hstep{.5pt}\n\\def\\vstep{.25pt}\n\n\\def\\pone{%s}\n\\def\\ptwo{%s}\n\\def\\site{%d}\n\\def\\X{%s}\n\\def\\minx{%f}\n\\def\\maxx{%f}\n\\def\\T{%s}\n\\def\\Tnorm{%s}\n\\def\\LAB{%s}\n\\def\\coloring{-1}\n\\def\\globalindex{%s}\n\\def\\globalindiceson{-1}\n\\def\\regrafttime{%f}\n\\def\\regrafttimenorm{%f}\n\\def\\recombinationtime{%f}\n\\def\\recombinationtimenorm{%f}\n\\def\\NL{%d}\n\\def\\NN{%d}\n\\def\\prune{%d}\n\\def\\regraft{%d}\n\n\\def\\prpar{-1}\n\\def\\rgpar{-1}\n\n\\def\\newnode{-1}\n\n\\input{./phytree.tex}\n\n\\end{tikzpicture}\n\\end{document}\n'
        # The output latex file should be compiled with 'phytree.tex' which
        # will be automatically created by print_tree() method

        def adjust_nodes(nodes,n,current_node): #Adjust the placement of nodes for better visualization
        # Idea is from TMCRA to bottom
            if current_node < n:
                return []
            for node in nodes:
                if node[2] == current_node:
                    return adjust_nodes(nodes, n, node[0]) + adjust_nodes(nodes, n, node[1]) + [node]
        d = {i:-1 for i in range(self.spot)} # Stores the x-axis position of spots
        currentmax = 0
        distance_between_bottom_two = 2
        nodes = self.nodes
        adjusted_nodes = adjust_nodes(nodes, self.nseq,nodes[-1][-1])
        for node in adjusted_nodes:
            if d[node[0]] == -1:
                d[node[0]] = currentmax
                currentmax += distance_between_bottom_two
            if d[node[1]] == -1:
                d[node[1]] = currentmax
                currentmax += distance_between_bottom_two
            d[node[2]] = (d[node[0]] + d[node[1]]) / 2
        xlist = list(d.values())
        pone =  '{' + str([tup[0] for tup in nodes])[1:-1] + '}'
        ptwo = '{' + str([tup[1] for tup in nodes])[1:-1] + '}'
        x = '{' + str(xlist)[1:-1] + '}'
        T = '{' + str(self.times)[1:-1] + '}'
        from math import sqrt
        scaled_times = [sqrt(i) * 10 for i in self.times] # For better visualization
        # if scaled_times[-1] < 5:
        #     scaled_times = [i*4 for i in scaled_times]
        Tnorm = '{' + str(scaled_times)[1:-1] + '}'
        LAB = '{' + str(list(range(self.spot)))[1:-1] + '}'
        globalindex = '{' + str(list(range(self.spot-self.nseq)))[1:-1] + '}'
        rec_to_time = -1
        rec_to_time_norm = -1
        rec_from_time = -1
        rec_from_time_norm = -1
        rec_from = -1
        rec_to = -1
        if self.recomb_exists:
            rec_to_time = self.recomb_to_time
            rec_to_time_norm = sqrt(rec_to_time) * 10 # Same scaling as above
            rec_from_time = self.recomb_from_time
            rec_from_time_norm = sqrt(rec_from_time) * 10 # Same scaling as above
            rec_from = self.recomb_from_individual
            rec_to = self.recomb_to_individual
        out = latex % (pone, ptwo, self.site, x, min(xlist), max(xlist), T, Tnorm,LAB, globalindex, \
                       rec_to_time, rec_to_time_norm, rec_from_time, rec_from_time_norm, \
                       self.nseq, self.spot, rec_from, rec_to)
        return out

    def print_tree(self, keep_latex = False, keep_log = False, keep_aux = False):
        # Creates pdf visualization for the tree by using Latex
        import os

        # path to be created
        path = "./trees"

        if not os.path.exists(path):
            try:
                os.mkdir(path)
            except OSError:
                print ("Creation of the directory %s failed" % path)
            else:
                print ("Successfully created the directory %s " % path)

        os.chdir(path)

        # Required tex file for visualizing the tree
        f = open("phytree.tex", 'w')
        phytree = '\\pgfmathsetmacro{\\xsite}{(\\minx+\\maxx)/2*\\hstep}\n\\pgfmathsetmacro{\\ysite}{-2*\\vstep}\n\\pgfmathsetmacro{\\sitelab}{\\site}\n\\node[scale=.5] at (\\xsite,\\ysite) {\\pgfmathprintnumber[fixed,precision=0]{\\sitelab}};\n\\pgfmathsetmacro{\\nb}{\\NL-2}\n\\foreach \\i in {0,...,\\nb} \n{\n\t\\pgfmathsetmacro{\\y}{\\Tnorm[\\NL+\\i]*\\vstep}\n\t\\pgfmathsetmacro{\\xs}{\\minx*\\hstep}\n\t\\pgfmathsetmacro{\\xe}{\\maxx*\\hstep}\n\t\\pgfmathsetmacro{\\xt}{(\\minx-1)*\\hstep}\n\t\\pgfmathsetmacro{\\tlab}{\\T[\\NL+\\i]}\n\t\n\t\\draw[-,dashed] (\\xs,\\y) -- (\\xe,\\y);\n\t\n\t\\node[scale=.5] at (\\xt,\\y) {\\pgfmathprintnumber[fixed,precision=3]{\\tlab}};\n}\n\\pgfmathsetmacro{\\maxit}{\\NN-1}\n\\foreach \\i in {0,...,\\maxit}\n{\n\t\\pgfmathsetmacro{\\x}{\\X[\\i]*\\hstep}\n\t\\pgfmathsetmacro{\\y}{\\Tnorm[\\i]*\\vstep}\n\t\n\t\\ifthenelse{\\i=\\NL \\OR \\i>\\NL}{\n\t\n\t\t\\pgfmathsetmacro{\\signL}{int(\\pone[\\i-\\NL])}\n\t\t\\pgfmathsetmacro{\\signR}{int(\\ptwo[\\i-\\NL])}\n\t\t\t\n\t\t\\ifthenelse{\\signL>-1 \\AND \\signR>-1}{\n\t\t\t\\pgfmathsetmacro{\\xl}{\\X[\\pone[\\i-\\NL]]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yl}{\\Tnorm[\\pone[\\i-\\NL]]*\\vstep}\n\t\t\t\\pgfmathsetmacro{\\xr}{\\X[\\ptwo[\\i-\\NL]]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yr}{\\Tnorm[\\ptwo[\\i-\\NL]]*\\vstep}\n\t\t\t\\draw[-,line width = 1pt] (\\xl,\\yl) -- (\\xl,\\y) -- (\\xr,\\y) -- (\\xr,\\yr);\t\t\n\t\t}{};\n\t\t\n\t\t\\ifthenelse{\\signL<0 \\AND \\signR>-1}{\n\t\t\t\\pgfmathsetmacro{\\xr}{\\X[\\ptwo[\\i-\\NL]]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yr}{\\Tnorm[\\ptwo[\\i-\\NL]]*\\vstep}\n\t\t\t\\draw[-,line width = 1pt] (\\xr,\\y) -- (\\xr,\\yr);\t\t\n\t\t}{};\n\t\t\\ifthenelse{\\signL>-1 \\AND \\signR<0}{\n\t\t\t\\pgfmathsetmacro{\\xl}{\\X[\\pone[\\i-\\NL]]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yl}{\\Tnorm[\\pone[\\i-\\NL]]*\\vstep}\n\t\t\t\\draw[-,line width = 1pt] (\\xl,\\yl) -- (\\xl,\\y);\t\t\n\t\t}{};\n\t}{};\t\t\n}\n\\ifthenelse{\\coloring>-1}{\n\\foreach \\i in {0,...,\\maxit}\n{\n\t\\pgfmathsetmacro{\\x}{\\X[\\i]*\\hstep}\n\t\\pgfmathsetmacro{\\y}{\\Tnorm[\\i]*\\vstep}\n\t\\pgfmathsetmacro{\\clr}{\\colors[\\i]}\n\t\\ifthenelse{\\clr=2}{\n\t\\node[draw,circle,fill=white, inner sep =0pt, minimum size = 2mm] at (\\x,\\y) {};\n\t}{\n\t\\ifthenelse{\\clr=0}{\n\t\\node[draw,circle,fill=black, inner sep =0pt, minimum size = 2mm] at (\\x,\\y) {};\n\t}{\n\t\\node[draw,circle,fill=gray, inner sep =0pt, minimum size = 2mm] at (\\x,\\y) {};\n\t};\n\t};\t\n}\n}{};\n\\foreach \\i in {0,...,\\maxit}\n{\t\n\t\\pgfmathsetmacro{\\x}{\\X[\\i]*\\hstep}\n\t\\pgfmathsetmacro{\\y}{\\Tnorm[\\i]*\\vstep-1*\\vstep}\n\t\\pgfmathsetmacro{\\lab}{\\LAB[\\i]}\t\n\n\t\\node[fill=white,inner sep =0pt, scale=.5] at (\\x,\\y) {\\pgfmathprintnumber[fixed,precision=0]{\\lab}};\n\t\n}\n\\ifthenelse{\\globalindiceson>0}{\n\\pgfmathsetmacro{\\nbmaxit}{\\NL-2}\n\\foreach \\i in {0,...,\\nbmaxit}\n{\t\n\t\\pgfmathsetmacro{\\x}{\\X[\\i+\\NL]*\\hstep+.5*\\hstep}\n\t\\pgfmathsetmacro{\\y}{\\Tnorm[\\i+\\NL]*\\vstep+.5*\\vstep}\n\t\\pgfmathsetmacro{\\glab}{\\globalindex[\\i]}\n\t\\node[fill=white,inner sep =0pt, scale=.5] at (\\x,\\y) {\\pgfmathprintnumber[fixed,precision=0]{\\glab}};\t\n}\n}{};\n\\ifthenelse{\\prune>-1}{\n\n\t\\pgfmathsetmacro\\rectimesign{int(round(\\recombinationtime))}\n\t\\ifthenelse{\\rectimesign<0}{\n\t\t\\pgfmathsetmacro{\\xstart}{\\X[\\prune]*\\hstep}\n\t\t\\pgfmathsetmacro{\\ystart}{\\Tnorm[\\prune]*\\vstep + (\\Tnorm[\\prpar] - \\Tnorm[\\prune])/2*\\vstep}\n\t}{\n\t\t\\pgfmathsetmacro{\\xstart}{\\X[\\prune]*\\hstep}\n\t\t\\pgfmathsetmacro{\\ystart}{\\recombinationtimenorm*\\vstep}\n\t};\n\t\\node[scale=0.5] at (\\xstart,\\ystart) {$\\times$};\t\n\t\\pgfmathsetmacro\\rgtimesign{int(round(\\regrafttime))}\n\t\\ifthenelse{\\rgtimesign<0}{\n\t\n\t\t\\ifthenelse{\\rgpar>-1}{\n\t\t\t\\pgfmathsetmacro{\\xend}{\\X[\\regraft]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yend}{\\Tnorm[\\regraft]*\\vstep + (\\Tnorm[\\rgpar] - \\Tnorm[\\regraft])/2*\\vstep}\n\t\t}{\n\t\t\t\\pgfmathsetmacro{\\xend}{\\X[\\regraft]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yend}{\\Tnorm[\\NN-1]*\\vstep}\n\t\t};\n\t}{\n\t\t\\pgfmathsetmacro{\\xend}{\\X[\\regraft]*\\hstep}\n\t\t\\pgfmathsetmacro{\\yend}{\\regrafttimenorm*\\vstep}\t\t\n\t};\n\t\\pgfmathsetmacro{\\xmid}{((\\xstart+\\xend)/2)}\n\t\\draw[color=red,->] (\\xstart,\\ystart) -- (\\xmid,\\ystart) -- (\\xmid,\\yend) -- (\\xend,\\yend);\t\n}{};'
        f.write(phytree)
        f.close()

        latex_file = "Tree_Site%d_Rec%s" % (self.site,self.recomb_exists)
        f = open(latex_file + ".tex", 'w')
        f.write(self.get_latex())
        f.close()

        try:
            os.system("pdflatex " + latex_file + ".tex")
        except:
            print("Problem with pdflatex command: Please check that you have the latest version of Latex")
        else:
            if not keep_aux:
                os.remove(latex_file + ".aux")
            if not keep_log:
                os.remove(latex_file + ".log")
            if not keep_latex:
                os.remove(latex_file + ".tex")

        os.remove("phytree.tex")
        os.chdir("..")
        print("Tree for site", self.site, "has been printed")

    def get_newick(self):
        def helper(mrca):
            index = mrca - self.nseq
            if index < 0:
                return ''
            l, r, _ = self.nodes[index]
            l_time = round(self.times[mrca] - self.times[l], 5)
            r_time = round(self.times[mrca] - self.times[r], 5)
            l_tree = helper(l)
            r_tree = helper(r)
            tree = "(%s%d:%f,%s%d:%f)" % (l_tree, l, l_time, r_tree, r, r_time)
            return tree
        return helper(len(self.times) - 1) + ';'

    def robinson_foulds(self, tr):
        from ete3 import Tree
        t1 = self.get_newick()
        t2 = tr.get_newick()
        t1 = Tree(t1)
        t2 = Tree(t2)
        return t1.robinson_foulds(t2)[0] # First element gives the difference

    def is_same(self, tr):
        """
        

        Parameters
        ----------
        tr : TYPE
            DESCRIPTION.

        Returns
        -------
        bool
            DESCRIPTION.

        """
        if self.nodes == tr.nodes and \
        self.times == tr.times:
            return True
        return False

    def create_ancestor_map(self, nseq, nodes):
        """ Creates a nseq * nseq matrix where index [i][j] represent the number
        of common ancestors individual i and j have. This matrix will be used to
        find recombination in extreme cases.

        Parameters
        ----------
        nseq : int
            number of individuals (n) at age 0
        nodes : (int * int * int) list
            list of nodes in format of (a, b, c) where a and b coalesces at c

        Returns
        -------
        common_ancestor_map : (int list) list
            nseq * nseq matrix where index [i][j] represent the number of
            common ancestors individual i and j have

        """
        tmrca_map = [] # TMCRA map for all n, index number represents n number
        for n in range(nseq):
            _map = [n]
            current = n
            while current != (nseq * 2 - 2): # TMCRA is (nseq * 2 - 2)
                for item in nodes:
                    if current in item:
                        _map.append(item[-1])
                        current = item[-1]
            tmrca_map.append(_map)

        # Stores the number of common ancestors between two individual
        common_ancestor_map = [[0] * nseq for i in range(nseq)]

        for n1 in range(nseq):
            for n2 in range(n1, nseq):
                common_ancestor_map[n1][n2] = len(set(tmrca_map[n1]).intersection(tmrca_map[n2]))
                common_ancestor_map[n2][n1] = common_ancestor_map[n1][n2]

        return common_ancestor_map
#%%

def find_recombination(trees):
    """Given list of localtree objects, this function finds recombinations
    and edits the localtree objects

    Parameters
    ----------
    trees : localtree list
        list of localtree objects

    Returns
    -------
    None
    """
    for i in range(len(trees) - 1):
        try:
            #if trees[i].robinson_foulds(trees[i + 1]) > 0:
            if trees[i].nodes != trees[i + 1].nodes: # If topology of the adjacent tree is not the same
                diff_set = set(trees[i].times) - set(trees[i + 1].times)
                if diff_set == set(): # In case the lost time is not unique
                    #print("Not unique time set")
                    diff_set = (set(trees[i + 1].times) - set(trees[i].times))
                    if diff_set == set(): # In case the lost time is not unique
                        #print("Times are completely the same")
                        map1 = trees[i].ancestor_map
                        map2 = trees[i + 1].ancestor_map
                        diff = 0 # Number of differences
                        diff_pos = 0 # The position of differences in case the difference is unique
                        for j in range(len(map1[0])): # Checking for number 0
                            if map1[0][j] != map2[0][j]:
                                diff += 1
                                diff_pos = j
                        if diff > 1:
                            start = 0
                        else: # diff == 1
                            start = diff_pos
                        time1 = trees[i].times[start]
                        for node in trees[i + 1].nodes:
                            if start in node and start != node[2]:
                                end_time = trees[i].times[node[2]]
                                if node[0] == start:
                                    end = node[1]
                                else:
                                    end = node[0]
                                break
                        time2 = trees[i + 1].times[end] # Time on tree i + 1
                        if time2 != 0:
                            end = trees[i].times.index(time2)
                        start_time = (time1 + end_time) / 2
                        trees[i].add_recomb(start, end, start_time, end_time)
                        continue

                    end_time = diff_set.pop()
                    for item in trees[i + 1].nodes: ####### Fix thisss = smc5
                        if trees[i + 1].times[item[2]] == end_time:
                            cand1 = item[0] # Canditate 1
                            cand2 = item[1] # Canditate 2
                            time1 = trees[i + 1].times[cand1] # Time on tree i + 1
                            time2 = trees[i + 1].times[cand2] # Time on tree i + 1
                            start = trees[i].times.index(time1) # Label on tree i
                            end = trees[i].times.index(time2) # Label on tree i

                    start_time = (time1 + end_time) / 2
                    trees[i].add_recomb(start, end, start_time, end_time)
                    continue

                lost_time = diff_set.pop()  # lost time in tree
                lost_index = trees[i].times.index(lost_time) # lost coalescence point

                # print("lost_time", lost_time, "lost_index", lost_index)
                for item in trees[i].nodes: # Find start node of recombination
                    if lost_index == item[2]:
                        temp = (item[0], item[1]) # Two points that coalesces at lost_index
                        if lost_time ==  trees[i].times[-1]: # If TMCRA is lost
                            start = item[0]
                            break
                    elif lost_index in item:
                        pretime = trees[i].times[item[2]] # where lost_index goes up in tree
                        pre = trees[i + 1].times.index(pretime) # finding the person at "pretime" in the next tree
                        for item in trees[i + 1].nodes: # Find the start point of recombination
                        ### Check whether this is true, mostly likely it is, kendinden kucukle recomb olmuyor, mantiksiz
                            if pre == item[2]:
                                if temp[0] in item:
                                    start = temp[1]
                                    break
                                else:
                                    start = temp[0]
                                    break
                        break

                start_time = (lost_time + trees[i].times[start]) / 2 # Starting time of the recombination

                # After recombination, coalescence index number might change
                start_in_nxt = start
                if start >= trees[i].nseq: # This is not the case if item < n.
                    start_in_nxt = trees[i + 1].times.index(trees[i].times[start])

                for item in trees[i + 1].nodes: # Find the end point of recombination
                    if start_in_nxt in item and start_in_nxt != item[2]:
                        if item[0] == start_in_nxt:
                            end = item[1]
                        else:
                            end = item[0]
                        # After recombination, coalescence index number might change
                        # This is not the case if item < n. Otherwise, we have to
                        # find the correct time
                        if end >= trees[i].nseq:
                            end = trees[i].times.index(trees[i + 1].times[end])
                        index = item[2] #To find the recombination time

                end_time = trees[i + 1].times[index] # End time of the recombination

                trees[i].add_recomb(start, end, start_time, end_time)
            elif not trees[i].is_same(trees[i + 1]):
                if trees[i].nodes == trees[i + 1].nodes: # Same nodes, different time
                    diff_set = (set(trees[i].times) - set(trees[i + 1].times)) # lost time in tree
                    if diff_set == set(): # In case the lost time is not unique
                        #print("Not unique time set")
                        diff_set = (set(trees[i + 1].times) - set(trees[i].times))
                        if diff_set == set(): # In case the lost time is not unique
                            print("This is a very strange case!")
                            sys.exit(1)
                            #trees[i].add_recomb(start, end, start_time, end_time)
                            continue

                        end_time = diff_set.pop()
                        lost_index = trees[i + 1].times.index(end_time)
                        lost_time = trees[i].times[lost_index]
                        for node in trees[i + 1].nodes:
                            if node[2] == lost_index:
                                if trees[i].times[node[0]] > trees[i].times[node[1]]:
                                    # If RCA goes up
                                    if end_time > lost_time:
                                        end = lost_index
                                        start = node[1]
                                        start_time = (lost_time + trees[i].times[start]) / 2
                                    else: # If RCA goes down ######
                                        end = node[0]
                                        start = node[1]
                                        start_time = (end_time + trees[i].times[start]) / 2
                                else:
                                    # If RCA goes up
                                    if end_time > lost_time:
                                        end = lost_index
                                        start = node[0]
                                        start_time = (lost_time + trees[i].times[start]) / 2
                                    else: # If RCA goes down ######
                                        end = node[1]
                                        start = node[0]
                                        start_time = (end_time + trees[i].times[start]) / 2

                        trees[i].add_recomb(start, end, start_time, end_time)
                        continue

                    lost_time = diff_set.pop()
                    lost_index = trees[i].times.index(lost_time) # lost index
                    end_time = trees[i + 1].times[lost_index]

                    for node in trees[i].nodes:
                        if lost_index == node[2]:
                            if trees[i].times[node[0]] > trees[i].times[node[1]]:
                                # If RCA goes up
                                if trees[i + 1].times[lost_index] > lost_time:
                                    end = lost_index
                                    start = node[1]
                                    ## Because it is going up, adjusting the start point,
                                    # so that it does not go above the line
                                    start_time = (lost_time + trees[i].times[start]) / 2
                                else: # If RCA goes down ######
                                    end = node[0]
                                    start = node[1]
                                    start_time = (end_time + trees[i].times[start]) / 2
                            else:
                                # If RCA goes up
                                if trees[i + 1].times[lost_index] > lost_time:
                                    end = lost_index
                                    start = node[0]
                                    start_time = (lost_time + trees[i].times[start]) / 2 ## Mesela evet, her zaman bu yapilabilir
                                else: # If RCA goes down ######
                                    end = node[1]
                                    start = node[0]
                                    start_time = (end_time + trees[i].times[start]) / 2
                            break

                    trees[i].add_recomb(start, end, start_time, end_time)

        except:
            trees[i].print_info()
            trees[i + 1].print_info()
            print("Problem is at tree", i)
            #sys.exit(1)
            continue


#%%
def get_trees_with_recomb(trees): #Returns a list localtree object that has unique topology
    lst = []
    for t in trees:
        if t.recomb_exists:
            lst.append(t)
    return lst

def get_trees_with_recomb_and_after_recomb(trees): #Returns a list localtree object that has recombination
    lst = []
    for i in range(len(trees)):
        if trees[i].recomb_exists:
            lst.append(trees[i])
            from copy import deepcopy
            nxt = deepcopy(trees[i + 1])
            nxt.recomb_exists = False
            lst.append(nxt)
    return lst

class ARG:
    '''A class used to represent ARG

    Attributes
    ----------
    posterior : float
        log posterior probabilty of the ARG
    localtrees : localtree list
        list of localtrees in the ARG
    nrecombs: int
        number of recombinations
    nseq: int
        number of individuals (n) at age 0
    '''
    def __init__(self, posterior, localtrees, algorithm):
        self.posterior = posterior
        self.localtrees = localtrees
        self.nrecombs = len(localtrees) - 1
        self.nseq = localtrees[0].nseq
        self.algorithm = algorithm


    def print_info(self):
        print(self.get_info())

    def get_info(self):
        info = "Log posterior: " + str(self.posterior) + '\n'
        info += "Number of localtrees: " + str(len(self.localtrees)) + '\n'
        info += "Number of recombs: " + str(self.nrecombs) + '\n'
        return info

    def get_latex(self, nrow = 0, max_row = 4):
        """Creates latex codes for local trees of an ARG. It creates latex
        by grouping local trees based on the input nrow and max_row.
        The result will look like a matrix of local trees.

        Parameters
        ----------
        nrow : int, optional
            Number of rows in the printed ARG. The default is 0.
        max_row : int, optional
            Maximum number of local trees allowed in one row. The default is 4.

        Returns
        -------
        latex : str list
            List of latex codes for parts of the ARG

        """
        import numpy as np
        trees = self.localtrees

        if nrow == 0: # Automatically adjust number of rows
            nrow = (len(trees) - 1) // max_row + 1

        n_trees = (len(trees) + (nrow - 1)) // nrow # Number of trees in one row

        latex = []

        for row in range(nrow):

            if row == nrow - 1:
                tr = trees[(row * n_trees):]
            else:
                tr = trees[(row * n_trees):((row + 1) * n_trees)]

            start = tr[0].get_latex()

            s = start.split("vstep{.25pt}")
            beginning = s[0] + "vstep{.25pt}"

            e = s[1].split("{./phytree.tex}")
            tree = ''
            end = e[1]

            nseq = tr[0].nseq
            compress_factor = 2 # Compress it so that nodes are closer
            shift_factor = nseq * 2 # vertical shift for the tree
            last_shift = (compress_factor**2 - 2) / compress_factor # to make sure distance between trees is 2
            for index in range(0, len(tr)):
                start = tr[index].get_latex()

                s = start.split("vstep{.25pt}")
                beginning = s[0] + "vstep{.25pt}"

                e = s[1].split("{./phytree.tex}")
                t1 = e[0] + "{./phytree.tex}"
                end = e[1]
                t = t1.split("\X{")
                t_start = t[0] + "\X{"

                i,xlist = find_matching_par('{', t[1], 0)

                xlist = np.array([xlist[1:-1].split(", ")], dtype=float)
                xlist = (xlist + shift_factor * index) / compress_factor # adjust positions
                xlist = xlist + last_shift * index

                xlist = [str(num) for num in xlist.tolist()[0]]

                xlist = '{' + ", ".join(xlist) + '}'

                temp = t[1][i:].split("minx")
                mid1 = temp[0] + "minx"

                j,minx = find_matching_par('{', temp[1], 0)

                minx = (float(minx[1:-1]) + (shift_factor * index)) / compress_factor
                minx = minx + last_shift * index
                minx = '{' + str(minx)  + '}'

                temp = temp[1][j:].split("maxx")
                mid2 = temp[0] + "maxx"

                k,maxx = find_matching_par('{', temp[1], 0)

                maxx = (float(maxx[1:-1]) + (shift_factor * index)) / compress_factor
                maxx = maxx + last_shift * index
                maxx = '{' + str(maxx)  + '}'

                temp[1][k:]

                new_tree = t_start + xlist + mid1 + minx + mid2 + maxx + temp[1][k:]

                tree += new_tree + "\n"

            final = beginning + tree + end
            latex.append(final)

        return latex

    def print_arg(self, dir_name = None, keep_latex = False, keep_log = False, keep_aux = False):
            # Creates pdf visualization for the ARG by using Latex
            import os

            # path to be created
            if not dir_name:
                dir_name = "./trees"

            if not os.path.exists(dir_name):
                try:
                    os.mkdir(dir_name)
                except OSError:
                    print ("Creation of the directory %s failed" % dir_name)
                else:
                    print ("Successfully created the directory %s " % dir_name)

            os.chdir(dir_name)

            # Required tex file for visualizing the tree
            f = open("phytree.tex", 'w')
            phytree = '\\pgfmathsetmacro{\\xsite}{(\\minx+\\maxx)/2*\\hstep}\n\\pgfmathsetmacro{\\ysite}{-2*\\vstep}\n\\pgfmathsetmacro{\\sitelab}{\\site}\n\\node[scale=.5] at (\\xsite,\\ysite) {\\pgfmathprintnumber[fixed,precision=0]{\\sitelab}};\n\\pgfmathsetmacro{\\nb}{\\NL-2}\n\\foreach \\i in {0,...,\\nb} \n{\n\t\\pgfmathsetmacro{\\y}{\\Tnorm[\\NL+\\i]*\\vstep}\n\t\\pgfmathsetmacro{\\xs}{\\minx*\\hstep}\n\t\\pgfmathsetmacro{\\xe}{\\maxx*\\hstep}\n\t\\pgfmathsetmacro{\\xt}{(\\minx-1)*\\hstep}\n\t\\pgfmathsetmacro{\\tlab}{\\T[\\NL+\\i]}\n\t\n\t\\draw[-,dashed] (\\xs,\\y) -- (\\xe,\\y);\n\t\n\t\\node[scale=.5] at (\\xt,\\y) {\\pgfmathprintnumber[fixed,precision=3]{\\tlab}};\n}\n\\pgfmathsetmacro{\\maxit}{\\NN-1}\n\\foreach \\i in {0,...,\\maxit}\n{\n\t\\pgfmathsetmacro{\\x}{\\X[\\i]*\\hstep}\n\t\\pgfmathsetmacro{\\y}{\\Tnorm[\\i]*\\vstep}\n\t\n\t\\ifthenelse{\\i=\\NL \\OR \\i>\\NL}{\n\t\n\t\t\\pgfmathsetmacro{\\signL}{int(\\pone[\\i-\\NL])}\n\t\t\\pgfmathsetmacro{\\signR}{int(\\ptwo[\\i-\\NL])}\n\t\t\t\n\t\t\\ifthenelse{\\signL>-1 \\AND \\signR>-1}{\n\t\t\t\\pgfmathsetmacro{\\xl}{\\X[\\pone[\\i-\\NL]]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yl}{\\Tnorm[\\pone[\\i-\\NL]]*\\vstep}\n\t\t\t\\pgfmathsetmacro{\\xr}{\\X[\\ptwo[\\i-\\NL]]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yr}{\\Tnorm[\\ptwo[\\i-\\NL]]*\\vstep}\n\t\t\t\\draw[-,line width = 1pt] (\\xl,\\yl) -- (\\xl,\\y) -- (\\xr,\\y) -- (\\xr,\\yr);\t\t\n\t\t}{};\n\t\t\n\t\t\\ifthenelse{\\signL<0 \\AND \\signR>-1}{\n\t\t\t\\pgfmathsetmacro{\\xr}{\\X[\\ptwo[\\i-\\NL]]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yr}{\\Tnorm[\\ptwo[\\i-\\NL]]*\\vstep}\n\t\t\t\\draw[-,line width = 1pt] (\\xr,\\y) -- (\\xr,\\yr);\t\t\n\t\t}{};\n\t\t\\ifthenelse{\\signL>-1 \\AND \\signR<0}{\n\t\t\t\\pgfmathsetmacro{\\xl}{\\X[\\pone[\\i-\\NL]]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yl}{\\Tnorm[\\pone[\\i-\\NL]]*\\vstep}\n\t\t\t\\draw[-,line width = 1pt] (\\xl,\\yl) -- (\\xl,\\y);\t\t\n\t\t}{};\n\t}{};\t\t\n}\n\\ifthenelse{\\coloring>-1}{\n\\foreach \\i in {0,...,\\maxit}\n{\n\t\\pgfmathsetmacro{\\x}{\\X[\\i]*\\hstep}\n\t\\pgfmathsetmacro{\\y}{\\Tnorm[\\i]*\\vstep}\n\t\\pgfmathsetmacro{\\clr}{\\colors[\\i]}\n\t\\ifthenelse{\\clr=2}{\n\t\\node[draw,circle,fill=white, inner sep =0pt, minimum size = 2mm] at (\\x,\\y) {};\n\t}{\n\t\\ifthenelse{\\clr=0}{\n\t\\node[draw,circle,fill=black, inner sep =0pt, minimum size = 2mm] at (\\x,\\y) {};\n\t}{\n\t\\node[draw,circle,fill=gray, inner sep =0pt, minimum size = 2mm] at (\\x,\\y) {};\n\t};\n\t};\t\n}\n}{};\n\\foreach \\i in {0,...,\\maxit}\n{\t\n\t\\pgfmathsetmacro{\\x}{\\X[\\i]*\\hstep}\n\t\\pgfmathsetmacro{\\y}{\\Tnorm[\\i]*\\vstep-1*\\vstep}\n\t\\pgfmathsetmacro{\\lab}{\\LAB[\\i]}\t\n\n\t\\node[fill=white,inner sep =0pt, scale=.5] at (\\x,\\y) {\\pgfmathprintnumber[fixed,precision=0]{\\lab}};\n\t\n}\n\\ifthenelse{\\globalindiceson>0}{\n\\pgfmathsetmacro{\\nbmaxit}{\\NL-2}\n\\foreach \\i in {0,...,\\nbmaxit}\n{\t\n\t\\pgfmathsetmacro{\\x}{\\X[\\i+\\NL]*\\hstep+.5*\\hstep}\n\t\\pgfmathsetmacro{\\y}{\\Tnorm[\\i+\\NL]*\\vstep+.5*\\vstep}\n\t\\pgfmathsetmacro{\\glab}{\\globalindex[\\i]}\n\t\\node[fill=white,inner sep =0pt, scale=.5] at (\\x,\\y) {\\pgfmathprintnumber[fixed,precision=0]{\\glab}};\t\n}\n}{};\n\\ifthenelse{\\prune>-1}{\n\n\t\\pgfmathsetmacro\\rectimesign{int(round(\\recombinationtime))}\n\t\\ifthenelse{\\rectimesign<0}{\n\t\t\\pgfmathsetmacro{\\xstart}{\\X[\\prune]*\\hstep}\n\t\t\\pgfmathsetmacro{\\ystart}{\\Tnorm[\\prune]*\\vstep + (\\Tnorm[\\prpar] - \\Tnorm[\\prune])/2*\\vstep}\n\t}{\n\t\t\\pgfmathsetmacro{\\xstart}{\\X[\\prune]*\\hstep}\n\t\t\\pgfmathsetmacro{\\ystart}{\\recombinationtimenorm*\\vstep}\n\t};\n\t\\node[scale=0.5] at (\\xstart,\\ystart) {$\\times$};\t\n\t\\pgfmathsetmacro\\rgtimesign{int(round(\\regrafttime))}\n\t\\ifthenelse{\\rgtimesign<0}{\n\t\n\t\t\\ifthenelse{\\rgpar>-1}{\n\t\t\t\\pgfmathsetmacro{\\xend}{\\X[\\regraft]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yend}{\\Tnorm[\\regraft]*\\vstep + (\\Tnorm[\\rgpar] - \\Tnorm[\\regraft])/2*\\vstep}\n\t\t}{\n\t\t\t\\pgfmathsetmacro{\\xend}{\\X[\\regraft]*\\hstep}\n\t\t\t\\pgfmathsetmacro{\\yend}{\\Tnorm[\\NN-1]*\\vstep}\n\t\t};\n\t}{\n\t\t\\pgfmathsetmacro{\\xend}{\\X[\\regraft]*\\hstep}\n\t\t\\pgfmathsetmacro{\\yend}{\\regrafttimenorm*\\vstep}\t\t\n\t};\n\t\\pgfmathsetmacro{\\xmid}{((\\xstart+\\xend)/2)}\n\t\\draw[color=red,->] (\\xstart,\\ystart) -- (\\xmid,\\ystart) -- (\\xmid,\\yend) -- (\\xend,\\yend);\t\n}{};'
            f.write(phytree)
            f.close()

            latex = self.get_latex()
            i = 1

            for file in latex:
                latex_file = "ARG_part" + str(i) # extensions will be added later
                f = open(latex_file + ".tex", 'w')
                f.write(file)
                f.close()

                try:
                    os.system("pdflatex " + latex_file + ".tex")
                except:
                    print("Problem with pdflatex command: Please check that you have the latest version of Latex")
                else:
                    if not keep_aux:
                        os.remove(latex_file + ".aux")
                    if not keep_log:
                        os.remove(latex_file + ".log")
                    if not keep_latex:
                        os.remove(latex_file + ".tex")

                print("Part", i, "is printed")
                i += 1

            os.remove("phytree.tex")
            os.chdir("..")
            print('ARG has successfully been printed!')
