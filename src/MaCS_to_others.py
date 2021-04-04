import os
#%%
create_arbores = True # Whether the data is also converted into arbores
#%%
for file in os.listdir('.'): # Getting the file which starts with haplotypes
    if file.startswith("haplotypes"):
        break
#%%
f = open(file)
inp = f.read()
f.close()
del f
#%%
file = inp.split('\n')[:-1]
del inp
#%%
sites = []
i = 1
line = file[-i]
while not line.startswith("positions"):
    sites.append(line)
    i += 1
    line = file[-i]
sites.reverse() # Because we are starting adding from the end
#%%
positions = file[-i]
del line, i
#%%
positions = [round(float(i) * 10000) for i in positions.split(' ')[1:]]
#%%
if create_arbores:
    start = positions[0] - 1
    positions = [str(pos - start) for pos in positions]
    del start
else:
    positions = [str(pos) for pos in positions]
#%%
res = " ".join(positions) + '\n'
n0 = ['0' for i in range(len(positions))]
res += " ".join(n0) + '\n'
for site in sites:
    res += " ".join(list(site)) + '\n'
#%%
print(res)
#%%
if create_arbores:
    f = open("Arbores_input.txt",'w')
    f.write(res)
    f.close()
#%% Creating it again to convert to SMARTree
del create_arbores, f, file, n0
res = " ".join(positions) + '\n'
# n0 = ['0' for i in range(len(positions))]
# res += " ".join(n0) + '\n'
for site in sites:
    res += " ".join(list(site)) + '\n'
#%%
inp2 = res[:-1].split("\n")
n = len(inp2) - 1
inp3 = [i.split(' ') for i in inp2]
seqlen = 10000
output = ["{%d, %d}" % (seqlen, n)]
#%%
for i in range(len(inp3[0])):
    col = '{' + inp3[0][i] + ", {"
    for j in range(n):
        col += inp3[j + 1][i] + ", "
    col = col[:-2] + "}, {0.5, 0.5}}"
    output.append(col)
#%%
txt = "\n".join(output)
#%%
print(txt)
#%%
f = open("SMARTree_input.txt", "w")
f.write(txt)
f.close()
#%% Converting it to ARGweaver
del col, f, i, inp2, inp3, j, n, output, res, site, txt
#%%
names = ["NAMES"] + ['n' + str(i) for i in range(len(sites))]
names = '\t'.join(names)
#%%
region = '\t'.join(["REGION", "chr", '1', str(seqlen)])
#%%
lines = []
for i in range(len(positions)):
    line = positions[i] + '\t'
    for site in sites:
        if site[i] == '1':
            line += 'C'
        else:
            line += 'G'
    lines.append(line)
#%%
res = names + '\n' + region + '\n' + '\n'.join(lines) + '\n'
#%%
print(res)
#%%
f = open("ARGweaver_input.sites", 'w')
f.write(res)
f.close()
