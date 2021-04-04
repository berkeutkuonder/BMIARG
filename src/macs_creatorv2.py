nseq = 8
mu = ".00048"
rho = str(float(mu) / 1)
seqlen = '1e4'

#%%
import os
it = 0
info = [("Seed","nsites","nrecomb","unique_n")]

lst = [15102,22480,43784,45118,112026,138939,151520,157202,168963,189812]

for it in lst:
# while it < 100000000:
#     it += 1
    command = "./macs %d %s -T -t %s -r %s -h 1e2 -s %d 2>trees.txt | ./msformatter > haplotypes.txt" % (nseq,seqlen,mu,rho,it)
    os.system(command)


    f = open("haplotypes.txt")
    inp = f.read()
    f.close()

    file = inp.split('\n')[:-1]

    if len(file) < 4:
        continue

    nsites = int(file[-(nseq+2)].split(": ")[1])
    sites = file[-nseq:]
    nrecombs = len(file) - nseq - 2 - 4 - 1 # 2 lines above seqs, first 4 lines, -1 for recombs
    narbsites = 0

    # if nrecombs < 4 or nrecombs > 10:
    #     continue

    if nsites <= nrecombs + 3:
        continue

    for i in range(nsites):
        freq = 0 # Frequency of 1s
        for j in range(nseq):
            if sites[j][i] == '1':
                freq += 1
        if freq > 1:
            narbsites += 1

    unique_n = len(set(sites))
    if len(sites) - unique_n > 2:
        continue

    if nsites == narbsites:
        info.append((str(it), str(nsites), str(nrecombs), str(unique_n)))

    print("Iter no: %s Found: %d" % (it,len(info)-1))

    if len(info)-1 == 50:
        print("Enough data found!")
        break

res = "\n".join(["\t".join(i) for i in info])
f = open("stats.txt", "w")
f.write(res)
f.close()
