import numpy as np
import os
from matplotlib import pyplot as plt


result_dirs = {"none":[], "medium":[], "high":[]}
for filename in os.listdir("results/"):
    filename = "results/"+filename
    if os.path.isdir(filename):
        print(filename)
        if filename.startswith("results"):
            runtype = filename.split("-")[1]
            result_dirs[runtype].append(filename)


for runtype in result_dirs.keys():
    if(len(result_dirs[runtype]) == 0):
        continue
    fig, axs = plt.subplots(1,len(result_dirs[runtype]), sharey="all")
    labels = []
    for i, dirname in enumerate(result_dirs[runtype]):
        with open(f"{dirname}/evolutionary-potts/out.txt", "r") as f:
            desired_lines = []
            fitnessdata = []
            for line in f:
                if line.startswith("Genfitness"):
                    words = [w for w in line.rstrip().split(" ") if w != '']
                    if words[0] == 'Genfitness':
                        words = words[1:]
                    words = [word.rstrip(':') for word in words]
                    words = [word for word in words if word != 'dev']
                    keys = words[::2]
                    values = list(map(float, words[1::2]))
                    linedict = dict(zip(keys, values))
                    fitnessdata.append(linedict)
            lcl_labels = []
            for key in fitnessdata[0].keys():
                lcl_labels.append(key)
                axs[i].plot(range(len(fitnessdata)), [gen[key] for gen in fitnessdata], label=key)
            labels = lcl_labels
            plt.gca().set_ylim([-500,9000])
    fig.legend(labels)
    fig.suptitle(runtype)
    plt.tight_layout()
    plt.savefig(f"results/{runtype}_evohistory.png")
    plt.show()
    