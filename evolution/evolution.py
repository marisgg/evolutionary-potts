import subprocess
import pandas as pd
import io
import random
from multiprocessing import Pool, Value, cpu_count
import time
import os
import json
import numpy as np

PARAM_DIR = "./params"
GENERATION_SIZE = 4*12
MUTATION_SCALE = 1
N_ELITE = 2

def mutate_cell(cell):
    cell = cell.copy()
    for attr in cell.keys():
        if attr in ["P", "V", "LAMBDA_P", "LAMBDA_V"]:
            # Do not mutate, keep as is
            continue
        if np.random.choice([True]):
            cell[attr] = max(2, np.round(cell[attr] + MUTATION_SCALE*np.random.standard_cauchy(), decimals=2))
    return cell

def fitness(history):
    df = history
    try:
        startpos = np.array((df["x"].iloc[0],df["y"].iloc[0]))
        endpos = np.array((df["x"].iloc[-1],df["y"].iloc[-1]))
    except IndexError:
        print("Sim failed")
        # Simulation failed - no output
        return -10000
    return np.linalg.norm(endpos-startpos)

def fitness_from_tuple(output):
    output = output.splitlines()[-1]
    def to_int_or_float(s):
        try:
            return int(s)
        except ValueError:
            return float(s)
    spl = output.rstrip("\n").split(",")
    int_tuple = tuple(map(to_int_or_float, spl))
    # fitness = time alive + livelihood left + (100 - distance to nearest food)
    fitness = int_tuple[0]
    return int_tuple[0]+ int_tuple[1] + (100+int_tuple[2])

def run_js_simulation(args):
    (modelname, paramname) = args
    cmd = f"node ./{modelname} {paramname}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    output = proc.communicate()[0]
    output = output.decode("utf-8")
    return output
    # outputIO = io.StringIO(output)
    # df = pd.read_csv(outputIO, sep="\t", header=None, names=["step","id","type","x","y"])
    # return df

def create_param_files(generation):
    j = json.loads('''
                    {
                        "conf":{
                            "MAX_ACT": [0, 0, NaN, 30],
                            "V": [0, 30, NaN, 500],
                            "P": [0, 5, NaN, 260],
                            "LAMBDA_ACT": [0, 0, NaN, 300],
                            "LAMBDA_V": [0, 1000, NaN, 5],
                            "LAMBDA_P": [0, 1, NaN, 2],
                            "LAMBDA_CH": [0, 0, NaN, 500]
                        }
                    }''')
    paramnames = []
    for i, cell in enumerate(generation):
        j["conf"]["MAX_ACT"][-1] = cell["MAX_ACT"]
        j["conf"]["V"][-1] = cell["V"]
        j["conf"]["P"][-1] = cell["P"]
        j["conf"]["LAMBDA_ACT"][-1] = cell["LAMBDA_ACT"]
        j["conf"]["LAMBDA_V"][-1] = cell["LAMBDA_V"]
        j["conf"]["LAMBDA_P"][-1] = cell["LAMBDA_P"]
        j["conf"]["LAMBDA_CH"][-1] = cell["LAMBDA_CH"]
        filename = f"{PARAM_DIR}/{i}.json"
        paramnames.append(filename)
        with open(filename, "w+") as f:
            f.write(json.dumps(j))
    return paramnames

def simulate_generation(generation, modelname, num_procs=12):
    paramnames = create_param_files(generation)
    args = list(map(lambda x: (modelname, x), paramnames))
    with Pool(num_procs) as p:
        sim_results = p.map(run_js_simulation, args)
    fitnesses = map(fitness_from_tuple, sim_results)
    gen_fitnesses = list(zip(generation, fitnesses))
    return gen_fitnesses

def init():
    os.makedirs(f"{PARAM_DIR}", exist_ok=True)

def next_gen_elites_only(generation_with_fitnesses):
    # Sort by increasing fitness
    gen_w_f = sorted(generation_with_fitnesses, key=lambda x: x[1], reverse=True)
    
    gen_w_f = list(map(lambda x: x[0], gen_w_f))
    gen_w_f = gen_w_f[:N_ELITE]

    i = 0
    while len(gen_w_f) < GENERATION_SIZE:
        gen_w_f.append(mutate_cell(gen_w_f[i%N_ELITE]))
        i += 1
    return gen_w_f

def next_generation_elitism_and_roulette(generation_with_fitnesses):
    # Sort by increasing fitness
    gen_w_f = sorted(generation_with_fitnesses, key=lambda x: x[1], reverse=True)
    print(gen_w_f[0])
    fitnesses = list(map(lambda x: x[1], gen_w_f))

    gen_w_f = list(map(lambda x: x[0], gen_w_f))
    elites = gen_w_f[:N_ELITE]

    i = 0
    # while len(gen_w_f) < GENERATION_SIZE:
    #     gen_w_f.append(mutate_cell(gen_w_f[i%N_ELITE]))
    #     i += 1
    sample_weights = np.array(fitnesses)
    sample_weights = sample_weights - np.min(sample_weights) + 100
    sample_weights = sample_weights/sum(sample_weights)
    print(sample_weights)
    sampled_cells = np.random.choice(gen_w_f, size=GENERATION_SIZE-N_ELITE, p=sample_weights)
    next_gen = elites + [mutate_cell(c) for c in sampled_cells]
    return next_gen

def next_generation_elitism_and_inverse_position_sample(generation_with_fitnesses):
    # Sort by increasing fitness
    gen_w_f = sorted(generation_with_fitnesses, key=lambda x: x[1], reverse=True)
    print(gen_w_f[0])
    gen_w_f = list(map(lambda x: x[0], gen_w_f))
    elites = gen_w_f[:N_ELITE]

    i = 0
    # while len(gen_w_f) < GENERATION_SIZE:
    #     gen_w_f.append(mutate_cell(gen_w_f[i%N_ELITE]))
    #     i += 1
    sample_weights = np.array([(1/(x+3))**2 for x in range(GENERATION_SIZE)])
    sample_weights = sample_weights/sum(sample_weights)
    sampled_cells = np.random.choice(gen_w_f, size=GENERATION_SIZE-N_ELITE, p=sample_weights)
    next_gen = elites + [mutate_cell(c) for c in sampled_cells]
    return next_gen

def evolve(filename, num_generations):
    init()
    generation = [{'MAX_ACT': 10, 'P': 250, 'V': 500, 'LAMBDA_ACT': 300, 'LAMBDA_P': 2, 'LAMBDA_V': 5, 'LAMBDA_CH': 20} for i in range(GENERATION_SIZE)]
    generation = list(map(mutate_cell, generation))
    for i in range(num_generations):
        gen_fitnesses = simulate_generation(generation, filename, num_procs=cpu_count())
        generation = next_generation_elitism_and_roulette(gen_fitnesses)

    