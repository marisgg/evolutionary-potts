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
MUTATION_SCALE = .5
N_ELITE = 2
POTTS_SEED = 1
CHANGE_POTTS_SEED_PER_GEN = True
GENERATION_NO = 1
np.random.seed(POTTS_SEED)

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
    # fitness = time alive-minimum lifetime + livelihood left + (200 - distance to nearest food)
    # fitness = int_tuple[0]
    return int_tuple[0]+ int_tuple[1] + (200+int_tuple[2])

def run_js_simulation(args):
    (modelname, paramname) = args
    cmd = f"node ./{modelname} {paramname}"
    proc = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)
    try:
        output = proc.communicate()[0]
        output = output.decode("utf-8")
        return fitness_from_tuple(output)
    except:
        # failed
        return 0
    # outputIO = io.StringIO(output)
    # df = pd.read_csv(outputIO, sep="\t", header=None, names=["step","id","type","x","y"])
    # return df

def create_param_files(generation, prefix="", seed=None):
    j = json.loads('''
                    {
                        "conf":{
                            "MAX_ACT": [0, 0, 0, 30],
                            "V": [0, 30, 0, 500],
                            "P": [0, 5, 0, 260],
                            "LAMBDA_ACT": [0, 0, 0, 300],
                            "LAMBDA_V": [0, 1000, 0, 5],
                            "LAMBDA_P": [0, 1, 0, 2],
                            "LAMBDA_CH": [0, 0, 0, 500],
                            "seed": 1
                        }
                    }''')
    paramnames = []
    global POTTS_SEED
    if seed is None:
        j["conf"]["seed"] = POTTS_SEED
        seed = POTTS_SEED
        print("seed:", POTTS_SEED)
    else:
        j["conf"]["seed"] = seed
        print("seed:", POTTS_SEED, "using previous seed for run 2:", seed)
    directory  =f"{PARAM_DIR}/gen{GENERATION_NO}/"
    os.makedirs(directory, exist_ok=True)
    for i, cell in enumerate(generation):
        j["conf"]["MAX_ACT"][-1] = cell["MAX_ACT"]
        j["conf"]["V"][-1] = cell["V"]
        j["conf"]["P"][-1] = cell["P"]
        j["conf"]["LAMBDA_ACT"][-1] = cell["LAMBDA_ACT"]
        j["conf"]["LAMBDA_V"][-1] = cell["LAMBDA_V"]
        j["conf"]["LAMBDA_P"][-1] = cell["LAMBDA_P"]
        j["conf"]["LAMBDA_CH"][-1] = cell["LAMBDA_CH"]
        filename = f"{directory}/{prefix}{i}.json"
        paramnames.append(filename)
        with open(filename, "w+") as f:
            f.write(json.dumps(j))
    return paramnames

def simulate_two(obj):
    # obj: (modelname, (param1, param2))
    (modelname, (param1, param2)) = obj
    return (run_js_simulation((modelname, param1))+run_js_simulation((modelname, param2)))/2

def simulate_generation(generation, modelname, num_procs=12):
    global POTTS_SEED
    paramnames = create_param_files(generation)
    # Also run with previous seed to smooth out errors
    paramnames2 = create_param_files(generation, "seed2_", POTTS_SEED-1)
    if CHANGE_POTTS_SEED_PER_GEN:
        POTTS_SEED = POTTS_SEED + 1
    paramnametuple = zip(paramnames, paramnames2)
    args = list(map(lambda tup: (modelname, tup), paramnametuple))
    with Pool(num_procs) as p:
        sim_results = p.map(simulate_two, args)
    fitnesses = sim_results
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

    print(gen_w_f)
    fitnesses = list(map(lambda x: x[1], gen_w_f))
    print(f"Genfitness min: ", np.min(fitnesses), " mean: ", np.mean(fitnesses), " median: ", np.median(fitnesses), " max: ", np.max(fitnesses), " std dev: ", np.std(fitnesses))
    gen_w_f = list(map(lambda x: x[0], gen_w_f))
    elites = gen_w_f[:N_ELITE]

    i = 0
    # while len(gen_w_f) < GENERATION_SIZE:
    #     gen_w_f.append(mutate_cell(gen_w_f[i%N_ELITE]))
    #     i += 1
    sample_weights = np.array(fitnesses)
    sample_weights = sample_weights - np.min(sample_weights) + 50
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

def init_individual():
    start = {'MAX_ACT': 2, 'P': 250, 'V': 500, 'LAMBDA_ACT': 5, 'LAMBDA_P': 2, 'LAMBDA_V': 5, 'LAMBDA_CH': 5}
    for p in ['MAX_ACT', 'LAMBDA_ACT', 'LAMBDA_CH']:
        start[p] = np.round(np.random.uniform(low=2, high=50), decimals=2)
    return start

def evolve(filename, num_generations, seed=None):
    global POTTS_SEED
    if seed is not None:
        POTTS_SEED = int(seed)
    print(f"Starting to simulate for {filename}, {num_generations}")
    print(f"Starting seed: {POTTS_SEED}")
    init()
    generation = [init_individual() for i in range(GENERATION_SIZE)]
    #generation = list(map(mutate_cell, generation))
    for i in range(num_generations):
        global GENERATION_NO
        print(f"Simulation generation: {GENERATION_NO}")
        gen_fitnesses = simulate_generation(generation, filename, num_procs=cpu_count())
        GENERATION_NO = GENERATION_NO+1
        generation = next_generation_elitism_and_roulette(gen_fitnesses)
