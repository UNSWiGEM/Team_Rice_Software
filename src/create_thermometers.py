from deap import base, creator, tools
from scoop import futures
from functools import reduce
import uuid
import subprocess
import os
import random
import pickle
import operator
import time

bases = ['A', 'T', 'G', 'C']

VERSION = "0.1.0"

creator.create("FitnessMulti", base.Fitness, weights = (1.0, 1.0))
creator.create("Thermometer", list, fitness = creator.FitnessMulti)

def run(context_before, context_after, variable_region,
        NGEN, CXPB, MUTPB, OTMUT,
        upper_temp, lower_temp,
        rbs_start, rbs_end,
        seed):
    random.seed(seed)

    toolbox = base.Toolbox()

    toolbox.register("thermometer", initThermometer, creator.Thermometer)
    toolbox.register("newPop", initPopulation, list, toolbox.thermometer, permute_var, bases = bases)

    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutShuffleIndexes, indpb = 0.3)
    toolbox.register("select", tools.selAutomaticEpsilonLexicase)
    toolbox.register("previous", tools.selBest)
    toolbox.register("evaluate", evaluate, context_before = context_before, context_after = context_after,
                                           cmd_upper = get_cmd_upper(upper_temp), cmd_lower = get_cmd_lower(lower_temp),
                                           rbs_start = rbs_start, rbs_end = rbs_end)

    toolbox.register("map", futures.map)

    pop = toolbox.newPop(seq = variable_region)
    fitness = toolbox.map(toolbox.evaluate, pop)
    for ind, fit in zip(pop, fitness):
        ind.fitness.values = fit

    best, _ = extremes(pop)
    print("---INITIAL BEST INDIVIDUAL---\n" + fmtInd(best, context_before, context_after))

    with open("history_reg_" + str(seed) + "_" + str(NGEN) + "_" + str(VERSION) + ".txt", "x") as history:
        history.write("VERSION = REG_" + str(VERSION) + "\n")
        history.write("SEED = " + str(seed) + "\n")
        history.write("NGEN = " + str(NGEN) + "\n")
        history.write("VARIABLE_START = " + variable_region + "\n")
        for g in range(NGEN):
            prevBest = toolbox.previous(pop, 20)
            chosen = toolbox.select(pop, 50)
            older = list(map(lambda ind: toolbox.clone(ind), prevBest))
            offspring = list(map(lambda ind: toolbox.clone(ind), chosen))

            best, worst = extremes(offspring + older)

            print("---BEST INDIVIDUAL---\n" + fmtInd(best, context_before, context_after))
            print("---WORST INDIVIDUAL---\n" + fmtInd(worst, context_before, context_after))
            print("Generation #" + str(g) + " completed!")

            history.write(fmtWrite(best, context_before, context_after))

            for child1, child2 in zip(offspring[::2], offspring[1::2]):
                if random.random() < CXPB:
                    toolbox.mate(child1, child2)
                    toolbox.mate(child1, child2)
                    del child1.fitness.values
                    del child2.fitness.values

            for mutant in offspring:
                change = random.random()
                if change < MUTPB:
                    toolbox.mutate(mutant)
                    del mutant.fitness.values
                elif change > MUTPB + OTMUT:
                    mutate(mutant, bases)
                    del mutant.fitness.values

            offspring = offspring + older

            invalid_ind = [ind for ind in offspring if not ind.fitness.valid]
            fitnesses = toolbox.map(toolbox.evaluate, invalid_ind)
            for ind, fit in zip(invalid_ind, fitnesses):
                ind.fitness.values = fit

            no_dup = []
            for ind in offspring:
                if ind not in no_dup:
                    no_dup.append(ind)

            pop[:] = no_dup

    best, _ = extremes(pop)
    print("---FINAL BEST INDIVIDUAL---\n" + fmtInd(best, context_before, context_after))

    with open("reg_" + str(seed) + "_" + str(NGEN) + "_" + str(VERSION) + ".obj", "xb") as binFile:
        pickle.dump(pop, binFile)

    with open("reg_" + str(seed) + "_" + str(NGEN) + "_" + str(VERSION) + ".txt", "x") as savefile:
        savefile.write("VERSION = REG_" + str(VERSION) + "\n")
        savefile.write("SEED = " + str(seed) + "\n")
        savefile.write("NGEN = " + str(NGEN) + "\n")
        savefile.write("VARIABLE_START = " + variable_region + "\n")
        for ind in pop:
            savefile.write(fmtWrite(ind, context_before, context_after))

def evaluate(individual, context_before, context_after, cmd_upper, cmd_lower, rbs_start, rbs_end):
    sequence = context_before + "".join(individual) + context_after
    assert len(sequence) == 60
    name = str(uuid.uuid4())
    input_name = name + ".in"
    output_temp_name = name + ".ocx-mfe"
    output_pair_name = name + ".ocx-ppairs"
    output_other1 = name + ".ocx"
    output_other2 = name + ".ocx-epairs"
    output_other3 = name + ".ocx-key"

    up_cmd = cmd_upper + name
    lo_cmd = cmd_lower + name

    with open(input_name, "x") as stage:
        stage.write("1\n" + sequence + "\n1\n")

    subprocess.run(up_cmd, shell=True)
    with open(output_temp_name, "r") as temp, open(output_pair_name, "r") as pair:
        res_t = temp.readlines()
        res_p = pair.readlines()
    up_erg = [line for line in res_t if "%" not in line][2].strip()
    up_lev = [line for line in res_t if "%" not in line][3].strip()
    up_pair = rbs_prob(res_p, rbs_start, rbs_end)

    os.remove(output_temp_name)
    os.remove(output_pair_name)
    os.remove(output_other1)
    os.remove(output_other2)
    os.remove(output_other3)

    subprocess.run(lo_cmd, shell=True)
    with open(output_temp_name, "r") as temp, open(output_pair_name, "r") as pair:
        res_t = temp.readlines()
        res_p = pair.readlines()
    low_erg = [line for line in res_t if "%" not in line][2].strip()
    low_lev = [line for line in res_t if "%" not in line][3].strip()
    low_pair = rbs_prob(res_p, rbs_start, rbs_end)

    os.remove(output_temp_name)
    os.remove(output_pair_name)
    os.remove(output_other1)
    os.remove(output_other2)
    os.remove(output_other3)
    os.remove(input_name)

    delta = float(up_erg) - float(low_erg)
    pair_diff = low_pair - up_pair
    distance = levenshteinDistance(up_lev, low_lev)

    return pair_diff, distance

def mutate(ind, bases):
    snp = random.randint(0, len(ind) - 1)
    char = ind[snp]

    while char == ind[snp]:
        char = random.choices(bases)[0]

    ind[snp] = char

    return ind,

def initPopulation(pplr, crtr, var_pop, seq, bases):
    pop = var_pop(seq, bases)
    return pplr(crtr(var = item) for item in pop)

def initThermometer(crtr, var):
    thermo = crtr(var)
    return thermo

def extremes(pop):
    best = pop[0]
    worst = pop[1]
    for ind in pop:
        if ind.fitness.valid and ind.fitness.dominates(best.fitness):
            best = ind
        if ind.fitness.valid and worst.fitness.dominates(ind.fitness):
            worst = ind
    return best, worst

def fmtInd(ind, context_before, context_after):
    return whole(ind, context_before, context_after) + ", " + \
        "{0:.3f}".format(ind.fitness.values[0]) + " base pairs, " + \
        "{0:.3f}".format(ind.fitness.values[1]) + " differences."

def fmtWrite(ind, context_before, context_after):
    fullseq = context_before + "".join(ind) + context_after + ", "
    fitnesses = map(lambda fitness: "{0:.3f}".format(fitness) + ", ", ind.fitness.values)
    return fullseq + "".join(fitnesses).strip()[:-1] + "\n"

def whole(ind, context_before, context_after):
    return context_before + "".join(ind) + context_after

def rbs_prob(ppairs_lines, rbs_start, rbs_end):
    n_plus_one = rbs_end + 1
    pairs_lines = [line.strip() for line in ppairs_lines if "%" not in line and line != "\n"][1:]

    pairings = []
    for line in pairs_lines:
        columns = line.split("\t")
        columns = tuple([item.strip() for item in columns][:3])

        pairings.append(columns)

    pairings = [(int(tup[0]), int(tup[1]), float(tup[2])) for tup in pairings]

    ispairings = [tup for tup in pairings if (tup[0] >= rbs_start and tup[0] <= rbs_end) or (tup[1] >= rbs_start and tup[1] <= rbs_end)]
    ispairings = [tup for tup in ispairings if (int(tup[0]) != n_plus_one) and (int(tup[1]) != n_plus_one)]

    unpairings = [tup for tup in pairings if tup[1] == n_plus_one and (tup[0] >= rbs_start and tup[0] <= rbs_end)]

    prob = reduce(operator.add, [tup[2] for tup in ispairings])
    unprob = reduce(operator.add, [tup[2] for tup in unpairings])

    total = prob - unprob

    return total

def permute_var(complement, bases):
    nested_result = []

    for index, base in enumerate(complement):
        other_bases = [elem for i, elem in enumerate(bases) if base != elem]
        nested_result.append([complement[:index] + elem + complement[index + 1:] for i, elem in enumerate(other_bases)])

    result = [item for sublist in nested_result for item in sublist]
    result.append(complement)

    return result

def levenshteinDistance(s1, s2):
    if len(s1) > len(s2):
        s1, s2 = s2, s1

    distances = range(len(s1) + 1)
    for i2, c2 in enumerate(s2):
        distances_ = [i2+1]
        for i1, c1 in enumerate(s1):
            if c1 == c2:
                distances_.append(distances[i1])
            else:
                distances_.append(1 + min((distances[i1], distances[i1 + 1], distances_[-1])))
        distances = distances_
    return distances[-1]

def get_cmd_upper(upper_temp):
    return "complexes -T " + str(upper_temp) + " -material rna -pairs -mfe -quiet "

def get_cmd_lower(lower_temp):
    return "complexes -T " + str(lower_temp) + " -material rna -pairs -mfe -quiet "
