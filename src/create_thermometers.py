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

VERSION = "1.0.0"

# These must be created before the other functions can be called. Otherwise each "thread" SCOOP starts 
# complains that they cannot find these objects.

creator.create("FitnessMulti", base.Fitness, weights = (1.0, 1.0))
creator.create("Thermometer", list, fitness = creator.FitnessMulti)

def run(context_before, context_after, variable_region,
        NGEN, CXPB, MUTPB, OTMUT,
        upper_temp, lower_temp,
        rbs_start, rbs_end,
        seed):
    """

    With the given parameters, runs a genetic algorithm to increase the change in conformation and binding of a ribosome binding site of a RNA sequence between two given
    temperatures. This should be an acceptable proxy for the performance of an RNA thermometer. 

    Components of the genetic algorithm are created using Distributed Evolutionary Algorithms in Python (DEAP). The genetic algorithm runs for a number of generations, 
    generating an initial population by point mutating each base of the variable region given to this function. It performs a two-point crossover between each sequence 
    of a population, randomly mutates some of them, then evaluates the fitness function for each one. Duplicate sequences are removed. The 20 best from each population are kept, 
    and a new 50 are chosen through automatic epsilon lexicase selecting as described in the DEAP documentation.

    Evaluation of the fitness function is parallelized using the Python 3 library Scalable COncurrent Operations in Python (SCOOP)

    Please read DEAP documentation to understand the use of the genetic algorithm tools:

    https://deap.readthedocs.io/en/master/

    Parameters:

        context_before - a string containing "A", "T", "G", or "C" representing the 5' (upstream) context of a ribosome binding site (RBS)
        context_after - a string containing "A", "T", "G", or "C" representing the 3' (downstream) context of the RBS
        variable_region - a string containing "A", "T", "G", or "C" representing the variable_region, which will be placed upstream of the 3' context

        NGEN - an integer representing the number of generations to evolve the RNA thermometers for
        CXPB - a float between 0 and 1 representing the probability of a two point crossover occuring
        MUTPB - an float between 0 and 1 representing the probability a shuffle mutation occurs (bases are randomly shuffled)
        OTMUT - an float between 0 and 1 representing the probability a point mutation occurs

        upper_temp - an integer representing the temperature in °C for the temperature at which (ideally) the RBS should be fully exposed
        lower_temp - an integer representing the temperature in °C for the temperature at which (ideally) the RBS should be fully bound by the variable_region

        rbs_start - an integer representing the the starting position in the full sequence (context_before + variable_region + context_after) the RBS is located. 
                    Use the location of the beginning of context_after if you're not sure where the RBS is located.

        rbs_end - an integer representing the the ending position in the full sequence (context_before + variable_region + context_after) the RBS is located.
                  Use the location of the end of context_after if you're not sure where the RBS is located.
        
        seed - an integer corresponding to the seed fed to the random number generator.

    Outputs:

        Doesn't return anything.

        Prints the best and worst sequence of each generation as well as their fitness functions.

        Writes to the following files:
            history_<seed>_<NGEN>_<VERSION>.txt - a text file containing a record of the best and worst thermometers (and their fitness functions)
            reg_<seed>_<NGEN>_<VERSION>.obj - a pickled (binary) Python object representing the population at the end of <NGEN> generations.
            reg_<seed>_<NGEN>_<VERSION>.txt - a text file containing the population at the end of <NGEN> generations.
    """
    random.seed(seed)

    toolbox = base.Toolbox()

    toolbox.register("thermometer", initThermometer, creator.Thermometer)
    toolbox.register("newPop", initPopulation, list, toolbox.thermometer, permute_var, bases = bases)

    toolbox.register("mate", tools.cxTwoPoint)
    toolbox.register("mutate", tools.mutShuffleIndexes, indpb = 0.3)
    toolbox.register("select", tools.selAutomaticEpsilonLexicase)
    toolbox.register("previous", tools.selBest)
    toolbox.register("evaluate", evaluate, context_before = context_before, context_after = context_after,
                                           cmd_upper = get_cmd(upper_temp), cmd_lower = get_cmd(lower_temp),
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
                elif change > MUTPB + (1 - (MUTPB + OTMUT)):
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
    """
    
    Evaluates the fitness function for a given individual. This consists of the RBS base pairing probability difference between the two temperatures, evaluated using
    NuPack, as well as the change in secondary structure, evaluated by finding the Levenshtein distance between the dot parentheses representing of the minimum free energy (MFE) 
    secondary structure.

    Generates a randomly named NuPack file, then runs a NuPack command on it for the two temperature bounds.

    Parameters:

        individual - a string representing the individual whose fitness is evaluated

        context_before - a string containing "A", "T", "G", or "C" representing the 5' (upstream) context of a ribosome binding site (RBS)
        context_after - a string containing "A", "T", "G", or "C" representing the 3' (downstream) context of the RBS

        cmd_upper - a string representing NuPack command to run for the upper temperature bound
        cmd_lower - a string representing NuPack command to run for the lower temperature bound

        rbs_start - an integer representing the the starting position in the full sequence (context_before + variable_region + context_after) the RBS is located. 
                    Use the location of the beginning of context_after if you're not sure where the RBS is located.

        rbs_end - an integer representing the the ending position in the full sequence (context_before + variable_region + context_after) the RBS is located.
                  Use the location of the end of context_after if you're not sure where the RBS is located.

    Outputs:

        Returns:
            pair_diff - a float representing the expected value of the base pairing probability difference between the two temperature bounds
            distance - an integer represeting the Levenshtein distance between the string representing the secondary structure at the upper temperature bound
                       and the string representing the secondary structure at the lower temperature bound.

    """
    sequence = context_before + "".join(individual) + context_after
    assert len(sequence) == 60
    # Generates unique file name for NuPack
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

    # Runs NuPack command at terminal

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

def rbs_prob(ppairs_lines, rbs_start, rbs_end):
    """

    Produces the expected value of the number of bases bound in the RBS at a given temperature. Parses the output of a NuPack command, extracts the probability
    that each base is bound, and subtracts from that the probability that that base is unbound.

    Parameters:

        ppairs_lines - a list of strings representing the .ocx-ppairs file outputted by the 'complexes' NuPack command
             
        rbs_start - an integer representing the the starting position in the full sequence (context_before + variable_region + context_after) the RBS is located. 
                    Use the location of the beginning of context_after if you're not sure where the RBS is located.

        rbs_end - an integer representing the the ending position in the full sequence (context_before + variable_region + context_after) the RBS is located.
                  Use the location of the end of context_after if you're not sure where the RBS is located.

    Outputs:

        total - a float representing the expected value of the probability that any base of the RBS is bound at a given temperature.

    """
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

def mutate(ind, bases):
    """
    
    Mutates an individual at a single base.

    Parameters:

        ind - a string representing an individual to be mutated

        bases - a list of characters representing the possible bases to substitute

    Outputs:

        ind - a string representing an individual which has undergone a single point mutation

    """
    snp = random.randint(0, len(ind) - 1)
    char = ind[snp]

    while char == ind[snp]:
        char = random.choices(bases)[0]

    ind[snp] = char

    return ind,

def initPopulation(pplr, crtr, var_pop, seq, bases):
    """

    A higher order function which calls a number of helper functions to create an initial population using DEAP from a starting variable region.

    Parameters:

        pplr - a wrapper function for the population
        crtr - a DEAP creator function for the individual members of a population
        var_pop - a creator function that permutes a sequence by mutating each base
        seq - a string representing the starting variable region
        bases - a list of characters representing the DNA (will be transcribed into RNA) bases

    Outputs:

        a list of strings representing the population in such a way that DEAP can understand and manipulate it


    """
    pop = var_pop(seq, bases)
    return pplr(crtr(var = item) for item in pop)

def initThermometer(crtr, var):
    """
    
    Wrapper function to initialize an individual using DEAP.

    Parameters:

        crtr - a DEAP creator function for the individual members of a population
        var - a object (a string in this case) representing the individual to be wrapped by a DEAP function call.

    Outputs:

        thermo - a DEAP-wrapped individual

    """
    thermo = crtr(var)
    return thermo

def extremes(pop):
    """

    Checks through a population (by using DEAP to check which individuals dominate all others) and finds the best
    and worst individuals (in terms of fitness functions)

    Parameters:
        pop - a list of strings representing the population

    Output:
        best, worst - a 2-tuple of strings repesenting the best and worst individuals

    """
    best = pop[0]
    worst = pop[1]
    for ind in pop:
        if ind.fitness.valid and ind.fitness.dominates(best.fitness):
            best = ind
        if ind.fitness.valid and worst.fitness.dominates(ind.fitness):
            worst = ind
    return best, worst

def fmtInd(ind, context_before, context_after):
    """

    Formats an individual by adding the 5' and 3' contexts (and fitnesses) so that it can be printed to the console.
    
    Parameters:

        individual - a string representing the individual whose fitness is evaluated

        context_before - a string containing "A", "T", "G", or "C" representing the 5' (upstream) context of a ribosome binding site (RBS)
        context_after - a string containing "A", "T", "G", or "C" representing the 3' (downstream) context of the RBS

    Outputs:

        a string representing the individual which will be printed to the console.

    """
    return whole(ind, context_before, context_after) + ", " + \
        "{0:.3f}".format(ind.fitness.values[0]) + " base pairs, " + \
        "{0:.3f}".format(ind.fitness.values[1]) + " differences."

def fmtWrite(ind, context_before, context_after):
    """

    Formats an individual by adding the 5' and 3' contexts (and fitnesses) so that it can be written to a text file.
    
    Parameters:

        individual - a string representing the individual whose fitness is evaluated

        context_before - a string containing "A", "T", "G", or "C" representing the 5' (upstream) context of a ribosome binding site (RBS)
        context_after - a string containing "A", "T", "G", or "C" representing the 3' (downstream) context of the RBS

    Outputs:

        a string representing the individual which will be written to a text file.

    """

    fullseq = context_before + "".join(ind) + context_after + ", "
    fitnesses = map(lambda fitness: "{0:.3f}".format(fitness) + ", ", ind.fitness.values)
    return fullseq + "".join(fitnesses).strip()[:-1] + "\n"

def whole(ind, context_before, context_after):
    """

    Parameters:

        individual - a string representing the individual whose fitness is evaluated

        context_before - a string containing "A", "T", "G", or "C" representing the 5' (upstream) context of a ribosome binding site (RBS)
        context_after - a string containing "A", "T", "G", or "C" representing the 3' (downstream) context of the RBS

    Outputs:
        a string representing the individual (has 5' and 3' contexts on each of the variable_region)

    """
    return context_before + "".join(ind) + context_after

def permute_var(complement, bases):
    """

    Given a variable region, returns a list of strings which consist of the variable_region mutated at every point.

    Parameters:
        complement - a string representing the sequence (usually complementary to the RBS) from which a permutation of point mutations is to be generated.

        bases - a list of characters representing the bases.

    """
    nested_result = []

    for index, base in enumerate(complement):
        other_bases = [elem for i, elem in enumerate(bases) if base != elem]
        nested_result.append([complement[:index] + elem + complement[index + 1:] for i, elem in enumerate(other_bases)])

    result = [item for sublist in nested_result for item in sublist]
    result.append(complement)

    return result

def levenshteinDistance(s1, s2):
    """

    Finds the Levenshtein (edit) distance between two strings.

    Parameters:
        s1 - a string
        s2 - another string

    Outputs:

        an integer representing the Levenshtein distance between s1 and s2

    """
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

def get_cmd(temp):
    """
    Given a temperature, returns the command line NuPack command to run

    Parameters:
        temp - an integer representing the temperature at which the NuPack command should be run

    Output:
        a string representing the command to be executed on the console by this program

    """
    return "complexes -T " + str(temp) + " -material rna -pairs -mfe -quiet "
