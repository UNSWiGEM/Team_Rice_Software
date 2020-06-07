import create_thermometers

import random
bases = ['A', 'C', 'T', 'G']

# Not sure what to do for the parameters :(

context_before = ''.join([random.choice(bases) for i in range(20)])
context_after = ''.join([random.choice(bases) for i in range(38)])
variable_region = 'AC'

ngen = 1000
cxpb = 0.1
mutpb = 0.1
otmut = 0.1

upper_temp = 28
lower_temp = -5

rbs_start = 5
rbs_end = 7

seed = 12345


create_thermometers.run(context_before, context_after, variable_region,
    ngen, cxpb, mutpb, otmut,
    upper_temp, lower_temp,
    rbs_start, rbs_end,
    seed)
