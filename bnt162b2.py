from dnachisel import *
import python_codon_tables
import mmap

codon_usage_table = python_codon_tables.get_codons_table(57486)
species = 'h_sapiens'
#TaxIDs: 10090,57486 (mus musculus) 9544 (macaca mulatta) 9606 (homo sapiens)

def mmap_io(filename):
    with open(filename, mode="r", encoding="utf8") as file_obj:
        with mmap.mmap(file_obj.fileno(), length=0, access=mmap.ACCESS_READ) as mmap_obj:
            return mmap_obj.read().decode("UTF-8").split("\n")[:-1]

input = mmap_io("codons.txt")

virus,vaccine = input[0],input[1]

def compute_match(one, two):
    num_mismatches = 0
    for base1, base2 in zip(one, two):
        if base1 != base2:
            num_mismatches = num_mismatches + 1
    return (1 - float(num_mismatches) / float(len(one)))

def optimize_virus():
    problem = DnaOptimizationProblem(
        sequence=virus,
        constraints=[
            EnforceTranslation(genetic_table='Standard', start_codon='ATG'),
            EnforceGCContent(mini=0.54, maxi=0.93, window=120),
        ],
        objectives=[
            CodonOptimize(method = "use_best_codon", codon_usage_table=codon_usage_table),
        ]
    )
    problem.resolve_constraints()
    problem.optimize()
    return compute_match(vaccine, problem.sequence)

def average_match(num_iterations=20):
    return max([optimize_virus() for i in range(num_iterations)])

score = average_match()
print(f"{score:.2%}")
