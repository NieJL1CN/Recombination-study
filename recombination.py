import numpy as np
import random

def generate_f1_genotype():
    return ['AB'] * 10

def recombine(chromosome1, chromosome2, cM=10):
    # 每10cM的距离均匀分布一个标记位点
    breakpoints = [i for i in range(1, 10)]
    # 发生重组的位点
    recomb_positions = [p for p in breakpoints if random.random() < cM / 100.0]
    
    new_chromosome1, new_chromosome2 = [], []
    current_chromosome = 0  # 0表示chromosome1，1表示chromosome2
    for i in range(10):
        if i in recomb_positions:
            current_chromosome = 1 - current_chromosome
        if current_chromosome == 0:
            new_chromosome1.append(chromosome1[i])
            new_chromosome2.append(chromosome2[i])
        else:
            new_chromosome1.append(chromosome2[i])
            new_chromosome2.append(chromosome1[i])
    
    return new_chromosome1, new_chromosome2

def self_cross(f1_genotype, generations=10, progeny_size=300):
    population = [f1_genotype] * progeny_size
    all_generations = [population]
    
    for gen in range(generations):
        next_gen = []
        for individual in population:
            chromosome1 = [allele[0] for allele in individual]
            chromosome2 = [allele[1] for allele in individual]
            for _ in range(progeny_size):
                recomb_chromosome1, recomb_chromosome2 = recombine(chromosome1, chromosome2)
                new_individual = [recomb_chromosome1[i] + recomb_chromosome2[i] for i in range(10)]
                next_gen.append(new_individual)
        population = next_gen
        all_generations.append(population)
    
    return all_generations

# 初始F1代基因型
f1_genotype = generate_f1_genotype()
# 模拟10代的自交过程
generations = self_cross(f1_genotype, generations=10, progeny_size=300)
# 打印每代的基因型
for gen_num, generation in enumerate(generations):
    print(f"Generation {gen_num}:")
    for individual in generation[:5]:  # 只显示前5个个体
        print(individual)
    print()