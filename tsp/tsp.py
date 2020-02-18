import numpy as np
from numba import njit
import random
import matplotlib.pyplot as plt
CORSSRATE=0.5
MUTATIONRATE=0.4
GENERATION=200
GENE_SIZE=1000
POPULATION_SIZE=30
SUPERGENERATE=0.5
a = np.loadtxt('tsp.txt')
CITY_POSITION = np.reshape(a, (1000, 2))
DOT=[]

def make_population():
    population=np.vstack([np.random.permutation(GENE_SIZE) for _ in range(POPULATION_SIZE)])
    for i in range(0,len(population),2):
        population[i]=super_man(population[i])
    return population

def super_man(individual):
    index1 = random.randint(0, GENE_SIZE - 100)  # randomly get the mutation position #
    #index2 = random.randint(index1, GENE_SIZE - 1)
    index2=index1+99
    crossgene = []
    superman=[]
    crossgene.append(individual[index1])
    for i in range(0,index2-index1):
        crossgene.append(individual[nearest(crossgene,individual,crossgene[i])])
    for flag, item in enumerate(individual):
        if flag == index1:
            superman.extend(crossgene)
        if item not in crossgene:
            superman.append(item)
    return superman

def super_man1(individual):
    index1 = random.randint(0, GENE_SIZE - 6)  # randomly get the mutation position #
    #index2 = random.randint(index1, GENE_SIZE - 1)
    index2=index1+5
    crossgene = []
    superman=[]
    crossgene.append(individual[index1])
    for i in range(0,index2-index1):
        crossgene.append(individual[nearest(crossgene,individual,crossgene[i])])
    for flag, item in enumerate(individual):
        if flag == index1:
            superman.extend(crossgene)
        if item not in crossgene:
            superman.append(item)
    return superman
@njit
def nearest(crossgene,individual,i):
    best_distance=100
    best_index=-1
    for index,item in enumerate(individual):
        if item not in crossgene:
            distance = ((CITY_POSITION[individual[index]][0] - CITY_POSITION[i][0]) ** 2 +
                        (CITY_POSITION[individual[index]][1] - CITY_POSITION[i][1]) ** 2) ** 0.5
            if distance < best_distance:
                best_distance = distance
                best_index = index
    return best_index

#@njit
def fitness(population,y):
    best_individual=population[0]
    best_fitness=getfitness(best_individual)
    fitness_list=[0]*POPULATION_SIZE
    sumfitness=0
    for i,individual in enumerate(population):
        fitness_list[i]=getfitness(individual)
        plt.scatter(y,1000/fitness_list[i] ,marker='.',s=1)
        plt.xlabel('generation')
        plt.ylabel('distance')
        if fitness_list[i]>best_fitness:
            best_individual=individual
            best_fitness=fitness_list[i]
        sumfitness+=fitness_list[i]
    return sumfitness,best_individual,fitness_list,best_fitness
@njit
def getfitness(individual):
    distance=0.000000
    for i in range(1,len(individual)):
        distance += ((CITY_POSITION[individual[i]][0] - CITY_POSITION[individual[i-1]][0]) ** 2 +
                     (CITY_POSITION[individual[i]][1] - CITY_POSITION[individual[i-1]][1]) ** 2)**0.5

    return 1000/distance
#@njit
def cross(mom,dad):
    index1 = random.randint(0, GENE_SIZE - 3)  # randomly get the mutation position #
    index2 = random.randint(index1+1, GENE_SIZE - 1)
    crossgene = []
    for i in range(index1, index2):
        crossgene.append(mom[i])
    child=[]
    for flag,item in enumerate(dad):
        if flag==index1:
            child.extend(crossgene)
        if item not in crossgene:
            child.append(item)
    return child

def mutaion(individual):
    index = random.sample(range(0, GENE_SIZE), 2)
    new_individual=individual[:]
    new_individual[index[0]],new_individual[index[1]]=new_individual[index[1]],new_individual[index[0]]
    return new_individual
def select(population,sumfitness):
    flag=random.uniform(0,sumfitness)
    for individual in population:
        flag-=getfitness(individual)
        if flag<=0:
            return individual

def born(population,sumfitness):
    dad=select(population,sumfitness)
    while dad is None:
        dad = select(population, sumfitness)
    mom = select(population, sumfitness)
    while dad is None:
        mom=select(population,sumfitness)
    rate = random.uniform(0, 1)
    if rate < CORSSRATE:
        child = cross(dad, mom)
    else:
        child = dad
    rate = random.uniform(0, 1)
    if rate < MUTATIONRATE:
        child = mutaion(child)

    return child
#@njit
def generation(population,y):
    sumfitness,best_individual,fitness_list,best_fitness=fitness(population,y)
    new_population=[]
    new_population.append(best_individual)
    while len(new_population) < POPULATION_SIZE:
        rate = random.uniform(0, 1)
        if rate<0.3:
            new_population.append(super_man(born(population,sumfitness)))
        else:
            new_population.append(super_man1(born(population, sumfitness)))
    return new_population
    #generation_times += 1

def draw(individual):
    x=[0]*(GENE_SIZE+1)
    y=[0]*(GENE_SIZE+1)
    for i in range(0,GENE_SIZE):
        x[i]=CITY_POSITION[individual[i]][0]
    x[GENE_SIZE]=x[0]
    for i in range(0,GENE_SIZE):
        y[i]=CITY_POSITION[individual[i]][1]
    y[GENE_SIZE] = y[0]
    fig_short = plt.figure()
    ax1 = fig_short.add_subplot(1, 1, 1)
    ax1.plot(x, y, lw=0.5, marker='o', markersize=1, mfc='w')
    ax1.set_title('shortest path(EA with superman)')
    plt.show()
if __name__=='__main__':
    population=make_population()
    i=GENERATION
    shortest_distance=[]
    while i>0:
        sumfitness, best_individual, fitness_list,best_fitness=fitness(population,GENERATION-i)
        print("distance:%f,  generation:%d"%(1000/best_fitness,GENERATION-i))
        shortest_distance.append(1000/best_fitness)
        population=generation(population,GENERATION-i)
        i-=1
    draw(best_individual)
    np.savetxt('best_individual.txt',best_individual)
    np.savetxt('short_distance.txt', shortest_distance)
    fig_shortest = plt.figure()
    ax1 = fig_shortest.add_subplot(1, 1, 1)
    times=[i for i in range(GENERATION)]
    ax1.plot(times, shortest_distance)
    ax1.set_title('shortest path(EA with superman)')
    ax1.set_xlabel('Evolutions')
    ax1.set_ylabel('Distance')
    plt.show()
