import random
import sympy
import numpy as np
import matplotlib.pyplot as plt
from math import sin,cos
from numba import njit
GENERATION=5
POPULATION_SIZE=30
GENE_SIZE=256
CROSS_RATE=0.8
MUTATION_RATE=0.3
operator_dic={1:'+',2:'-',3:'*',4:'/',5:'x',6:'c',7:'sin',8:'cos'}
operator=['+','-','*','/']
trangle=['sin','cos']
x_or_const=['x','c']
a = np.loadtxt('data.txt')
POINT_POSITION = np.reshape(a, (1000, 2))
new_position=[]
for i in range(0,999,5):
    new_position.append(POINT_POSITION[i])
new_position = np.array(new_position)
def make_population():
    population=[]
    for _ in range(POPULATION_SIZE):
        individual=generate()
        error1 = 'zoo'
        error2 = 'nan'
        string = str(sympy.expand(to_string(individual, 1)))
        while error1 in string or error2 in string:
            individual = generate()
            string = str(sympy.expand(to_string(individual, 1)))
        population.append(individual)
    return population
def fitness(population,y):
    best_individual = population[0]
    best_fitness = get_fitness(best_individual)
    fitness_list = [0] * POPULATION_SIZE
    sumfitness = 0
    a = []
    for i, individual in enumerate(population):
        fitness_list[i] = get_fitness(individual)


        if fitness_list[i] > best_fitness:
            best_individual = individual
            best_fitness = fitness_list[i]
        sumfitness += fitness_list[i]
    return sumfitness, best_individual, fitness_list, best_fitness,a
'''def is_stright_line(gene):
    gene = str(sympy.expand(to_string(gene, 1)))
    x = 0.01
    y1=eval(gene)
    x=0.02
    y2 = eval(gene)
    x=0.03
    y3 = eval(gene)
    if y1==y2 and y1==y3:
        a=generate()
        return a
    else:
        return gene'''

def cross(mom,dad):
    child1=mom[:]
    child2=dad[:]
    index= random.randint(2,255)# 随机生成突变起始位置 #
    while child1[index]==0 or child2[index]==0:
        index = random.randint(2, 255)  # 随机生成突变起始位置 #
    child1,child2=change_subtree(index,child1,child2)
    error1 = 'zoo'
    error2 = 'nan'
    child1_gene = str(sympy.expand(to_string(child1, 1)))
    child2_gene = str(sympy.expand(to_string(child2, 1)))
    if (error1 in child1_gene or error2 in child1_gene) or (error1 in child2_gene or error2 in child2_gene) :
        return mom
    else: return child1 if get_fitness(child1)>get_fitness(child2) else child2
def change_subtree(index,childa,childb):
    child1=childa[:]
    child2=childb[:]
    child1[index],child2[index]=child2[index],child1[index]
    i=2
    while 2*index<256:
        index*=2
        for j in range(i):
            child1[index+j],child2[index+j]=child2[index+j],child1[index+j]
        i*=2
    return child1,child2
def mutation(a):
    child = a[:]
    index = random.randint(1, 255)
    while child[index] == 0:
        index = random.randint(1, 255)
    if child[index] in operator:
        rate = random.uniform(0, 1)
        if rate < 0.1:
            child[index] = operator_dic[random.randint(1, 4)]
        elif rate>0.1 and rate<0.2:
            child[index]=operator_dic[random.randint(7, 8)]
            child[2 * index + 1] = 0
            index=2*index+1
            while 2*index+1<=256:
                index=2*index
                child[index]=0
                child[index+1]=0
    elif child[index] in trangle:
        rate = random.uniform(0, 1)
        if rate<0.4:
            child[index]=operator_dic[random.randint(7,8)]

        else:
            child[index]=operator_dic[random.randint(1, 4)]
            child[2*index+1]=str(random.uniform(-10,10))
    elif child[index] == 'x':
        if 2*index<256:
            rate = random.uniform(0, 1)
            if rate<0.5:
                child[index]=child[index]=operator_dic[random.randint(7,8)]
                child[index*2]='x'
        else:
            child[index] = str(random.uniform(-10, 10))
    else:
        rate = random.uniform(0, 1)
        if rate < 0.4:
            child[index] = str(float(child[index]) + 0.1)
        elif rate < 0.8 and rate>=0.4:
            child[index] = str(float(child[index]) - 0.1)
        else:
            rate = random.uniform(0, 1)
            if 2*index<256 and rate<0.5:
                child[index]=operator_dic[random.randint(7,8)]
                child[2*index]='x'
            else:
                child[index] = 'x'
    error1 = 'zoo'
    error2 = 'nan'
    this_gene = str(sympy.expand(to_string(child, 1)))
    if error1 in this_gene or error2 in this_gene:
        return a
    else: return child
def tox(a,index):
    a[index]='x'
    while 2*index<=255:
        index=index*2
        a[index]=0
        a[index+1]=0
    return a
def select(population,sumfitness):
    flag=random.uniform(0,sumfitness)
    for individual in population:
        flag-=get_fitness(individual)
        if flag<=0:
            return individual
def born(population,sumfitness):
    dad=select(population,sumfitness)
    while dad is None:
        dad = select(population, sumfitness)
    mom = select(population, sumfitness)
    while mom is None:
        mom=select(population,sumfitness)
    rate = random.uniform(0, 1)
    if rate < CROSS_RATE:
        child = cross(dad, mom)
    else:
        child = dad
    rate = random.uniform(0, 1)
    if rate < MUTATION_RATE:
        child=mutation(child)
    return child
def generate():
    gene=[0]*256
    gene[1]=operator_dic[random.randint(1,4)]
    for i in range(1,7):
        for j in range(2**i,2**(i+1)):
            if gene[int(j/2)] in operator:
                gene[j]=operator_dic[random.randint(1,8)]
            if gene[int(j/2)] in trangle:
                if j%2==0:
                    gene[j]=operator_dic[random.randint(5,8)]
    for i in range(2**7,2**8):
        if gene[int(i/2)] in trangle:
            if i%2==0:
                gene[i] = operator_dic[random.randint(5, 6)]
        if gene[int(i/2)] in operator:
            gene[i]=operator_dic[random.randint(5,6)]
    error1 = 'zoo'
    error2 = 'nan'
    string = str(sympy.expand(to_string(gene, 1)))
    while error1 in string or error2 in string:
        gene = generate()
        string = str(sympy.expand(to_string(gene, 1)))
    set_c(gene)
    return gene
def generation(population,y):
    sumfitness,best_individual,fitness_list,best_fitness,a=fitness(population,y)
    s = []
    while len(s) < len(a):
        s.append(y)
    plt.scatter(s, a, marker='.', s=1)
    plt.xlabel('generation(30 individuals)')
    plt.ylabel('Mean distance')
    new_population=[]
    new_population.append(best_individual[:])
    while len(new_population) < POPULATION_SIZE:
        #child=born(population,sumfitness)
        #child=is_stright_line(child)
        new_population.append(born(population,sumfitness))
    return new_population

def to_string(gene,x):
    result=""
    if 2*x>=256 or gene[2*x]==0:
        result+=gene[x]
    elif gene[x] in trangle:
        result=result+gene[x]+'('+to_string(gene,2*x)+')'
    else :
        result=result+'('+to_string(gene,2*x)+gene[x]+to_string(gene,2*x+1)+')'
    return result
def set_c(a):
    for index,item in enumerate(a):
        if item=='c' :
            a[index]=str(random.uniform(-10,10))
    return a
#@njit
def get_fitness(gene):
    sum=0
    gene=str(sympy.expand(to_string(gene,1)))
    for i in range(0,len(POINT_POSITION),50):
        x=POINT_POSITION[i][0]
        sum+=abs(POINT_POSITION[i][1]-eval(gene))
    return 100/sum
def draw(shortest_distance):
    generation=[i for i in range(GENERATION)]
    fig_short = plt.figure()
    ax1 = fig_short.add_subplot(1, 1, 1)
    ax1.plot(generation,shortest_distance, lw=0.5, marker='o', markersize=1, mfc='w')
    ax1.set_title('function learning curve(GP)')
    ax1.set_xlabel('Evolutions')
    ax1.set_ylabel('Mean distance')
    plt.show()
def draw_graph(gene):
    point_x=[]
    point_y=[]
    for i in range(len(POINT_POSITION)):
        point_x.append(POINT_POSITION[i][0])
        point_y.append(POINT_POSITION[i][1])
    npgene=''
    for index in range(len(gene)):
        if (gene[index]=='s' and gene[index+1]=='i') or gene[index]=='c':
            npgene+='np.'
        npgene+=gene[index]
    print(npgene)
    x = np.arange(0, 10, 0.005)
    npgene+='+(x-x)'
    print(npgene)
    y = eval(npgene)
    plt.xlim(0, 10)
    plt.plot(x, y,c='r')
    plt.scatter(point_x,point_y,marker='o',s=1)
    plt.title('Curve fitting(GP)')
    plt.show()
if __name__=="__main__":
    population = make_population()
    i = GENERATION
    shortest_distance = []
    while i > 0:
        sumfitness, best_individual, fitness_list, best_fitness,a = fitness(population,GENERATION-1)
        print(str(sympy.expand(to_string(best_individual, 1))))
        print("mean_distance:%f,  generation:%d" % (100/best_fitness/20,GENERATION-i))
        shortest_distance.append(100/best_fitness/20)
        population = generation(population, GENERATION - i)
        i -= 1
    best_string = str(sympy.expand(to_string(best_individual, 1)))
    draw_graph(best_string)
    #np.savetxt('best_individual.txt', best_individual)
    #np.savetxt('short_distance.txt', shortest_distance)
    draw(shortest_distance)
    #np.savetxt('ea data1.txt', shortest_distance)
