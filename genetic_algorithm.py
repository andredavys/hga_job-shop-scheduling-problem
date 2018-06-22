import numpy as np
import random
import time
import operator
import math
import glob

def buildMatrices(file_name):
	count_jobs = 0
	count_mach = 0
	jobs_list  = []
	machines = []
	processing_times = []
	with open(file_name, 'r') as file:
		i = 0 # index file
		j = 0 # index job
		for line in file:
			line_split = line.split()
			if i == 0:
				count_jobs = int(line_split[0])
				count_mach = int(line_split[1])
				print(count_jobs, count_mach)
				machines = np.zeros((count_jobs, count_mach))
				processing_times = np.zeros((count_jobs, count_mach))
			else:
				im = 0 # index machine
				ip = 0 # index processing time
				for k in range(len(line_split)):
					if k%2 == 0:
						machines[j,im] = line_split[k]
						im += 1
					else:
						processing_times[j,ip] = line_split[k]
						ip += 1
				j += 1
			i += 1

	return machines, processing_times, count_jobs, count_mach

def evaluate(chromossome, machines, processing_times, count_jobs, count_mach):
	m_start = np.zeros(count_mach)
	j_start = np.zeros(count_jobs)
	t_next  = np.zeros(count_jobs) # control number of operations
	s_next  = np.ones(count_mach)
	
	for k in range(len(chromossome)):
		# print('+++++++++ iteraçao +++++++++', k)
		i = chromossome[k]
		l = int(t_next[i])
		j = int(machines[i, int(t_next[i])])
		# print('job', i, 'machine', j)
		t_next[i] += 1
		s_next[j] += 1
		# print('t_next', t_next)
		# print('s_next', s_next)
		start = max(j_start[i], m_start[j])
		# print('start', start)
		# print('processing', processing_times[i][l])
		j_start[i] = start + processing_times[i][l]
		m_start[j] = start + processing_times[i][l]
		# print('j start\n', j_start)
		# print('m start\n', m_start)

	# print('job start', j_start)
	# print('machine start', m_start)
	return max(max(j_start), max(m_start))

def build_initial_population(size, count_jobs, count_mach):
	population = []
	for i in range(size):
		population.append(build_random_chromossome(count_jobs, count_mach))
	return population

def calculate_index(p, count_jobs):
	p_index = []
	count_task = [0 for j in range(count_jobs)]
	for i in p:
		p_index.append((i, count_task[i]))
		count_task[i] += 1
	return p_index

def decode_chromossome(genotype):
	phenotype = []
	for task, idx in genotype:
		phenotype.append(task)
	return phenotype


def general_order_crossover(p1, p2, count_jobs):
	p1_index = calculate_index(p1, count_jobs)
	p2_index = calculate_index(p2, count_jobs)
	n = len(p1_index) #total of tasks (lenght of chromossome)
	tam = random.randint(int(n/3), int(n/2))
	j = random.randint(0, n-tam)
	# print(j, tam)
	genes_p2 = p2_index[j:j+tam]
	# print(genes_p2)
	child = []
	index_to_insert = 0

	new_list = [x for x in p1_index if x not in set(genes_p2)] # take p1 tasks to insert it

	# verify where I should insert p2 tasks
	for gene_idx in p1_index:
		if gene_idx not in genes_p2:
			index_to_insert+=1
		else:
			break

	#generate child with taks of two parents
	child =  new_list[0:index_to_insert] + genes_p2 + new_list[index_to_insert:]
	
	# print('index', index_to_insert)
	# print('parent 1', p1_index)
	# print('parent 2', p2_index)
	# print('gene_add', genes_p2)
	# print('new_list', new_list)
	return child

def build_random_chromossome(count_jobs, count_mach):
	chromossome = []
	to_insert = []
	for i in range(count_mach*count_jobs):
		to_insert.append(i%count_jobs)
	while len(to_insert) > 0:
		index = random.randint(0,len(to_insert)-1)
		value = to_insert.pop(index)
		chromossome.append(value)
	return chromossome

def swap_two_genes(ind, k, machines, processing_times, count_jobs, count_mach, dict_chromossomes):
	best_value = evaluate(ind, machines, processing_times, count_jobs, count_mach)
	best_ind = [i for i in ind] 
	new_ind = [i for i in ind]
	actual_value = 99999
	for i in range(k-1, max(0,k - count_jobs*count_mach*0.3), -1):
		if ind[i] != ind[k]:
			new_ind [k], new_ind[i] = new_ind[i], new_ind[k]
			if tuple(new_ind) not in dict_chromossomes:
				actual_value = evaluate(new_ind, machines, processing_times, count_jobs, count_mach)
				dict_chromossomes[tuple(new_ind)] = actual_value
			else:
				actual_value = dict_chromossomes[tuple(new_ind)]
			i += 1
			if actual_value < best_value:
				best_value = actual_value
				best_ind = new_ind
	return best_ind

def local_search(ind, rate, machines, processing_times, count_jobs, count_mach, dict_chromossomes):
		k = len(ind)
		# TODO: backtracking to local search calls

		new_ind_aux = [i for i in ind] # copy by value
		# for j in range(k-1, k-math.ceil(k*rate)-1, -1):
		for j in range(k-1, k-1-rate):
			new_ind = swap_two_genes(ind, j, machines, processing_times, count_jobs, count_mach, dict_chromossomes)
		return new_ind_aux

if __name__ == '__main__':
	file_name = 'instances/la35_instance.jsp'

	for file in glob.glob('instances/*.jsp'):
		time_exec = time.time()
		file_results = open('results.txt', 'a')
		file_times   = open('times.txt', 'a')
		file_name = file
		print('tests to file: ', file)
		
		machines, processing_times, count_jobs, count_mach = buildMatrices(file_name)
		# print('matrix machines\n', machines)
		# print('matrix processing times\n', processing_times)
		# chromossome1 = [0,2,1,2,1,1,0,0]
		# chromossome2 = [0,2,1,0,1,1,0,2]
		# a = local_search(chromossome1, 1, machines, processing_times, count_jobs, count_mach)
		# print(evaluate(chromossome1, machines, processing_times, count_jobs, count_mach), 
		# 	evaluate(a,machines, processing_times, count_jobs, count_mach))


		# TODO: test this rates

		n = 200	
		elitism_rate  = 0.05 # rate of individual will survivor to next generation
		survival_rate = 0.7
		random_rate   = 0.25
		lcl_srch_rate = count_jobs*count_mach*0.1
		num_it = 600

		file_results.write('configuração:    '+str(n)+' '+str(num_it)+' '+str(lcl_srch_rate)+'    ')
		file_results.write( str(file)+'     ')
		file_times.write('***********'+ str(file) +'***********\n')
		file_times.write('configuração:   '+str(n)+' '+str(num_it)+' '+str(lcl_srch_rate)+'\n')
		# print(chromossome)
		# print(calculate_index(chromossome, count_jobs))
		# print('child 1:', general_order_crossover(chromossome1, chromossome2, count_jobs))
		# print('child 2:', general_order_crossover(chromossome2, chromossome1, count_jobs))
		sum_total = []
		for aux_k in range(5):
			# ******************************* BUILD INITIAL POPULATION STEP ******************************* 
			t_build_pop = time.time()
			population = build_initial_population(n, count_jobs, count_mach)
			t_build_pop = time.time() - t_build_pop

			n = 200	
			elitism_rate  = 0.05 # rate of individual will survivor to next generation
			survival_rate = 0.7
			random_rate   = 0.25
			lcl_srch_rate = count_jobs*count_mach*0.1
			num_it = 600

			dict_chromossomes = {} # store evaluate to each chromossome (will be global)
			for it in range(num_it):
				print('iteração ', it)
				t_iteration = time.time()
			# ******************************* CROSSOVER STEP *******************************
				childs = []
				t_cross_pop = time.time()
				# print('popular before gox:', len(population))
				for i in range(int(n/2)):
					p1 = population[i]
					p2 = population[i + int(n/2)]
					child = general_order_crossover(p1, p2, count_jobs)
					childs.append(child)
					child = general_order_crossover(p1, p2, count_jobs)
					childs.append(child)
				for ch in childs:
					population.append(decode_chromossome(ch))
				t_cross_pop = time.time() - t_cross_pop
				# print('population after gox: ', len(population))

			# ******************************* SURVIVOR STEP *******************************

				t_survivor = time.time()
				dict_actual_pop = {}
				dict_ind_pop = {}
				i = 0
				for ind in population:
					if tuple(ind) not in dict_chromossomes:
						value = evaluate(ind, machines, processing_times, count_jobs, count_mach)
						dict_chromossomes[tuple(ind)] = value
						dict_actual_pop[i] = value
						# dict_ind_pop[tuple(ind)] = value
					else:
						dict_actual_pop[i] = dict_chromossomes[tuple(ind)]
						# dict_ind_pop[tuple(ind)] = dict_chromossomes[tuple(ind)]
					i += 1

				sorted_actual_pop = sorted(dict_actual_pop.items(), key=operator.itemgetter(1))
				next_generation = []
				idx_list = []
				survivors_by_elitism = []
				# choice individuals to survive by elitism
				for k in range(math.ceil(n*elitism_rate)):
					idx = int(sorted_actual_pop[k][0])
					idx_list.append(idx)
					survivors_by_elitism.append(population[idx])

				t_survivor = time.time() - t_survivor
				# print('sorted', sorted_actual_child)
				# print('next_generation after elitism\n', next_generation)

				# print(len(population))
				# 'cause the problem that remove idx i before j, such as i < j
				idx_list = sorted(idx_list, reverse=True) 
				for idx in idx_list:
					del population[idx] # delete individuals that survivor by elitism

				# TODO: roulette method
				# choice with equal probability get some individual in the population
				j = 0
				while j < math.floor(n*survival_rate) - len(idx_list):
					idx = random.randint(0,len(population)-1)
					next_generation.append(population[idx])
					del population[idx]
					j += 1

				population = [i for i in next_generation]

				# time_ls = time.time()
				# ******************************* LOCAL SEARCH STEP *******************************
				# print(len(survivors_by_elitism))
				for i in range(len(survivors_by_elitism)):
					survivors_by_elitism[i] = local_search(survivors_by_elitism[i],math.ceil(lcl_srch_rate),machines, processing_times, count_jobs, count_mach,  dict_chromossomes)
				
				# ******************************* MUTATION STEP *******************************
				population = population + survivors_by_elitism

				while len(population) < n:
					ind = build_random_chromossome(count_jobs, count_mach)
					population.append(ind)


				random.shuffle(population)

				# for ind in population:
				# 	if tuple(ind) not in dict_chromossomes:
				# 		value = evaluate(ind, machines, processing_times, count_jobs, count_mach)
				# 		# dict_chromossomes[tuple(ind)] = value
				# 		# dict_actual_pop[i] = value
				# 		dict_ind_pop[tuple(ind)] = value
				# 	else:
				# 		# dict_actual_pop[i] = dict_chromossomes[tuple(ind)]
				# 		dict_ind_pop[tuple(ind)] = dict_chromossomes[tuple(ind)]
				# 	# i += 1
				# print('len', len(dict_ind_pop.keys()))

				# sorted_population = sorted(dict_test.items(), key=operator.itemgetter(1))
				# print(sorted_population)
				# print('best solution of iteration', it, ':\n', sorted_population[0])
				# print('time survivor step:   ', t_survivor)
				# print('iteration', it, '\n', sorted_population)
				t_iteration = time.time() - t_iteration
				if it == 0:
					file_times.write('time build population: '+str(t_build_pop)+'\n')
					file_times.write('time crossover step:   '+ str(t_cross_pop)+'\n')
					file_times.write('time by iteration:     '+str(t_iteration)+'\n')


			time_exec = time.time() - time_exec

			print('time to execute', time_exec)
			file_times.write('time to execute        '+str(time_exec)+'\n')
			file_times.write('**********************************\n\n')
			# file_results.write('best in last population: ', population)
			sorted_history = sorted(dict_chromossomes.items(), key=operator.itemgetter(1))
			file_results.write('best in history: '+str((sorted_history[0][1]))+'\n')
			file_results.write('\n\n')
			print('best in history', sorted_history[0])
			# print('new generation', len(next_generation))
			if it == num_it-1:
				sum_total.append(sorted_history[0][1])
		mean = np.mean(sum_total)
		std = np.std(sum_total)
		file_results.write('\naverage:  '+str(mean)+'  std:  '+str(std))

		file_results.close()
		file_times.close()

			# if it == 0:
			# 	