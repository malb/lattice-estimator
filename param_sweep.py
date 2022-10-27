import argparse
from estimator import LWE
from estimator.lwe_parameters import LWEParameters
from estimator.nd import NoiseDistribution
from estimator.schemes import TFHE630
import pickle
import time
import math
import multiprocessing
import numpy as np

# To run: sage --python param_sweep.py 600:1201:20 6:20:1

keys = ['arora-gb', 'bkw', 'usvp', 'bdd', 'bdd_hybrid', 'bdd_mitm_hybrid', 'dual', 'dual_hybrid', 'dual_mitm_hybrid']

def performCalculation(input_params, results):
	# input_params looks like (n, xe)
	n, xe = input_params
	print(f'Running for {n=}, {xe=}')
	LWEParams = LWEParameters(
		n=n,
		q = 2 ** 32,
	    Xs=NoiseDistribution.UniformMod(2),
	    Xe=NoiseDistribution.DiscreteGaussian(stddev=2**xe),
	    tag=f'N{n}_Xe{xe}',
	)
	attacks = LWE.estimate(LWEParams)
	results[(n, xe)] = min([math.log(attacks[i].rop, 2) for i in keys if i in attacks])
	print(f'Finished running {n=}, {xe=}')

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description='Iterate over possible options.')
	parser.add_argument('nrange', type=str, help='range for N, like 500:600:1', default="500:600:1")
	parser.add_argument('xerange', type=str, help='range for Xe: log(sigma), like 10:20:1. Logarithmic powers', default="10:20:1")

	args = parser.parse_args()

	nrange = [int(i) for i in args.nrange.split(':')]
	xerange = [float(i) for i in args.xerange.split(':')]
	assert len(nrange) == 3, 'n range needs to be 3 numbers separated by a colon'
	assert len(xerange) == 3, 'Xe range needs to be 3 numbers separated by a colon'

	result_dict = multiprocessing.Manager().dict()
	work = []
	for n in range(*nrange):
		for xe in np.arange(*xerange):
			work.append((n, xe))

	pool = multiprocessing.Pool(processes=8)
	for i in range(len(work)):
		pool.apply_async(performCalculation, (work[i], result_dict))
	pool.close()
	pool.join()

	name = str(int(time.time()))
	name = f'{name}_{args.nrange}_{args.xerange}_result.pickle'
	pickle.dump(dict(result_dict), open(name, 'wb'))
