from numpy.random import beta

def hamming_distance(p,q):
	dist=0
	for a, b in zip(list(Integer(p).binary()),list(Integer(q).binary())):
		dist+=abs(int(a)-int(b))
	return dist

def generate_prime(lower_bound,upper_bound):
	try:
		random_number = lower_bound+floor((upper_bound-lower_bound)*beta(0.03,0.03))

		prime = next_prime(random_number)

		if prime < upper_bound:
			return prime
		else:
			return previous_prime(random_number)

	except ValueError:
		print "Error"
		return generate_prime(lower_bound,upper_bound)

def get_prime(lower_bound,upper_bound,existing_primes):
	p=generate_prime(lower_bound,upper_bound)
	while p in existing_primes:
		p=generate_prime(lower_bound,upper_bound)
	return p

def check_completition_criteria(bins):
	if False in map(lambda x: len(x)>=5, bins):
		return True
	else:
		return False

def find_pairs(primes_list,bins):
	newest_prime=primes_list[-1]
	for prime in primes_list[:-1]:
		bins[ceil(float(hamming_distance(newest_prime,prime))/5)-1].append((prime,newest_prime))


def generate_factors(n):
	number_of_bins=int(n/5)
	lower_bound=next_prime(2**(n-1))
	upper_bound=previous_prime(2**n)

	n_prime = next_prime(lower_bound)
	p_prime = previous_prime(upper_bound)

	primes_list=[lower_bound,upper_bound]
	bins=[[] for i in range(number_of_bins)]

	while check_completition_criteria(bins):
		find_pairs(primes_list,bins)

		new_prime = get_prime(lower_bound,upper_bound,primes_list)

		if new_prime == n_prime:
			lower_bound = new_prime
			n_prime = next_prime(lower_bound)
		elif new_prime == p_prime:
			upper_bound = new_prime
			p_prime = previous_prime(upper_bound)

		primes_list.append(new_prime)
		
	return bins

def pick_ten(factor_list):
	return [x[:10] for x in factor_list]


for n in range(20,120,10):
	factors=pick_ten(generate_factors(n))

	f=open("/Users/emileokada/Desktop/semi_primes_"+str(n)+".txt",'w')
	for x, index in zip(factors,range(len(factors))):
		f.write('Hamming distance in range: ['+str(index*5+1)+','+str((index+1)*5)+']\n')
		for p, q in x:
			f.write('('+str(p*q)+','+str(hamming_distance(p,q))+','+str(p)+','+str(q)+')'+'\n')
		f.write('\n')
	f.close()