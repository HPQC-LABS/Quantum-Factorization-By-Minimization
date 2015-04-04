import re
import sys
from multiprocessing import Process, Manager
import time

##########################Change these if using on new computer########################

INPUT_FILE_DIRECTORY = "/home/tanburn/Quantum-Factorization-By-Minimization/semiprimes_hamming/"
OUTPUT_DIRECTORY = '/home/tanburn/Quantum-Factorization-By-Minimization/sage_timings/'

#######################################################################################

manager=Manager()

def extract(str_):
    match_ = re.search(r'(\d+\.\d+) ?user', str_)
    if match_:
        return match_.groups()[0]
    else:
        return str_

def run_qsieve(semiprime,ans):
    (factors,cpu_time)=qsieve(semiprime,time=True)
    ans['cpu_time']=cpu_time

def take_time(semiprime, previous_time):
    ans = manager.dict()
    p = Process(target=run_qsieve,args=(semiprime,ans,))
    p.start()
    if previous_time is None or previous_time is 'N/A':
        p.join()
    else:
        p.join(float(previous_time)*10)

    if p.is_alive():
        print "killing qsieve for " + str(semiprime)
        # Terminate
        p.terminate()
        p.join()
    if 'cpu_time' in ans:
        return extract(ans['cpu_time'])
    else:
        return 'N/A'

def import_data(file_name):
    input_file=open(file_name,'r')
    data= [eval(x.replace('\n','')) for x in input_file.readlines() if x.replace('\n','') is not '']
    input_file.close()
    return data

def run_factorizations(N):
    try:
        input_file_name = INPUT_FILE_DIRECTORY+"hamming_primes_"+str(N)+"x"+str(N)+".txt"
        output_file_name = OUTPUT_DIRECTORY+"timings_"+str(N)+"x"+str(N)+".txt"

        semi_primes = import_data(input_file_name)

        output_file=open(output_file_name,'w',0)

        sys.stdout = output_file

        previous_time = None
        
        for (semiprime,hamming_distance,p,q) in semi_primes:
            previous_time = take_time(semiprime,previous_time)

            print '('+previous_time+','+str(semiprime)+','+str(hamming_distance)+','+str(p)+','+str(q)+')'+'\n'
        
        sys.stdout = sys.__stdout__
        output_file.flush()
        output_file.close()
    except KeyboardInterrupt:
        sys.stdout = sys.__stdout__
        output_file.flush()
        output_file.close()


#Run script:


for N in range(70,530,10):
    run_factorizations(N)
