This version doesn't make any assumptions but instead keeps track of how many qubits are possible and how many are solved, then reports this and writes the final equations.

.\PrimeFactorizationEquationGenerator.py  4 4 3454 filename.txt

4 4 - the length in binary of p and q
3454 - the number to be factored
filename.txt - the file to append the output to

=======================================================================
The format has to be like this:

index_i index_j index_k ...    coefficient_i_j_k_....
…

--------------------

For example to minimize the function of binary variables s1, s2, s3:

min_s1,s2,s3 (5+40*s1 + 30*s1*s2 + 100*s1*s2*s3 - 15*s2*s3 - 20*s3)

the input file is:

5
1 40
1 2 30
1 2 3 100
2 3 -15
3 -20
