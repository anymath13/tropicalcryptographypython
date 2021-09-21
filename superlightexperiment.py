from numpy import linalg as LA
from time import process_time
import numpy as np
import secrets
from time import process_time
import tropicalalgebra as tp
from superlightattack import *
from superlightdiscretelog import *
from threecomponentsmatrices import *




def test_attack(reps, max_terms, max_period,
                matrix_order, min,max):
    term_indexes = []
    orders = []
    times = []
    errors = []

    for i in range(reps):
        #M= [[-317, 937, -963, -705, 260, -82, 143, -74, 265, -135], [589, 837, 522, 589, -817, 135, 146, -783, -79, 318], [398, 455, -379, 399, 280, 82, -800, 802, 824, -95], [-709, 743, 971, -454, -332, 406, 938, 331, -178, -965], [207, -126, 844, -275, 412, 820, 763, 843, 19, 885], [-323, -220, -199, 903, -891, -272, -428, -949, -586, 894], [-490, -914, -520, 422, 829, 262, 981, -1, -115, -29], [383, 641, -718, -940, 370, -657, 144, -232, -466, 27], [-20, 445, 223, 469, 686, 614, 653, 978, 7, -789], [-369, -850, 957, -60, 628, 904, 641, 864, -8, 433]]
        #H= [[-166, -542, 599, 559, -520, 7, 355, 589, 906, 399], [976, -414, 508, -36, 847, -24, 75, 308, -405, -535], [-261, 653, 100, -932, 849, 36, 331, 998, 304, -417], [433, 303, 517, 38, -547, 192, 315, -167, -852, -967], [269, -619, -155, 996, -150, -361, 861, -897, -248, -168], [-393, 700, 905, -208, 477, -700, -736, -945, 605, 476], [151, -166, 114, -466, 464, -505, 521, 620, 159, -665], [-51, 725, -307, -470, -858, 271, 199, 580, 527, 926], [-904, 594, 243, -328, 562, 938, 377, -3, -251, -545], [-871, 707, 463, 228, 73, 218, -37, -602, -403, -471]]
        #ma= 104243141396510889934473744799720789817242410794263334808119
        #nb= 1241786337007531973597923782253776357075970508879224881387292
        M=tp.generate_random_matrix(matrix_order, matrix_order,min,max)
        print('M=',M)
        H=tp.generate_random_matrix(matrix_order, matrix_order,min,max)
        #H=generatespecialmatrix(matrix_order)
        #H= [[-438.0, -26.0, 9.0, -916.0, -569.0, -588.0, -188.0, -738.0, -676.0, -919.0], [-383.0, -61.0, -18.0, -548.0, -80.0, -817.0, -71.0, -282.0, -380.0, -782.0], [-449.0, -18.0, 17.0, -513.0, -682.0, -579.0, -240.0, -176.0, -589.0, -565.0], [-251.0, -868.0, -914.0, -685.0, -499.0, -687.0, 33.0, 54.0, 66.0, -799.0], [-383.0, -783.0, -423.0, -477.0, -428.0, -461.0, -26.0, -527.0, 7.0, -107.0], [3.0, -636.0, 21.0, -171.0, -343.0, 175.0, -402.0, 102.0, 12.0, -784.0], [-429.0, -1014.0, -567.0, -720.0, -26.0, 81.0, -13.0, -287.0, 20.0, -607.0], [-915.0, -495.0, -254.0, -327.0, -217.0, 102.0, 8.0, -764.0, -257.0, -190.0], [-379.0, -307.0, -481.0, -162.0, 7.0, -81.0, -566.0, -97.0, -624.0, -481.0], [-670.0, -38.0, -501.0, -270.0, -678.0, -210.0, -29.0, -587.0, -862.0, 19.0]]
        print('H=',H)
        ma=tp.generate_exponent(max_period)
        print('ma=',ma)
        nb=tp.generate_exponent(max_period)
        print('nb=',nb)
        A1=tp.adjointpowerfast(M,H,ma)
        A=A1[0]
        B1=tp.adjointpowerfast(M,H,nb)
        B=B1[0] 
        start_time = process_time()
        K_attack = superlightattack(M,H,A,B,ma,nb)
        end_time = process_time()
        times.append(end_time-start_time)
        #print('K_attack[0]-ma',K_attack[0]-ma)
        if (K_attack[0]-ma==0)&(K_attack[1]-nb==0):
            
            term_indexes.append(K_attack[0])
            #print(term_indexes)
            orders.append(K_attack[1])
            #print(orders)
        else:
            errors.append((M,H,A,B,ma,nb))
            print("Attack could not find exponent")
        if (i + 1) % 100 == 0:
            print('Iteration: ', i + 1)
            print('Percentage broken: ', (len(term_indexes) * 100) / (i + 1), "%")
            #print('Max terms searched: ', max(term_indexes))
            #print("Max cycle length: ", max(orders))
   
    print("Matrix order:", matrix_order)
    print("Repetitions:", reps)
    print("Number broken:", len(term_indexes))
    #print("Max terms searched:", max(term_indexes))
    #print("Max cycle length:", max(orders))
    #print("Cycle occurrences:")
    #print('times',(times))
    
    average_times=sum(times)/len(times)
    print('average_times',average_times)
    #for i in set(orders):
        #print(i, ":", orders.count(i))

    return term_indexes, orders, times, errors
  
test_attack(1000,10000,200,7,-1000,1000)