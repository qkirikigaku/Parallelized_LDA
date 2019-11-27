import numpy as np
import sys
import os.path

def main():
    args =  sys.argv
    # args[1]:data_type [2]:D [3]:K
    
    D = int(args[2]) #D : max_iteration
    K = int(args[3]) #K : max_num_topic
    VLB_mat = np.zeros([D,K-1])
    for d in range(1,D+1):
        for k in range(2,K+1):
            if(k <= 9):
                num_topic = '0' + str(k)
            else:
                num_topic = str(k)
            FILE = 'result/data' + args[1] + '_' + str(d) + '/result_k' + num_topic + '.txt'
            if os.path.exists(FILE):
                lines = open(FILE, 'r').readlines()
                VLB = lines[0]
            else:
                VLB = '-1000000000000000.0\n'
            VLB_mat[d-1,k-2] = float(VLB[:-1])
    index_list = VLB_mat.argmax(0)
    print('index_list : ' + '\n')
    print(index_list)
    index = VLB_mat.argmax()
    print('index : ' + str(index))
    max_d = (index // int(K-1)) + 1
    max_k = (index % int(K-1)) + 2
    print('max iteration : ' + str(max_d) + '\n')
    print('max topic num : ' + str(max_k) + '\n')
    output = open('buf/' + args[1] + '.txt', 'w')
    output.write(str(max_d) + '\n')
    output.write(str(max_k) + '\n')
    for d in range(K-1):
        output.write(str(index_list[d] + 1) + '\n')

if(__name__ == '__main__'):
    main()
