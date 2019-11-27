import sys
import numpy as np

def main():
    args = sys.argv
    ## [1]:number of topic
    K = int(args[1])
    if(K <= 9):
        topic = '0' + args[1]
    else:
        topic = args[1]

    input_file = 'result/data2/result_k' + topic + '.txt'

    n2_data = load_data(input_file, K)
    n1_data = to_n1(n2_data, K)

    output_file = 'result/data2/figure/' + args[1] + '_signature/k' + topic + '.txt'
    write_data(output_file, n1_data, K)

def load_data(input_file, K):
    result = open(input_file,'r')
    lines = result.readlines()
    count = -2
    for line in lines:
        if(count ==  0):
            probability = line.split(' ')
            data = np.zeros([K,len(probability) - 1])
            for i in range(len(probability) - 1):
                data[count,i] = float(probability[i])
        elif((count > 0) and (count < K)):
            probability = line.split(' ')
            for i in range(len(probability) - 1):
                data[count,i] = float(probability[i])
        count += 1
    result.close()
    return data

def to_n1(n2_data, K):
    V = 1536
    index = [0 for i in range(V)]
    for v in range(V):
        forward = (v % 384) // 96
        substitution = ((v % 384) % 96) //16
        backward = (((v % 384) % 96) % 16) //4
        index[v] = 24 * forward + 4 * substitution + backward
    data = np.zeros([K,96])
    for k in range(K):
        for i in range(V):
            data[k,index[i]] += n2_data[k,i]
    return data

def write_data(output_file, n1_data, K):
    output = open(output_file, 'w')
    output.write('0\n')
    output.write('0\n')
    for k in range(K):
        for i in range(96):
            if(i != 95):
                string = str(n1_data[k,i]) + ' '
            else:
                string = str(n1_data[k,i]) + '\n'
            output.write(string)
    output.close()

if __name__ == '__main__':
    main()
