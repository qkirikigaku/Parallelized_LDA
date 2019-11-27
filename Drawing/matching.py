import sys
import numpy as np
from scipy import stats
from scipy.spatial.distance import cosine


def swap(i,list):
    temp = list[i]
    list[i] = list[i+1]
    list[i+1] = temp

def main():
    args = sys.argv
    """args[1]:data_type [2]:number_of_topic"""

    K = int(args[2])
    if(K <= 9):
        topic = '0' + args[2]
    else:
        topic = args[2]

    input_file = 'result/data' + args[1] + '/result_k' + topic + '.txt'

    ex = open(input_file)
    co = open('data/signature_probability.txt')

    data_ex = ex.readlines()
    count_ex = 0
    p_ex = np.zeros([K,96])
    for line in data_ex:
        if((count_ex > 1) and (count_ex < K + 2)):
            words = line.split()
            for signature in range(96):
                p_ex[count_ex-2,signature] = float(words[signature])
        count_ex += 1

    data_co = co.readlines()
    count_co = 0
    p_co = np.zeros([30,96])
    for line in data_co:
        if(count_co != 0):
            words = line.split()
            for signature in range(33):
                if(signature >= 3):
                    p_co[signature-3,count_co-1] = float(words[signature])
        count_co += 1

    ex.close()
    co.close()

    JS_ex = np.zeros([K,30])
    JS_co = np.zeros([30,K])
    for i in range(K):
        for j in range(30):
            #qx_ex = (p_ex[i]+p_co[j])/2
            #JS_ex[i,j] = 0.5*(stats.entropy(p_ex[i],qx_ex,2) + stats.entropy(p_co[j],qx_ex,2))
            JS_ex[i,j] = cosine(p_ex[i],p_co[j])

    for i in range(30):
        for j in range(K):
            #qx_co = (p_co[i]+p_ex[j])/2
            #JS_co[i,j] = 0.5*(stats.entropy(p_co[i],qx_co,2) + stats.entropy(p_ex[j],qx_co,2))
            JS_co[i,j] = cosine(p_co[i], p_ex[j])

    print(JS_ex)

    Pre_ex = np.argsort(JS_ex)
    Pre_co = np.argsort(JS_co)

    min_JS = np.zeros([K])
    min_sig = np.zeros([K])

    match = '{'
    for i in range(K):
        min_j = 0
        for j in range(1,30):
            if(JS_ex[i,j] < JS_ex[i,min_j]):
                min_j = j
        min_JS[i] = JS_ex[i,min_j]
        min_sig[i] = min_j
        match += str(i) + ': ' + str(min_j)
        if(i != K-1):
            match += ', '
    match += '}'

    output_file = "result/data" + args[1] + '/figure/' + args[2] + '_signature/matching.txt' 

    output = open(output_file,'w')
    output.write(match)
    output.close()

    output_right_sig = open('result/data' + args[1] + '/figure/' + args[2] + '_signature/match_right_sig.txt' , 'w')

    right_sig_zero_two = 0
    right_sig_zero_one_five = 0
    right_sig_zero_one = 0
    right_sig_zero_zero_five = 0
    for i in range(K):
        if (min_JS[i] < 0.05):
            right_sig_zero_two += 1
            right_sig_zero_one_five += 1
            right_sig_zero_one += 1
            right_sig_zero_zero_five += 1
        elif(min_JS[i] < 0.1):
            right_sig_zero_two += 1
            right_sig_zero_one_five+= 1
            right_sig_zero_one += 1
        elif(min_JS[i] < 0.15):
            right_sig_zero_two += 1
            right_sig_zero_one_five+= 1
        elif(min_JS[i] < 0.2):
            right_sig_zero_two += 1
    output_right_sig.write('0.2以下 : ' + str(right_sig_zero_two) + '\n')
    output_right_sig.write('0.15以下 : ' + str(right_sig_zero_one_five) + '\n')
    output_right_sig.write('0.1以下 : ' + str(right_sig_zero_one) + '\n')
    output_right_sig.write('0.05以下 : ' + str(right_sig_zero_zero_five) + '\n')

    for i in range(K):
        print('Predicted signature ' + str(i) + ' : ' + str(min_JS[i]) + '\n')

    output_sig = open('result/data' + args[1] + '/figure/' + args[2] + '_signature/match_sigs.txt', 'w')
    output_sig.write('0.20以下 : ')
    out_zero_two_zero = list()
    for i in range(K):
        if(min_JS[i] < 0.2):
            out_zero_two_zero.append(min_sig[i])
    out_zero_two_zero.sort()
    for i in range(len(out_zero_two_zero)):
        output_sig.write(str(int(out_zero_two_zero[i] + 1)) + ', ')

    out_zero_two_five = list()
    output_sig.write('\n' + '\n' + '0.25以下 : ')
    for i in range(K):
        if(min_JS[i] < 0.25):
            out_zero_two_five.append(min_sig[i])
    out_zero_two_five.sort()
    for i in range(len(out_zero_two_five)):
        output_sig.write(str(int(out_zero_two_five[i] + 1)) + ', ')

    out_zero_three_zero = list()
    output_sig.write('\n' + '\n' + '0.30以下 : ')
    for i in range(K):
        if(min_JS[i] < 0.30):
            out_zero_three_zero.append(min_sig[i])
    out_zero_three_zero.sort()
    for i in range(len(out_zero_three_zero)):
        output_sig.write(str(int(out_zero_three_zero[i] + 1)) + ', ')

    output_not_sig = open('result/data' + args[1] + '/figure/' + args[2] + '_signature/not_match_predicted.txt', 'w')
    output_not_sig.write('0.50以上 :')
    for i in range(K):
        if(min_JS[i] >= 0.5):
            output_not_sig.write(str(int(i + 1)) + ', ')

    min_cos = open('result/data' + args[1] + '/figure/' + args[2] + '_signature/minimum_cos.txt', 'w')
    for i in range(K):
        min_cos.write(str(min_JS[i]) + '\n')

if __name__ == '__main__':
    main()
