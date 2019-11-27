import sys
import re
import pandas as pd
import numpy as np
import linecache

from get_M1 import select_data
from get_M1 import cut_out_threshold

def main():
    args = sys.argv
    """args[1]:threshold number of word"""

    dictionary = load_dictionary()
    temp_dictionary = list()

    K = len(dictionary)
    threshold = int(args[1])
    cancer_type_list = ['adrenal_gland', 'biliary_tract', 'bone', 'breast', 'central_nervous_system', 'cervix', 'endometrium',\
            'haematopoietic_and_lymphoid_tissue', 'kidney', 'large_intestine', 'liver', 'lung', 'oesophagus', 'ovary', 'pancreas', 'parathyroid',\
            'prostate', 'skin', 'small_intestine', 'soft_tissue', 'stomach', 'thyroid', 'upper_aerodigestive_tract',\
            'urinary_tract', 'vulva']

    Mutation_file = '~/Data/Mutation_Signature/CosmicMutantExport.tsv'
    pre_data = pd.read_csv(Mutation_file, delimiter = '\t')
    pre_data = select_data(pre_data)
    pre_data.reset_index(drop=True, inplace=True)

    pre_data = cut_out_threshold(pre_data, threshold)
    pre_data.reset_index(drop=True, inplace=True)

    for cancer_type in cancer_type_list:
        data = select_cancer_type(pre_data, threshold, cancer_type)
        data.reset_index(drop=True, inplace=True)
        name_list = data['Sample name']
        last_document = 'xxxxxxxxxx'
        number_of_document = 0
        for name in name_list:
            if(name != str(last_document)):
                last_document = name
                number_of_document += 1
        if (number_of_document >= 10):
            output_to_file(data, cancer_type, number_of_document, K, dictionary, temp_dictionary)

def output_to_file(pre_data, cancer_type, number_of_document, K, dictionary, temp_dictionary):
    output_file = '../data/data2_' + cancer_type + '.txt'
    output = open(output_file, 'w')

    name_list = pre_data['Sample name']
    last_document = name_list[0]
    data_mat = np.zeros([number_of_document, K])
    index = 0
    error = 0
    document = 0
    for name in name_list:
        if(name != str(last_document)):
            last_document = name
            document += 1
        selected = calc_word_M2(pre_data['Mutation CDS'][index],pre_data['Mutation genome position'][index],pre_data['Mutation strand'][index])
        if(selected == -1):
            error += 1
        else:
            temp_dictionary = write_temp_dictionary(temp_dictionary, dictionary, selected)
            data_mat[document, selected] += 1
        index += 1
        print(index)
    
    """
    #check
    drop_list = list()
    for i in range(number_of_document):
        sum_words = 0
        for j in range(K):
            sum_words += data_mat[i,j]
        if(sum_words == 0):
            drop_list.append(i)
    number_of_document -= len(drop_list)
    data_mat = np.delete(data_mat,drop_list,0)

    print(str(len(drop_list)))

    re_data = rewrite_data_mat(data_mat,dictionary,temp_dictionary,number_of_document)
    """

    re_data = data_mat.copy()
    number_of_vocab = 1536
    output.write(str(number_of_document) + ' ')
    output.write(str(number_of_vocab) + '\n')

    for i in range(number_of_document):
        for j in range(number_of_vocab):
            if(j != number_of_vocab - 1):
                output.write(str(int(re_data[i][j])) + ' ')
            else:
                output.write(str(int(re_data[i][j])))
        output.write('\n')

    """
    temp_dict_file = '../dictionary/data2_' + cancer_type +'.txt'
    temp_dict = open(temp_dict_file, 'w')
    for i in range(number_of_vocab):
        temp_dict.write(temp_dictionary[i] + '\n')
    """

def select_cancer_type(pre_data, threshold, cancer_type):
    if(cancer_type != 'all'):
        pre_data = pre_data.loc[pre_data['Primary site'] == cancer_type]
    print(pre_data)
    pre_data.sort_values(by='Sample name', inplace=True)
    pre_data.reset_index(drop=True, inplace=True)
    name_list = pre_data['Sample name']
    last_document = name_list[0]
    drop_list = list()
    sum_of_words = 0
    temp_index = 0
    index_list = list()
    for name in name_list:
        if(name != last_document):
            if(sum_of_words < threshold):
                drop_list.extend(index_list)
            sum_of_words = 0
            last_document = name
            index_list = list()
        sum_of_words += 1
        index_list.append(temp_index)
        temp_index += 1
    if(sum_of_words < threshold):
        drop_list.extend(index_list)
    pre_data.drop(drop_list, inplace=True)
    return pre_data

def calc_word_M2(mutation, position, strand):
    before = mutation[len(mutation)-3]
    after = mutation[len(mutation)-1]
    position_list = re.split(r'[:-]',position)
    if(len(position_list) != 3):
        print('position error')
        return -1
    if(int(position_list[0]) == 23):
        chromosome = 'X'
    elif(int(position_list[0]) == 24):
        chromosome = 'Y'
    elif(int(position_list[0]) == 25):
        chromosome = 'M'
    else:
        chromosome = int(position_list[0])
    start = int(position_list[1])
    num = int(position_list[2]) - int(position_list[1]) + 1
    GRCh_file = '/Users/tarom/Data/Mutation_Signature/chr' + str(chromosome) + '.fa'
    quotient = start // 50
    surplus = start % 50

    if(surplus != 0):
        target_index = int(surplus) - 1
    else:
        quotient -= 1
        target_index = 49
    target_line = linecache.getline(GRCh_file, int(quotient)+1)
    
    if(((target_line[target_index] != before) and (strand == '+')) or ((target_line[target_index] != swap(before))and(strand == '-'))):
        print('error: ' + mutation)
        print('target: ' + target_line[target_index])
        print('strand: ' + strand)
        strand = swap(strand)
        if(((target_line[target_index] != before) and (strand == '+')) or ((target_line[target_index] != swap(before))and(strand == '-'))):
            print('still error')
            return -1

    if((target_index >= 2) and (target_index <= 47)):
        pattern = 1
    elif(target_index == 0):
        pattern = 2
    elif(target_index == 1):
        pattern = 3
    elif(target_index == 48):
        pattern = 4
    else:
        pattern = 5

    if(pattern == 1):
        forward = target_line[target_index - 1]
        for_forward = target_line[target_index - 2]
        backward = target_line[target_index + 1]
        back_backward = target_line[target_index + 2]
    elif(pattern == 2):
        pre_line = linecache.getline(GRCh_file, int(quotient))
        forward = pre_line[49]
        for_forward = pre_line[48]
        backward = target_line[target_index + 1]
        back_backward = target_line[target_index + 2]
    elif(pattern == 3):
        pre_line = linecache.getline(GRCh_file, int(quotient))
        for_forward = pre_line[49]
        forward = target_line[target_index - 1]
        backward = target_line[target_index + 1]
        back_backward = target_line[target_index + 2]
    elif(pattern == 4):
        post_line = linecache.getline(GRCh_file, int(quotient)+2)
        back_backward = post_line[0]
        forward = target_line[target_index - 1]
        for_forward = target_line[target_index - 2]
        backward = target_line[target_index + 1]
    if(pattern == 5):
        post_line = linecache.getline(GRCh_file, int(quotient)+2)
        backward = post_line[0]
        back_backward =post_line[1]
        forward = target_line[target_index - 1]
        for_forward = target_line[target_index - 2]

    if(((strand == '+') and (before in ['A', 'G'])) or ((strand == '-') and (before in ['C', 'T']))):
        buf_f = swap(forward)
        buf_ff = swap(for_forward)
        forward = swap(backward)
        for_forward = swap(back_backward)
        backward = buf_f
        back_backward = buf_ff
    if(before in ['A', 'G']):
        before = swap(before)
        after = swap(after)

    if(for_forward == 'A'):
        first = 0
    elif(for_forward == 'C'):
        first = 1
    elif(for_forward == 'G'):
        first = 2
    else:
        first = 3

    if(forward == 'A'):
        second = 0
    elif(forward == 'C'):
        second = 1
    elif(forward == 'G'):
        second = 2
    else:
        second = 3

    if(before == 'C'):
        if(after == 'A'):
            third = 0
        elif(after == 'G'):
            third = 1
        else:
            third = 2
    elif(before == 'T'):
        if(after == 'A'):
            third = 3
        elif(after == 'C'):
            third = 4
        else:
            third = 5
    elif(before == 'G'):
        if(after == 'T'):
            third = 0
        elif(after == 'C'):
            third = 1
        else:
            third = 2
    else:
        if(after == 'T'):
            third = 3
        elif(after == 'G'):
            third = 4
        else:
            third = 5

    if(back_backward == 'A'):
        fourth = 0
    elif(back_backward == 'C'):
        fourth = 1
    elif(back_backward == 'G'):
        fourth = 2
    else:
        fourth = 3

    if(backward == 'A'):
        fifth = 0
    elif(backward == 'C'):
        fifth = 1
    elif(backward == 'G'):
        fifth = 2
    else:
        fifth = 3
    answer = 384*first + 96*second + 16*third + 4*fourth + fifth
    return(answer)


def swap(base):
    if(base == 'A'):
        return('T')
    elif(base == 'C'):
        return('G')
    elif(base == 'G'):
        return('C')
    elif(base == 'T'):
        return('A')
    elif(base == '+'):
        return('-')
    elif(base == '-'):
        return('+')
    else:
        return(base)

def load_dictionary():
    dictionary_file = '../dictionary/dictionary.txt'
    file = open(dictionary_file, 'r')
    lines = file.readlines()
    dictionary = list()
    for line in lines:
        text = line[:-1]
        dictionary.append(text)
    return dictionary

def write_temp_dictionary(temp_dictionary, dictionary, pre_index):
    flag = 0
    for i in temp_dictionary:
        if(i == dictionary[pre_index]):
            flag = 1
            break
    if(flag == 0):
        temp_dictionary.append(dictionary[pre_index])
    return temp_dictionary

def rewrite_data_mat(data_mat, dictionary, temp_dictionary, number_of_document):
    number_of_vocab = len(temp_dictionary)
    re_data = [[0 for i in range(number_of_vocab)] for j in range(number_of_document)]
    for i in range(number_of_document):
        for j in range(len(dictionary)):
            if(data_mat[i,j] != 0):
                index = temp_dictionary.index(dictionary[j])
                re_data[i][index] = data_mat[i,j]
    return re_data

if __name__ == '__main__':
    main()
