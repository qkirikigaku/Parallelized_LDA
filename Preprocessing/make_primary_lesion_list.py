import os

def main():
    files = os.listdir('../data/')
    output_file1 = '../data/PL_data1_list.txt'
    output1 = open(output_file1, 'w')
    data1_list = []
    output_file2 = '../data/PL_data2_list.txt'
    output2 = open(output_file2, 'w')
    data2_list = []
    for x in files:
        if (x.startswith('data1_')):
            end_index = x.index('.txt')
            data1_list.append(x[6:end_index])
        if (x.startswith('data2_')):
            end_index = x.index('.txt')
            data2_list.append(x[6:end_index])
    output1.write(str(len(data1_list)) + '\n')
    output2.write(str(len(data2_list)) + '\n')
    for x in data1_list:
        output1.write(x + '\n')
    for x in data2_list:
        output2.write(x + '\n')
    output1.close()
    output2.close()

if __name__ == '__main__':
    main()
