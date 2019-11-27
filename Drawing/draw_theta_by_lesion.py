import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

def main():
    args = sys.argv ## args[1] : data_type, [2] : num_topic
    
    K = int(args[2])
    cancer_type_list = load_types(args[1])
    cancer_type_num = []
    for x in cancer_type_list:
        File = open('data/data' + args[1] + '_' + x + '.txt', 'r')
        cancer_type_num.append(int(File.readline().split()[0]))
        File.close()
    doc_num = sum(cancer_type_num)

    if(K >= 10):
        topic = str(K)
    else:
        topic = '0' + str(K)
    File = open('result/data' + args[1] + '/result_k' + topic + '.txt', 'r')
    lines = File.readlines()[(3 + K):]
    index = 0
    heights = np.zeros([len(cancer_type_list), K])
    for i,x in enumerate(cancer_type_list):
        for j in range(cancer_type_num[i]):
            line = lines[index].split()
            index += 1
            for k in range(K):
                heights[i,k] += float(line[k])
        for k in range(K):
            heights[i,k] /= cancer_type_num[i]

    for i,x in enumerate(cancer_type_list):
        left = range(1, K+1)
        temp_height = {}
        for k in range(K):
            temp_height.update({k:heights[i,k]})
        height = list()
        labels = list()
        for key, v in sorted(temp_height.items(), key=lambda x: -x[1]):
            height.append(v)
            labels.append('Signature ' + str(key+1))
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.bar(left, height, width=1, align='center')
        ax1.set_ylim(0, max(height))
        ax1.set_xlim(0, K+1)
        ax1.set_xticks(left)
        ax1.set_xticklabels(labels)
        for tick in ax1.get_xticklabels():
            tick.set_rotation(90)
        ax1.tick_params(labelsize=10)
        ax1.set_xlabel('Signature', fontsize=10)
        ax1.set_ylabel('Average appearance', fontsize=10)
        Title = x
        ax1.set_title(Title, fontsize=10)
        fig.tight_layout()
        name = 'result/data' + args[1] + '/figure/' + args[2] + '_lesion/' + x + '_appearance.png'
        fig.savefig(name, dpi=200)
        plt.close(1)

def load_types(data_type):
    File = open('data/PL_data' + data_type + '.txt', 'r')
    File.readline()
    cancer_type_list = []
    for line in File.readlines():
        cancer_type_list.append(line[:-1])
    File.close()
    return cancer_type_list

if __name__ == '__main__':
    main()
