import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns

PI = 3.1415926535897

def main():
    K = int(args[2])
    cancer_type_list = load_types(args[1])
    cancer_type_num = []
    for x in cancer_type_list:
        File = open('data/data' + args[1] + '_' + x + '.txt', 'r')
        cancer_type_num.append(int(File.readline().split()[0]))
        File.close()
    doc_num = sum(cancer_type_num)

    doc_arrange, alpha_list, Average_arrange = load_result(doc_num, \
        cancer_type_list, cancer_type_num)

    for i,x in enumerate(cancer_type_list):
        left = range(1, K+1)
        height = alpha_list[i]
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.bar(left, height, width=1, align="center")
        ax1.set_xlim(0, K+1)
        ax1.set_ylim(0, 5.0)
        ax1.set_xticks(left)
        ax1.tick_params(labelsize=10)
        ax1.set_xlabel("Signature", fontsize=10)
        ax1.set_ylabel("Alpha", fontsize=10)
        Title = x + " alpha"
        ax1.set_title(Title, fontsize=10)
        fig.tight_layout()
        name = "result/data" + args[1] + "/figure/" + args[2] + "_arrangement/" +\
                x + "_alpha.png"
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

def load_result(doc_num, cancer_type_list, cancer_type_num):
    K = int(args[2])
    doc_arrange = []
    alpha_list = []
    if(K <= 9):
        topic = '0' + args[2]
    else:
        topic = args[2]
    File = open('result/data' + args[1] + '/result_k' + topic + '.txt', 'r')
    File.readline(); File.readline();
    for i in range(K):
        File.readline()
    for i in range(doc_num):
        doc_arrange.append([])
        temp_list = File.readline().split()
        for j in range(K):
            doc_arrange[i].append(float(temp_list[j]))
    for i in range(len(cancer_type_list)):
        alpha_list.append([])
        temp_list = File.readline().split()
        for j in range(K):
            alpha_list[i].append(float(temp_list[j]))
    File.close()
    Average_arrange = []
    index = 0
    for i in range(len(cancer_type_list)):
        Average_arrange.append([0 for k in range(K)])
        for j in range(cancer_type_num[i]):
            for k in range(K):
                Average_arrange[i][k] += doc_arrange[index][k]
            index += 1
        for k in range(K):
            Average_arrange[i][k] /= cancer_type_num[i]
    return doc_arrange, alpha_list, Average_arrange

if __name__ == '__main__':
    args = sys.argv ## args[1] : data_type, [2] : num_topic
    main()
