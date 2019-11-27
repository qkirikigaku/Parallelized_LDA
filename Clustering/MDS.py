import sys
import numpy as np
from sklearn.manifold import MDS
from scipy.stats import entropy
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    cancer_type_list = load_types(args[1])
    cancer_type_num = []
    for x in cancer_type_list:
        File = open("data/data" + args[1] + "_" + x + ".txt", "r")
        cancer_type_num.append(int(File.readline().split()[0]))
        File.close()
    doc_num = sum(cancer_type_num)
    doc_arrange, alpha_list, Average_arrange = load_result(doc_num, cancer_type_list,
            cancer_type_num)
    markers = ["o", "^", "s"]
    c_palette = sns.color_palette("hls", len(cancer_type_list))
    clf = MDS(n_components=2, dissimilarity="precomputed", n_init=1, max_iter=100)
    JS_dis = np.zeros([len(doc_arrange), len(doc_arrange)])
    for i in range(len(doc_arrange)):
        for j in range(len(doc_arrange)):
            m = list()
            for k in range(K):
                m.append((doc_arrange[i][k] + doc_arrange[j][k])/2)
            JS_dis[i,j] = (entropy(doc_arrange[i], m)+ entropy(doc_arrange[j], m))/2
    transformed = clf.fit_transform(JS_dis)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    index = 0
    for i,x in enumerate(cancer_type_list):
        temp_first_components = []
        temp_second_components = []
        for j in range(cancer_type_num[i]):
            temp_first_components.append(transformed[index][0])
            temp_second_components.append(transformed[index][1])
            index +=  1
        ax.scatter(temp_first_components, temp_second_components, label=x, \
                c=c_palette[i], marker= markers[i%3], s=5)
    ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1), borderaxespad=0)
    fig.subplots_adjust(right=0.65)
    plt.title("MDS of activities (JS divergence)")
    name = "result/data" + args[1] + "/figure/" + args[2]\
            + "_arrangement/mds_JS.png"
    plt.savefig(name, dpi=300)
    plt.close(1)

def load_types(data_type):
    File = open('data/PL_data' + data_type + '_list.txt', 'r')
    File.readline()
    cancer_type_list = []
    for line in File.readlines():
        cancer_type_list.append(line[:-1])
    File.close()
    return cancer_type_list

def load_result(doc_num, cancer_type_list, cancer_type_num):
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

if __name__ == "__main__":
    args = sys.argv ## args[1]: data_type, [2] : num_topic
    K = int(args[2])
    main()
