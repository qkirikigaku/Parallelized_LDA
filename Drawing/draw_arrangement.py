import sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
from sklearn.cluster import KMeans
import pylab

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

    Co = draw_poly()
    
    doc_arrange, alpha_list, Average_arrange = load_result(doc_num, \
        cancer_type_list, cancer_type_num)

    markers = ['o', '^', 's']

    fig = plt.figure()
    ax = fig.add_subplot(111)
    poly = matplotlib.patches.Polygon(Co, fill=False)
    ax.add_patch(poly)

    plt.xlim(-1.1, 2.5)
    plt.ylim(-1.1, 1.5)

    for i,x in enumerate(Co):
        plt.text(x[0], x[1], str(i+1), fontsize=15)
    for i,x in enumerate(Co):
        left = [0, x[0]]
        height = [0, x[1]]
        plt.plot(left, height, marker='o', linestyle='dashed', markersize='1', \
                markeredgecolor='k', linewidth='0.5', color='k')

    c_palette = sns.color_palette('hls', len(cancer_type_list))
    for i,x in enumerate(cancer_type_list):
        temp_left_list = []
        temp_height_list = []
        for j in range(50):
            theta = np.random.dirichlet(alpha_list[i])
            temp_left = 0
            temp_height = 0
            for k in range(K):
                temp_left += theta[k] * Co[k][0]
                temp_height += theta[k] * Co[k][1]
            temp_left_list.append(temp_left)
            temp_height_list.append(temp_height)
        plt.scatter(temp_left_list, temp_height_list, label=x, c=c_palette[i], \
                marker=markers[i%3], s=5)
    plt.legend()

    name = 'result/data' + args[1] + '/figure/' + args[2] \
            + '_arrangement/alpha_sampling.png'
    fig.tight_layout()
    fig.savefig(name,dpi=300)
    plt.close(1)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    poly = matplotlib.patches.Polygon(Co, fill=False)
    ax.add_patch(poly)

    plt.xlim(-1.1, 2.5)
    plt.ylim(-1.1, 1.5)

    for i,x in enumerate(Co):
        plt.text(x[0], x[1], str(i+1), fontsize=15)
    for i,x in enumerate(Co):
        left = [0, x[0]]
        height = [0, x[1]]
        plt.plot(left, height, marker='o', linestyle='dashed', markersize='1', \
                markeredgecolor='k', linewidth='0.5', color='k')

    index = 0
    for i,x in enumerate(cancer_type_list):
        temp_left_list = []
        temp_height_list = []
        for j in range(cancer_type_num[i]):
            theta = doc_arrange[index]
            temp_left = 0
            temp_height = 0
            for k in range(K):
                temp_left += theta[k] * Co[k][0]
                temp_height += theta[k] * Co[k][1]
            temp_left_list.append(temp_left)
            temp_height_list.append(temp_height)
            index += 1
        plt.scatter(temp_left_list, temp_height_list, label=x, c=c_palette[i], \
                marker=markers[i%3], s=5)
    plt.legend()

    name = 'result/data' + args[1] + '/figure/' + args[2] \
            + '_arrangement/arrangement.png'
    fig.tight_layout()
    fig.savefig(name,dpi=300)
    plt.close(1)

    pca = PCA(n_components=6)
    pca.fit(doc_arrange)
    transformed = pca.fit_transform(doc_arrange)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    index = 0
    for i,x in enumerate(cancer_type_list):
        temp_first_components = []
        temp_second_components = []
        for j in range(cancer_type_num[i]):
            temp_first_components.append(transformed[index][0])
            temp_second_components.append(transformed[index][1])
            index += 1
        ax.scatter(temp_first_components, temp_second_components, label=x, \
                c=c_palette[i], marker=markers[i%3], s=5)
    ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1), borderaxespad=0)
    fig.subplots_adjust(right=0.65)
    name = 'result/data' + args[1] + '/figure/' + args[2] \
            + '_arrangement/pca_1st_2nd.png'
    plt.savefig(name, dpi=300)
    plt.close(1)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    index = 0
    for i,x in enumerate(cancer_type_list):
        temp_first_components = []
        temp_second_components = []
        for j in range(cancer_type_num[i]):
            temp_first_components.append(transformed[index][1])
            temp_second_components.append(transformed[index][2])
            index += 1
        ax.scatter(temp_first_components, temp_second_components, label=x, \
                c=c_palette[i], marker=markers[i%3], s=5)
    ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1), borderaxespad=0)
    fig.subplots_adjust(right=0.65)
    name = 'result/data' + args[1] + '/figure/' + args[2] \
            + '_arrangement/pca_2nd_3rd.png'
    plt.savefig(name, dpi=300)
    plt.close(1)

    for perp in [50, 100]:
        output_file = "result/data" + args[1] + "/figure/" + args[2] +\
                      "_arrangement/tSNE_" + str(perp) + ".txt"
        output = open(output_file, "w")
        transformed = TSNE(n_components=2, n_iter=2000,\
                perplexity=perp).fit_transform(doc_arrange)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        index = 0
        for i,x in enumerate(cancer_type_list):
            temp_first_components = []
            temp_second_components = []
            for j in range(cancer_type_num[i]):
                temp_first_components.append(transformed[index][0])
                temp_second_components.append(transformed[index][1])
                index += 1
            for y,z in zip(temp_first_components, temp_second_components):
                output.write(str(y) + "\t" + str(z) + "\t" + x + "\n")
            ax.scatter(temp_first_components, temp_second_components, label=x,\
                    c=c_palette[i], marker=markers[i%3], s=5)
        ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1), borderaxespad=0)
        fig.subplots_adjust(right=0.65)
        plt.title('tSNE (2D)')
        name = 'result/data' + args[1] + '/figure/' + args[2] \
            + '_arrangement/tSNE_2D_perp_' + str(perp) + '.png'
        plt.savefig(name, dpi=300)
        plt.close(1)

    """
    transformed = TSNE(n_components=3, perplexity=49).fit_transform(doc_arrange)
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    index = 0
    for i,x in enumerate(cancer_type_list):
        temp_first_components = []
        temp_second_components = []
        temp_third_components = []
        for j in range(cancer_type_num[i]):
            temp_first_components.append(transformed[index][0])
            temp_second_components.append(transformed[index][1])
            temp_third_components.append(transformed[index][2])
            index += 1
        ax.scatter3D(temp_first_components, temp_second_components, \
                temp_third_components, label=x, c=c_palette[i], \
                marker=markers[i%3], s=5)
    ax.legend(loc='upper left', bbox_to_anchor=(1.0, 1), borderaxespad=0)
    fig.subplots_adjust(right=0.65)
    name = 'result/data' + args[1] + '/figure/' + args[2] \
            + '_arrangement/tSNE_3D.png'
    plt.savefig(name, dpi=600)
    plt.close(1)
   """

    fig = plt.figure()
    ax1 = fig.add_subplot(111)
    temp_left = range(1, len(cancer_type_list)+1)
    type_dict = {}
    for i,cancer_type in enumerate(cancer_type_list):
        type_dict.update({cancer_type:cancer_type_num[i]})
    height = []
    labels = []
    for key, v in sorted(type_dict.items(), key=lambda x: -x[1]):
        height.append(v)
        labels.append(key)
    ax1.bar(temp_left, height, width=1, align='center')
    ax1.set_ylim(0, max(height))
    ax1.set_xlim(0, len(cancer_type_list)+1)
    ax1.set_xticks(temp_left)
    ax1.set_xticklabels(labels)
    for tick in ax1.get_xticklabels():
        tick.set_rotation(90)
    ax1.tick_params(labelsize=10)
    ax1.set_xlabel('Primary Lesion', fontsize=10)
    Title = '# samples of each primary lesion'
    ax1.set_title(Title, fontsize=10)
    fig.tight_layout()
    name = 'result/data' + args[1] + '/figure/' + args[2] + '_arrangement/' \
            + 'Sample_num.png'
    fig.savefig(name, dpi=200)
    plt.close(1)

    for k in range(K):
        temp_left = range(1, len(cancer_type_list)+1)
        temp_height = {}
        for i,cancer_type in enumerate(cancer_type_list):
            temp_height.update({cancer_type:Average_arrange[i][k]})
        height = []
        labels = []
        for key, v in sorted(temp_height.items(), key=lambda x: -x[1]):
            height.append(v)
            labels.append(key)
        fig = plt.figure()
        ax1 = fig.add_subplot(111)
        ax1.bar(temp_left, height,width=1, align='center')
        ax1.set_ylim(0, max(height))
        ax1.set_xlim(0, len(cancer_type_list)+1)
        ax1.set_xticks(temp_left)
        ax1.set_xticklabels(labels)
        for tick in ax1.get_xticklabels():
            tick.set_rotation(90)
        ax1.tick_params(labelsize=10)
        ax1.set_xlabel('Primary Lesion', fontsize=10)
        ax1.set_ylabel('Average apperance', fontsize=10)
        Title = 'Predicted Signature ' + str(k+1)
        ax1.set_title(Title, fontsize=10)
        fig.tight_layout()
        name = 'result/data' + args[1] + '/figure/' + args[2] \
                + '_arrangement/' + str(k+1) + '_appearance.png'
        fig.savefig(name, dpi=200)
        plt.close(1)

def load_types(data_type):
    File = open('data/PL_data' + data_type + '_list.txt', 'r')
    File.readline()
    cancer_type_list = []
    for line in File.readlines():
        cancer_type_list.append(line[:-1])
    File.close()
    return cancer_type_list

def draw_poly():
    K = int(args[2])
    omega = float(2 * PI / K)
    Co = []
    Co.append([0.0, 1.0])
    for i in range(1, K):
        Co.append([np.cos(omega) * Co[i-1][0] + np.sin(omega) * Co[i-1][1], \
            -np.sin(omega) * Co[i-1][0] + np.cos(omega) * Co[i-1][1]])
    return Co

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
