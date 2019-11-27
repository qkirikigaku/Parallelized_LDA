import sys
import numpy as np
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    mat,labels = load_data()
    pred = KMeans(n_clusters=num_clusters).fit_predict(mat)
    draw_kmeans(mat, pred)
    draw_kmeans_cluster_lesions(pred, labels)

def load_data():
    input_file = "result/data" + data_type + "/figure/" +\
                 num_sig + "_arrangement/tSNE_50.txt"
    lines = open(input_file, "r").readlines()
    mat = np.zeros([len(lines), 2])
    labels = list()
    for i,line in enumerate(lines):
        temp_list = line.split()
        mat[i,0] = float(temp_list[0])
        mat[i,1] = float(temp_list[1])
        labels.append(temp_list[2])
    return mat, labels

def draw_kmeans(mat, pred):
    fig = plt.figure()
    c_palette = sns.color_palette("hls", num_clusters)
    for i,x in enumerate(mat):
        plt.scatter(x[0],x[1],c=c_palette[pred[i]])
    plt.title("tSNE result (KMeans clustering)")
    fig_name = "result/data" + args[1] + "/figure/" + args[2] +\
               "_arrangement/tSNE_50_KMeans.png"
    plt.savefig(fig_name, dpi=300)    
    plt.close(1)

def draw_kmeans_cluster_lesions(pred, labels):
    lesions = [{} for i in range(num_clusters)]
    sum_sample = np.zeros([num_clusters])
    for i,x in enumerate(pred):
        if(labels[i] not in lesions[x].keys()):
            lesions[x].update({labels[i]:1})
        else:
            lesions[x].update({labels[i]:(lesions[x][labels[i]]+1)})
        sum_sample[x] += 1
    for i,x in enumerate(lesions):
        fig = plt.figure()
        ax = fig.add_subplot(111)
        left = range(len(x)); les = list(); height = list()
        for k,v in sorted(x.items(), key=lambda y: -y[1]):
            les.append(k)
            height.append(v)
        ax.bar(left, height, width=1, align="center")
        ax.set_ylim(0, max(height))
        ax.set_xlim(-1, len(x))
        ax.set_xticks(left)
        ax.set_xticklabels(les)
        for tick in ax.get_xticklabels():
            tick.set_rotation(90)
        ax.tick_params(labelsize=10)
        ax.set_xlabel("Primary Lesion", fontsize=10)
        ax.set_ylabel("# samples", fontsize=10)
        Title = "Cluster " + str(i+1) + " (" + str(int(sum_sample[i])) + " samples)"
        ax.set_title(Title, fontsize=10)
        fig.tight_layout()
        name = "result/data" + data_type + "/figure/" + num_sig +\
               "_arrangement/kmeans/Cluster" + str(i+1) +\
               "_lesion.png"
        fig.savefig(name,dpi=200)
        plt.close(1)

if __name__ == "__main__":
    args = sys.argv
    data_type = args[1]
    num_sig = args[2]
    num_clusters = 6
    main()
