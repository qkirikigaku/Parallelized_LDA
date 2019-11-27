import sys
import numpy as np
import matplotlib.pyplot as plt

def main():
    args = sys.argv
    """args[1]:data_type"""
    K = int(args[2])

    files = (K-1)*[0]
    i = 2
    while(i < K+1):
        string = 'result/data' + args[1] + '/result_k'
        if(i <= 9):
            string += '0' + str(i) + '.txt'
            files[i-2] = string
        else:
            string += str(i) + '.txt'
            files[i-2] = string
        i += 1;

    ELBO = np.zeros([K-1])

    labels = (K-1)*[0]
    for i in range(K-1):
        data = open(files[i])
        elbo = data.readline()
        elbo = elbo.replace("\n","")
        ELBO[i] = float(elbo)
        labels[i] = 'topic' + str(i+2)

    max_el = 0
    semi_el = 0
    for i in range(K-1):
        if(ELBO[i] > ELBO[max_el]):
            semi_el = max_el
            max_el = i
        elif(ELBO[i] > ELBO[semi_el]):
            semi_el = i

    el_colors = (K-1)*["b"]
    el_colors[max_el] = "r"
    el_colors[semi_el] = "g"

    left = np.arange(1,K,1)
    height = ELBO.copy()

    plt.bar(left,height, align = "center", color = el_colors)
    plt.title("Variational lower bound of each signature number")
    plt.xlabel("signature numbers")
    plt.ylabel("Variational lower bound")
    plt.xticks(left, labels, rotation = 90, fontsize = "small")
    plt.tight_layout()
    name = 'result/data' + args[1] + '/figure/ELBO.png'
    plt.savefig(name)
    plt.close(1)

    plt.figure()

    plt.bar(left,height, align = "center", color = el_colors)
    plt.title("Variational lower bound of each signature number")
    plt.xlabel("signature numbers")
    plt.ylabel("Variational lower bound")
    plt.xticks(left, labels, rotation = 90, fontsize = "small")
    y_max =max(height) * 0.99
    y_min =min(height) * 1.01
    plt.ylim(ymax = y_max, ymin = y_min)
    plt.tight_layout()
    name = 'result/data' + args[1] + '/figure/ELBO_a.png'
    plt.savefig(name)
    plt.close(1)

    plt.figure()
    plt.bar(left[max_el-10:], height[max_el-10:], align = 'center', color = el_colors[max_el-10:])
    plt.title("Variational lower bound of each signature number")
    plt.xlabel("signature numbers")
    plt.ylabel("Variational lower bound")
    plt.xticks(left[max_el-10:], labels[max_el-10:], rotation = 90, fontsize = "small")
    y_max =max(height) + 100
    y_min =min(height[max_el-10:]) - 100
    plt.ylim(ymax = y_max, ymin = y_min)
    plt.tight_layout()
    name = 'result/data' + args[1] + '/figure/ELBO_b.png'
    plt.savefig(name)
    plt.close(1)


if __name__ == '__main__':
    main()
