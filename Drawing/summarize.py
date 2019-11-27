import os
import sys
import matplotlib.pyplot as plt
import seaborn as sns

def main():
    dic = {}
    for i in range(first, last+1):
        path = "result/data" + str(i) + "/figure/"
        files = os.listdir(path)
        files_dir = [f for f in files if os.path.isdir(os.path.join(path, f))]
        file_name = files_dir[0]
        k = int(file_name[:file_name.index("_")])
        if k in dic:
            dic[k] += 1
        else:
            dic.update({k:1})
    result_path = "result/data" + str(first) + "-" + str(last) + "/"
    if not os.path.isdir(result_path):
        os.mkdir(result_path)
    fig = plt.figure()
    ax = fig.add_subplot(111)
    left = []; height = []
    for k,v in sorted(dic.items(), key=lambda x:x[0]):
        left.append(k)
        height.append(v)
    ax.bar(left, height, width=1, align="center")
    ax.set_xlim(left[0]-1, left[len(left)-1]+1)
    ax.set_xticks(left)
    plt.tick_params(labelsize=18)
    plt.xlabel("The number of signatures: K", fontsize=18)
    Title = "Model selection (with PLDA)"
    ax.set_title(Title, fontsize=18)
    fig.tight_layout()
    name = result_path + "Predicted_K.png"
    fig.savefig(name, dpi=200)
    plt.close(1)

if __name__ == "__main__":
    args = sys.argv
    first = int(args[1])
    last = int(args[2])
    main()
