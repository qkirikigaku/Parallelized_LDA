import sys
import os

K = 5

def main():
    cos_dist = []
    for i in range(start, end+1):
        path = "result/data" + str(i) + "/figure/"
        files = os.listdir(path)
        files_dir = [f for f in files if os.path.isdir(os.path.join(path, f))]
        file_name = files_dir[0]
        k = int(file_name[:file_name.index("_")])
        if(k == 10):
            result_path = path + str(k) + "_signature/"
            cos_file = result_path + "minimum_cos.txt"
            lines = open(cos_file, "r").readlines()
            for line in lines:
                cos_dist.append(float(line))
    mean_val = sum(cos_dist)/len(cos_dist)
    print("PLDA: " + str(mean_val))

if __name__ == "__main__":
    args = sys.argv
    start = int(args[1])
    end = int(args[2])
    main()
