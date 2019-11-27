import sys
import os
from multiprocessing import Pool
import subprocess


def main(): 
    result_path = "result/data" + data_type + "_" + experiment
    if(os.path.exists(result_path) == False): os.mkdir(result_path)
    arguments = []
    for j in range(1, K): # j shows the number of topics.
        arguments.append((str(j), data_type, experiment))
    pool = Pool()
    _ = pool.starmap(execute, arguments)


def execute(num_topic, data_type, experiment):
    cmd = 'bin/MS ' + num_topic + ' ' + data_type + ' ' + experiment
    subprocess.call(cmd.split())


if __name__ == '__main__':
    args = sys.argv
    data_type = args[1]
    experiment = args[2]
    K = int(args[3])
    main()
