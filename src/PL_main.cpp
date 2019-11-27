#include "PL_LDA.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

void run_VB_LDA(int num_topic, int data_type, int iteration);

int main(int argc,char *argv[]){
    if(argc != 4){
        cout << "The number of argument is invalid." << endl;
        return(0);
    }
    int number_of_topic = atoi(argv[1])+1;
    int data_type = atoi(argv[2]); // = (int) 1,2,3,4.
    int iteration = atoi(argv[3]); // = (int) num. To avoid local minimum.
    run_VB_LDA(number_of_topic, data_type, iteration);
}

void run_VB_LDA(int num_topic, int data_type, int iteration){
    PL_LDA temp_object(num_topic, data_type, iteration);
    temp_object.load_cancer_types();
    temp_object.run_VB();
    temp_object.write_data();
}
