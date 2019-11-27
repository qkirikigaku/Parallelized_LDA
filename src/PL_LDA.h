//
// Created by Taro Matsutani on 2019/11/27.
//

#ifndef LDA_LDA_H
#define LDA_LDA_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <random>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/gamma.hpp>


using namespace std;

class PL_LDA {
    int num_topic;
    int num_vocabulary;
    vector<int> num_document;
    int num_PL;

    //document[L:lth lesion][d:dth document][i:ith word]
    vector<vector<vector<int> > > train_document;

    int data_type; // 1,2,3,4.
    vector<string> cancer_type_list; // primary lesion.
    int iteration;

    vector<vector<vector<vector<double> > > > log_responsibility;
    //q(z)　= responsibility[l:lth lesion][D:Dth document][I:Ith word][K:Kth topic]

    vector<vector<vector<double> > > xi_theta;
    //xi_theta[l:lth lesion][D:Dth document][K::Kth topic]
    vector<vector<double> > xi_phi;
    //xi_phi[K:Kth topic][V:Vth word]

    vector<vector<double> > alpha;
    //alpha[l:lth lesion][K:Kth topic]
    vector<double> beta;
    //beta[V:Vth word]

    vector<vector<vector<double> > > log_sum_res_for_i;
    //E[n_{ldk}]
    vector<vector<double> > log_sum_res_for_di;
    //E[n_{kv}]

    double temp_variational_lower_bound;
    double old_variational_lower_bound;

public:
    PL_LDA(int x, int y, int z);
    
    void Normalize(vector<double> &vector, int &length);
    double log_sum_exp(vector<double> &vec, int &length);
    
    void run_VB();
    void initialize();
    
    void Update_log_responsibility();
    void Update_parameter();
    void Update_hyperparameter();
    void calc_variational_lower_bound();
    void calc_n_dk();
    void calc_n_kv();

    void load_cancer_types();
    void load_data(string cancer_type, int index);
    void write_data();
    void show_ELBO();
};
    
PL_LDA::PL_LDA(int x, int y, int z) {
    num_topic = x;
	data_type = y;
	iteration = z;
}

void PL_LDA::Normalize(vector<double> &vec, int &length){
    double sum = log_sum_exp(vec, length);
    for (int i=0; i < length; i++){
        vec[i] -= sum;
    }
}

double PL_LDA::log_sum_exp(vector<double> &vec, int &length){
    double max_iter = 0;
    for (int i=1; i < length; i++){
        if(vec[i] > vec[max_iter])max_iter = i;
    }
    double sum = 0;
    for (int i=0; i < length; i++){
        sum += exp(vec[i] - vec[max_iter]);
    }
    double return_value = vec[max_iter] + log(sum);
    return(return_value);
};

void PL_LDA::run_VB() {
    
    for (int i=0; i < num_PL; i++){
        load_data(cancer_type_list[i], i);
    }
    initialize();
    old_variational_lower_bound = -10000000000000;
    Update_parameter();
    for (int i = 0; i < 1000; i++) {
        cout << "iter:" << i << endl;
        Update_log_responsibility();
        Update_parameter();
        Update_hyperparameter();
        calc_variational_lower_bound();
        if (fabs(temp_variational_lower_bound - old_variational_lower_bound) < 1) {
            show_ELBO();
            cout << endl;
            break;
        } else {
            show_ELBO();
        }
        cout << endl;
    }
};

void PL_LDA::initialize(){
    int l, d, i, k, v;

	random_device rnd;
    log_responsibility.resize(num_PL);
    for (l=0; l < num_PL; l++){
        log_responsibility[l].resize(num_document[l]);
        for (d = 0; d < num_document[l]; d++) {
            log_responsibility[l][d].resize(train_document[l][d].size());
            for (i = 0; i < train_document[l][d].size(); i++) {
                log_responsibility[l][d][i].resize(num_topic, 0);
                for (k=0; k < num_topic; k++){
                    log_responsibility[l][d][i][k] = (double) -rnd();
                }
                Normalize(log_responsibility[l][d][i], num_topic);
            }
        }
    }

    mt19937 mt(rnd());
    uniform_real_distribution<double> Uniform(0.0, 1.0);
    alpha.resize(num_PL);
    for (l=0; l < num_PL; l++){
        alpha[l].resize(num_topic);
        for (k = 0; k < num_topic; k++) {
            alpha[l][k] = Uniform(mt);
        }
    }

    beta.resize(num_vocabulary);
    for (v = 0; v < num_vocabulary; v++) {
        beta[v] = Uniform(mt);
    }

    xi_theta.resize(num_PL);
    for (l=0; l < num_PL; l++){
        xi_theta[l].resize(num_document[l]);
        for (d = 0; d < num_document[l]; d++) {
            xi_theta[l][d].resize(num_topic, 0);
        }
    }

    xi_phi.resize(num_topic);
    for (k = 0; k < num_topic; k++) {
        xi_phi[k].resize(num_vocabulary, 0);
    }

    log_sum_res_for_i.resize(num_PL);
    for (l=0; l < num_PL; l++){
        log_sum_res_for_i[l].resize(num_document[l]);
        for (d = 0; d < num_document[l]; d++) {
            log_sum_res_for_i[l][d].resize(num_topic, 0);
        }
    }

    log_sum_res_for_di.resize(num_topic);
    for (k = 0; k < num_topic; k++) {
        log_sum_res_for_di[k].resize(num_vocabulary, -730);
    }
};

void PL_LDA::Update_log_responsibility() {
    int l, d, v, i, k;

    vector<vector<double> > sum_xi_theta;
    sum_xi_theta.resize(num_PL);
    for (l=0; l < num_PL; l++){
        sum_xi_theta[l].resize(num_document[l], 0);
        for (d = 0; d < num_document[l]; d++) {
            for (k = 0; k < num_topic; k++) {
                sum_xi_theta[l][d] += xi_theta[l][d][k];
            }
        }
    }

    vector<double> sum_xi_phi;
    sum_xi_phi.resize(num_topic, 0);
    for (k = 0; k < num_topic; k++) {
        for (v = 0; v < num_vocabulary; v++) {
            sum_xi_phi[k] += xi_phi[k][v];
        }
    }

    for (l=0; l < num_PL; l++){
        for (d = 0; d < num_document[l]; d++) {
            for (i = 0; i < train_document[l][d].size(); i++) {
                for (k = 0; k < num_topic; k++) {
                    log_responsibility[l][d][i][k] =
                        boost::math::digamma(xi_phi[k][train_document[l][d][i]]) + boost::math::digamma(xi_theta[l][d][k])
                        - boost::math::digamma(sum_xi_phi[k]) - boost::math::digamma(sum_xi_theta[l][d]);
                }
            }
        }
    }

    for (l=0; l < num_PL; l++){
        for (d=0; d < num_document[l]; d++){
            for (i=0; i < train_document[l][d].size(); i++){
                Normalize(log_responsibility[l][d][i], num_topic);
            }
        }
    }
};

void PL_LDA::Update_parameter() {
    int l, d, v, i, k;
    //Update xi_theta
    calc_n_dk();
    for (l=0; l < num_PL; l++){
        for (d = 0; d < num_document[l]; d++) {
            for (k = 0; k < num_topic; k++) {
                xi_theta[l][d][k] = exp(log_sum_res_for_i[l][d][k]) + alpha[l][k];
            }
        }
    }
    //Update xi_phi
    calc_n_kv();
    for (k = 0; k < num_topic; k++) {
        for (v = 0; v < num_vocabulary; v++) {
            xi_phi[k][v] = exp(log_sum_res_for_di[k][v]) + beta[v];
        }
    }
};

void PL_LDA::Update_hyperparameter() {
    int l, d, v, i, k;

    calc_n_dk();

    //Update alpha
    vector<vector<double> > new_numerator_alpha;
    new_numerator_alpha.resize(num_PL);
    for (l=0; l < num_PL; l++){
        new_numerator_alpha[l].resize(num_topic, 0);
        for (k = 0; k < num_topic; k++) {
            for (d = 0; d < num_document[l]; d++) {
                new_numerator_alpha[l][k] += (boost::math::digamma(exp(log_sum_res_for_i[l][d][k]) + alpha[l][k]) - boost::math::digamma(alpha[l][k])) * alpha[l][k];
            }
        }
    }

    vector<double> new_denominator_alpha,sum_alpha;
    new_denominator_alpha.resize(num_PL,0);
    sum_alpha.resize(num_PL, 0);
    for (l=0; l < num_PL; l++){
        for (k = 0; k < num_topic; k++) {
            sum_alpha[l] += alpha[l][k];
        }
    }
    for (l=0; l < num_PL; l++){
        for (d = 0; d < num_document[l]; d++) {
            new_denominator_alpha[l] += boost::math::digamma(train_document[l][d].size() + sum_alpha[l]) - boost::math::digamma(sum_alpha[l]);
        }
    }

    vector<int> count;
    count.resize(num_PL, 0);
    for (l=0; l < num_PL; l++){
        for (k = 0; k < num_topic; k++) {
            if (new_numerator_alpha[l][k] == 0) {
                count[l] += 1;
            }
        }
    }
    for (l=0; l < num_PL; l++){
        if (count[l] == 0) {
            for (k = 0; k < num_topic; k++) {
                alpha[l][k] = new_numerator_alpha[l][k] / new_denominator_alpha[l];
            }
        }
    }

    //Update beta

    //Update beta for each word
    vector<double> new_numerator_beta;
    vector<double> new_denominator_beta;
    new_numerator_beta.resize(num_vocabulary,0);
    new_denominator_beta.resize(num_vocabulary,0);

    calc_n_kv();
    for(v=0; v<num_vocabulary; v++){
        for(k=0; k<num_topic; k++){
            new_numerator_beta[v] += (boost::math::digamma(exp(log_sum_res_for_di[k][v])+beta[v])-boost::math::digamma(beta[v]))*beta[v];
        }
    }

    vector<double> sum_res_for_div;
    double sum_beta = 0;
    sum_res_for_div.resize(num_topic,0);
    for(k=0; k<num_topic; k++){
        for(v=0; v<num_vocabulary; v++){
            sum_res_for_div[k] += exp(log_sum_res_for_di[k][v])+beta[v];
        }
    }
    for(v=0; v<num_vocabulary; v++){
        sum_beta += beta[v];
    }
    for(v=0; v<num_vocabulary; v++){
        for(k=0; k<num_topic; k++){
            new_denominator_beta[v] += boost::math::digamma(sum_res_for_div[k])-boost::math::digamma(sum_beta);
        }
    }
    for(v=0; v<num_vocabulary; v++){
        beta[v] = new_numerator_beta[v] / new_denominator_beta[v];
    }
};

void PL_LDA::calc_variational_lower_bound() {
    int l, d, v, i, k;
    double first_component = 0, second_component = 0, third_component = 0, fourth_component = 0, fifth_component = 0;

    //calc_first_component
    double Vbeta = 0;
    double sum_lgamma_beta = 0;
    for(v=0; v<num_vocabulary; v++){
        Vbeta += beta[v];
        sum_lgamma_beta += boost::math::lgamma(beta[v]);
    }
    for (k = 0; k < num_topic; k++) {
        double sum_xi_phi_for_v = 0, log_pi_gamma_xi_phi_for_v = 0;
        for (v = 0; v < num_vocabulary; v++) {
            sum_xi_phi_for_v += xi_phi[k][v];
            log_pi_gamma_xi_phi_for_v += boost::math::lgamma(xi_phi[k][v]);
        }
        first_component += boost::math::lgamma(Vbeta)
                        - sum_lgamma_beta
                        - boost::math::lgamma(sum_xi_phi_for_v)
                        + log_pi_gamma_xi_phi_for_v;
    }

    //calc_second_component
    calc_n_kv();
    for (k = 0; k < num_topic; k++) {
        double sum_xi_phi_for_v = 0;
        for (v = 0; v < num_vocabulary; v++) {
            sum_xi_phi_for_v += xi_phi[k][v];
        }
        for (v = 0; v < num_vocabulary; v++) {
            second_component += (exp(log_sum_res_for_di[k][v]) + beta[v] - xi_phi[k][v]) *
                                (boost::math::digamma(xi_phi[k][v]) - boost::math::digamma(sum_xi_phi_for_v));
        }
    }

    //calc_third_component
    vector<double> sum_alpha,log_sum_gamma_alpha;
    sum_alpha.resize(num_PL, 0); log_sum_gamma_alpha.resize(num_PL, 0);
    
    for (l=0; l < num_PL; l++){
        for (k = 0; k < num_topic; k++) {
            sum_alpha[l] += alpha[l][k];
            log_sum_gamma_alpha[l] += boost::math::lgamma(alpha[l][k]);
        }
    }
    for (l=0; l < num_PL; l++){
        for (d = 0; d < num_document[l]; d++) {
            double sum_xi_theta_for_k = 0, log_sum_gamma_xi_theta_for_k = 0;
            for (k = 0; k < num_topic; k++) {
                sum_xi_theta_for_k += xi_theta[l][d][k];
                log_sum_gamma_xi_theta_for_k += boost::math::lgamma(xi_theta[l][d][k]);
            }
            third_component += boost::math::lgamma(sum_alpha[l]) - log_sum_gamma_alpha[l] - boost::math::lgamma(sum_xi_theta_for_k) + log_sum_gamma_xi_theta_for_k;
        }
    }

    //calc_fourth_component
    calc_n_dk();
    for (l=0; l < num_PL; l++){
        for (d = 0; d < num_document[l]; d++) {
            double sum_xi_theta_for_k = 0;
            for (k = 0; k < num_topic; k++) {
                sum_xi_theta_for_k += xi_theta[l][d][k];
            }
            for (k = 0; k < num_topic; k++) {
                fourth_component += (exp(log_sum_res_for_i[l][d][k]) + alpha[l][k] - xi_theta[l][d][k]) *
                                    (boost::math::digamma(xi_theta[l][d][k]) - boost::math::digamma(sum_xi_theta_for_k));
            }
        }
    }

    //calc_fifth_component
    for (l=0; l < num_PL; l++){
        for (d = 0; d < num_document[l]; d++) {
            for (i = 0; i < train_document[l][d].size(); i++) {
                for (k = 0; k < num_topic; k++) {
                    fifth_component += exp(log_responsibility[l][d][i][k]) * log_responsibility[l][d][i][k];
                }
            }
        }
    }

    temp_variational_lower_bound = first_component + second_component + third_component + fourth_component - fifth_component;
};
    
void PL_LDA::calc_n_dk() {
    int l, d, v, i, k;
    vector<vector<vector<vector<double> > > > temp_log_res;
    temp_log_res.resize(num_PL);
    for (l=0; l < num_PL; l++){
        temp_log_res[l].resize(num_document[l]);
        for (d=0; d < num_document[l]; d++){
            temp_log_res[l][d].resize(num_topic);
            for (k=0; k < num_topic; k++){
                temp_log_res[l][d][k].resize(train_document[l][d].size(),0);
            }
            for (i=0; i < train_document[l][d].size(); i++){
                for (k=0; k < num_topic; k++){
                    temp_log_res[l][d][k][i] = log_responsibility[l][d][i][k];
                }
            }
        }
    }
    for (l=0; l < num_PL; l++){
        for (d = 0; d < num_document[l]; d++) {
            int doc_length = train_document[l][d].size();
            for (k=0; k < num_topic; k++){
                log_sum_res_for_i[l][d][k] = log_sum_exp(temp_log_res[l][d][k], doc_length);
            }    
        }
    }
};

void PL_LDA::calc_n_kv() {
    int l, d, v, i, k;
    
    vector<vector<vector<double> > > temp_log_res;
    temp_log_res.resize(num_topic);
    for (k=0; k < num_topic; k++){
        temp_log_res[k].resize(num_vocabulary);
        for (l=0; l < num_PL; l++){
            for (d=0; d < num_document[l]; d++){
                for (i=0; i < train_document[l][d].size(); i++){
                    temp_log_res[k][train_document[l][d][i]].push_back(log_responsibility[l][d][i][k]);
                }
            }
        }
    }

    for (k=0; k < num_topic; k++){
        for(v=0; v < num_vocabulary; v++){
            int vec_size = temp_log_res[k][v].size();
            log_sum_res_for_di[k][v] = log_sum_exp(temp_log_res[k][v], vec_size);
        }
    }
};

void PL_LDA::load_cancer_types(){
    ifstream ifs;
    string input_filename = "data/PL_data" + to_string(data_type) + "_list.txt";  
    ifs.open(input_filename.c_str(), ios::in);
    if(!ifs) {
        cout << "Cannot open " + input_filename << endl;
        exit(1);
    }
    char buf[100000];
    ifs.getline(buf, 100000);
    num_PL = atoi(buf);
    cancer_type_list.resize(num_PL);
    for (int i=0; i < num_PL; i++){
        ifs.getline(buf, 100000);
        cancer_type_list[i] = buf;
    }
    ifs.close();
    num_document.resize(num_PL);
    train_document.resize(num_PL);
};

void PL_LDA::load_data(string cancer_type, int index) {
    ifstream ifs;
    string input_file_name = "data/data" + to_string(data_type) + "_" + cancer_type + ".txt";
    ifs.open(input_file_name.c_str(), ios::in);
    if (!ifs) {
        cout << "Cannnot open " + input_file_name << endl;
        exit(1);
    }
    char buf[1000000];
    char *temp;
    vector<vector<int> > raw_document;
    vector<int> words_number;
    ifs.getline(buf, 1000000);
    temp = strtok(buf, " ");
    num_document[index] = atoi(temp);
    raw_document.resize(num_document[index]);
    words_number.resize(num_document[index],0);
    train_document[index].resize(num_document[index]);
    temp = strtok(NULL, " ");
    num_vocabulary = atoi(temp);
    int temp_word_number;

    for (int d = 0; d < num_document[index]; d++) {
        ifs.getline(buf, 1000000);
        for (int v = 0; v < num_vocabulary; v++) {
            if (v == 0) temp_word_number = atoi(strtok(buf, " "));
            else temp_word_number = atoi(strtok(NULL, " "));
            for(int i = 0; i < temp_word_number; i++){
                raw_document[d].push_back(v);
                words_number[d]++;
            }
        }
    }
	
    for (int d = 0; d < num_document[index]; d++){
	    train_document[index][d].resize(words_number[d]);
        for(int i = 0; i < words_number[d]; i++){
	        train_document[index][d][i] = raw_document[d][i];
        }
    }
    ifs.close();
};

void PL_LDA::write_data() {
    ofstream ofs;
    string output_file_name = "result/data" + to_string(data_type) + "_" + to_string(iteration) + "/result_k";
    if(num_topic < 10){
        output_file_name += "0" + to_string(num_topic) + ".txt";
    }
    else{
        output_file_name += to_string(num_topic) + ".txt";
    }
    ofs.open(output_file_name, ios::out);
    ofs << to_string(temp_variational_lower_bound) << "\n";
    ofs << "0" << "\n";
    
    int l, d, i, k, v;
    vector<vector<double> > Enkv;
    vector<double> sum_output;
    Enkv.resize(num_topic);
    sum_output.resize(num_topic, 0);
    for (k = 0; k < num_topic; k++) {
        Enkv[k].resize(num_vocabulary);
        for (v = 0; v < num_vocabulary; v++){
            sum_output[k] += exp(log_sum_res_for_di[k][v]);
        }
        for (v = 0; v < num_vocabulary; v++) {
            Enkv[k][v] = exp(log_sum_res_for_di[k][v]) / sum_output[k];
            ofs << to_string(Enkv[k][v]) << " ";
        }
        ofs << "\n";
    }

	vector<vector<vector<double> > > Endk;
    vector<vector<double> > sum_output_a;
    Endk.resize(num_PL);
    sum_output_a.resize(num_PL);
    for(l=0; l < num_PL; l++){
        Endk[l].resize(num_document[l]);
        sum_output_a[l].resize(num_document[l], 0);
        for(d = 0; d < num_document[l]; d++){
            Endk[l][d].resize(num_topic, 0);
            for(k = 0; k < num_topic; k++){
                sum_output_a[l][d] += exp(log_sum_res_for_i[l][d][k]);
            }
            for(k = 0; k < num_topic; k++){
                Endk[l][d][k] = exp(log_sum_res_for_i[l][d][k]) / sum_output_a[l][d];
                ofs << to_string(Endk[l][d][k]) << " ";
            }
            ofs << "\n";
        }
    }

    for (l=0; l < num_PL; l++){
	    for (k = 0; k < num_topic; k++){
	        ofs << alpha[l][k] << " ";
	    }
	    ofs << "\n";
    }
    ofs.close();
};

void PL_LDA::show_ELBO() {
    cout << "ELBO: " << temp_variational_lower_bound << endl;
    cout << "Improvement point: " << temp_variational_lower_bound - old_variational_lower_bound << endl;
    old_variational_lower_bound = temp_variational_lower_bound;
};

#endif //LDA_LDA_H
