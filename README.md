Please refer to our paper: [Parallelized Latent Dirichlet Allocation Provides a Novel Interpretability of Mutation Signatures in Cancer Genomes](https://www.mdpi.com/2073-4425/11/10/1127) for more detailed information.

# Parallelized_LDA
Parallelized-LDA can extract mutation signatures contained in samples from mutation catalogs of cancer genome with Bayesian model.
We named our model as Parallelized-LDA, an extended model of Latent Dirichelt Alloctaion (LDA).
This model can detect all of the mutation signatures simultaneously across the primary lesions.
Please refer to our paper for details about our model.

## Requirements
The source codes are written in C++ and Python3.
* g++
* boost (C++ Library)
* make

## Compile
Execute bellow command to compile.
```
make compile
```
Then, your executable file `MS` is created in `bin/` directory.

## Inputs
We have prepared a pre-processed mutation catalog for analysis in `data/` directory and you need not touch it.
If you want to use your own mutation catalog, you should prepare the dataset that meets the following standards:
* `data/PL_data${data_type}_list.txt`
```
L
${primary_lesion1}
${primary_lesion2}
...
${primari_lesionL}
```
${data_type} is corresponding to which type of mutation dictionary using. About mutation dictionary, please refer to our paper.
*L* shows the number of primary lesions.

* `data/data${data_type}_${primary_lesion}.txt`
This file is necessary for every primary lesion.
```
N V
M_11 M_12 ...  M_1V
M_21 M_22 ...  M_2V
...
M_N1 M_N2 ...  M_NV
```
*N* and *V* show # samples and # mutational types, respectively.
*M_ij* shows the number of *j*th mutation in *i*th sample.


## Discovering mutation signatures
Execute bellow command to detect mutation signatures under `Parallelized_LDA/` directory.
```
bash scripts/MS_real.sh ${data_type}
```
`${data_type}` is the command line argument and it shows the type of mutation dictionary (we have implemented for `data_type = 1 or 2`).
About mutation dictionary, please refer to our paper in detail.

By executing the above command, you can visualized the results all at once.

## Outputs
You can see the results of predicted signatures in `result/data${data_type}` directory.
* `figure/` : all of the figures are put in here.
	* `ELBO.png`, `ELBO_a.png`, `ELBO_b.png` : these figures show the variational lower bound which is the objective function of Variational Bayes method.
	* `${k}_arrangement/` : you can see the results about activity of each signature.
	* `${k}_lesion/` : you can see the results about activity of each primary lesion.
	* `${k}_signature/` : you can see the results about predicted signatures.
* `result_k${k}.txt` : result with text format when setting the number of topic to ${k}.

## Changing the options
By editing the variables in `scripts/MS_real.sh`, you can change the options to suit your needs.
* `${Iter}` (default : 100)
	* *Iter* variable is the numebr of iterations of experiments. As you know, Variational Bayes methods depends on initial values of the parameters to learn the space of latent variables, so we re-allocate it for *Iter* times. The computational time of `MS_real.sh` takes about linear order of variable *Iter*.

* `${K}` (default : 30)
	* *K* variable is the max number of topics in experiments. If you expect more topics than default value (30) in the entire catalogs, please set this variable to the large value. The computational time of `MS_real.sh` takes about linear order of variable *K*.
