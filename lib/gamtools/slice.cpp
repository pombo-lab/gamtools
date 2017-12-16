#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<string.h>
#include<cmath>
#include<sstream>
#include "slice.h"

#define m 204
#define n_p 1
#define b 30000
#define h 0.22
#define R 4.5
#define L 5000000000

#define file_tube  "matrix.out"
#define file_chr_names "chr_names.out" 
#define file_chr_indeces "chr_indeces.out" 

#define file_pi_out "./output/pi/pi_pairs"
#define file_out_pi_thr "./output/pi_thresholds/pi_thresholds"

#define under 0
#define over 100
#define tol_ext 5e-2
#define thr_u0 100
#define quant_value 0.95 
#define N_draws 1e4
#define single_chr -1
#define pi_out sqrt(-2) 

using namespace std;

int run_slice(){
    
#pragma mark DEFINE VARIABLES

    ifstream is;
    ofstream os;
    stringstream ss;
    ofstream os_pi;            

    int n_columns, n_rows, temp_int, n_chromosomes, n_loci_temp, n_chromo_to_analyse, index, index_chromo, ind_locus;
    int locus;
    int first_locus, last_locus;
    int n2, n1a, n1b, n1, n0;    
    bool flag;

    long int *nused, nused_tot, nloci_tot;
    long int max_dist_u0, g_dist;   
    long int *n_dist, *n_tot_dist, *n_chromo_dist, max_dist;
    long double *avg_m1s, *var_m1s, avg_m1s_tot, var_m1s_tot; 
    long double *avg_m2_dist, *avg_m1a_dist, *avg_m1b_dist, *avg_m0_dist, *avg_ratio_couples, *var_m2_dist, *var_m1a_dist, *var_m1b_dist, *var_m1_dist, *var_m0_dist, *var_ratio_couples;
    long double *eps;
    long double eps_value, pi, diff_exp_theo, ratio, u0_value;

    char temp_file_name [1000];

    vector< vector<bool> > matrix;
    vector<string> chromosomes, chromo_to_analyse;
    vector<unsigned int> start_ind, end_ind, chromo_length;
    vector<int> ind_chromo_to_analyse;
    vector<long double> u0, eps_vector, quantiles, pi_thr;
    vector<unsigned int> distances_u0;

    string temp_string;
    string temp_file_name_2;  

    double eps_fix = 0.;
    double tol = 1e-10;  
    double eps_fix_pi = 0.9;
    double tol_pi = 1e-6;     
      
#pragma CALCULATE PARAMETER VALUES	

    const long double bL=(long double)b/(long double)L;
    const long double heff_R=((long double)(h)/(long double)(R))+2*pow((long double)(bL),(long double)(1./3.));
    const long double v0=2./(2.+(heff_R));
 
    const double v0_value=2./(2.+(heff_R));

#pragma mark LOAD THE MATRIX OF TUBES 
    
    cout << "Loading matrix: " << file_tube << endl;    
    matrix=load_matrix_tubes(file_tube, n_rows, n_columns, m);    
    cout << "Rows (loci) and columns (tubes) counted (" << n_rows << " " << n_columns << ")" << endl;
    cout << "Matrix loaded!" << endl;
    cout << "************" << endl;
    
#pragma mark LOAD CHROMOSOMES INFO 
    
    cout << "Loading chromosomes info..." << endl;
    load_chromo_info(chromosomes, chromo_to_analyse, start_ind, end_ind, chromo_length, ind_chromo_to_analyse, n_chromo_to_analyse, file_chr_names, file_chr_indeces, single_chr);
    cout << "Analysing chromosomes: " << endl;    
    for(int i=0; i<ind_chromo_to_analyse.size(); i++){    
        cout << chromo_to_analyse[i] << " " << ind_chromo_to_analyse[i] << endl;
    }    
    cout << "Chromosome information loaded!" << endl;
    cout << "**********" << endl;
        
#pragma mark COMPUTE m1s and eps 
    
    avg_m1s=(long double *)(calloc(n_chromo_to_analyse, sizeof(long double)));
    var_m1s=(long double *)(calloc(n_chromo_to_analyse, sizeof(long double)));
    nused=(long int *)(calloc(n_chromo_to_analyse, sizeof(long int)));
    eps=(long double *)(calloc(n_chromo_to_analyse+1, sizeof(long double)));
    
    compute_m1s_eps(avg_m1s, var_m1s, nused, eps, avg_m1s_tot, var_m1s_tot,nused_tot, nloci_tot, n_chromo_to_analyse, ind_chromo_to_analyse,    start_ind,  chromo_length,  m, matrix, heff_R, n_p);
    
    free(avg_m1s);
    free(var_m1s);
      
#pragma mark WRITE OUTPUT FILE FOR DETECTION EFFICIENCY 
       
    for(int i=0; i<n_chromo_to_analyse; i++){
        eps_vector.push_back(eps[n_chromo_to_analyse]);
    }        
    free(nused);

#pragma mark FIND STATISTICS OF PAIRS
  
    for(int i=0; i<n_chromo_to_analyse; i++){ 
     
        index=ind_chromo_to_analyse[i];
        max_dist=chromo_length[index]-1;
        
        cout << "Starting to analyse chromosome " << chromosomes[index] << endl;
       
        avg_m2_dist=(long double *)(calloc(max_dist, sizeof(long double)));
        avg_m1a_dist=(long double *)(calloc(max_dist, sizeof(long double)));
        avg_m1b_dist=(long double *)(calloc(max_dist, sizeof(long double)));
        avg_m0_dist=(long double *)(calloc(max_dist, sizeof(long double)));
        avg_ratio_couples=(long double *)(calloc(max_dist, sizeof(long double)));
              
        var_m2_dist=(long double *)(calloc(max_dist, sizeof(long double)));
        var_m1a_dist=(long double *)(calloc(max_dist, sizeof(long double)));
        var_m1b_dist=(long double *)(calloc(max_dist, sizeof(long double)));
        var_m1_dist=(long double *)(calloc(max_dist, sizeof(long double)));
        var_m0_dist=(long double *)(calloc(max_dist, sizeof(long double)));
        var_ratio_couples=(long double *)(calloc(max_dist, sizeof(long double)));
               
        n_dist=(long int *)(calloc(max_dist, sizeof(long int)));
        n_tot_dist=(long int *)(calloc(max_dist, sizeof(long int)));
        n_chromo_dist=(long int *)(calloc(max_dist, sizeof(long int)));
        
        compute_stat_pairs_single(avg_m2_dist,avg_m1a_dist, avg_m1b_dist, avg_m0_dist,  avg_ratio_couples,  var_m2_dist, var_m1a_dist, var_m1b_dist, var_m1_dist,  var_m0_dist, var_ratio_couples,  n_dist,n_tot_dist, n_chromo_dist,   index,  start_ind,  end_ind, chromo_length[index],  m, matrix,  under,  over, chromosomes[index], max_dist);
                       
#pragma mark FIND U0 AS FUNCTION OF GENOMIC DISTANCE 
                 
        cout << "Computing u0 for chromosome " << chromosomes[index] << endl;
        compute_u0(u0, distances_u0, eps_vector[i], v0, n_p, max_dist, avg_m2_dist, avg_m1a_dist, avg_m1b_dist, n_dist, thr_u0, tol, tol_ext,chromosomes[index]);

        max_dist_u0 = distances_u0[distances_u0.size()-1];
        cout << "Maximal distance u0 for chromosome " << chromosomes[index] << "is:"<< max_dist_u0 << endl;

#pragma mark FIND THE QUANTILES AND THE THRESHOLDS AS FUNCTION OF THE GENOMIC DISTANCE  
                        
        cout << "Computing quantiles for chromosome " << chromosomes[index] << endl;
        compute_quantiles(quantiles, quant_value, u0, distances_u0, eps_vector[i], v0, n_p, avg_m2_dist, avg_m1a_dist, avg_m1b_dist, n_dist,chromosomes[index], N_draws, m);
                        
        if(quantiles.size()!=distances_u0.size()){
            cerr << "**********" << endl;
            cerr << "The sizes of distances and quantile vectors do not match!" << endl;
            cerr << "Size of distances: " << distances_u0.size() << endl;
            cerr << "Size of quantiles: " << quantiles.size() << endl;
            cerr << "Chromosome: " << chromosomes[index] << endl;
            
            exit(1);
        }
                       
#pragma mark FIND  PI_THRESHOLDS AS FUNCTION OF THE GENOMIC DISTANCE  
        
        cout << "Computing the thresholds for pi's for chromosome " << chromosomes[index] << endl;
        compute_pi_thresholds(pi_thr,  quantiles, u0,  eps_vector[i], v0,  n_p,  m,  tol);
        
        if(pi_thr.size()!=distances_u0.size()){
            cerr << "**********" << endl;
            cerr << "The sizes of distances and pi_thr vectors do not match!" << endl;
            cerr << "Size of distances: " << distances_u0.size() << endl;
            cerr << "Size of pi_thr: " << pi_thr.size() << endl;
            cerr << "Chromosome: " << chromosomes[index] << endl;
            
            exit(1);
        }  
                
        temp_string=("_"+chromosomes[index]+".out");
        strcpy(temp_file_name, file_out_pi_thr);
        strcat(temp_file_name, temp_string.c_str());
        
        cout << "Writing output file with thresholds for pi's for chromosome " << chromosomes[index] << endl;
        output_pi_thr(temp_file_name,  m,  n_p,  b,  R,  eps_vector[i],  h,  distances_u0,  pi_thr);
        cout << "Output written in " << temp_file_name << endl;
        cout << "Done!" << endl;
        cout << "******" << endl;

#pragma mark PI ESTIMATION

        cout << "\n\n\n\npi estimation\n\n\n\n" << endl;

        index_chromo=ind_chromo_to_analyse[i];
        eps_value=eps_vector[i]; 
        eps_value = eps[n_chromo_to_analyse];     
 
        cout << "Analysing chromosome " << chromosomes[index_chromo] << "(length: " << chromo_length[index_chromo] << ", det. efficiency: " <<  eps_value << ")." << endl;
        
#pragma mark OPEN THE OUTPUT FILE

        temp_file_name_2=file_pi_out;
        temp_file_name_2=temp_file_name_2+"_"+chromosomes[index_chromo]+".out";
        os_pi.open(temp_file_name_2.c_str(), ios::out);
        if(os_pi==0){
            cerr << "File " << temp_file_name_2 << " impossible to open!" << endl;
        }
 
#pragma mark DEFINE FIRST AND LAST LOCUS
        
        first_locus=start_ind[index_chromo];
        last_locus=end_ind[index_chromo];
        cout << "Analysing from " << 0 << "bp to " << b*(chromo_length[index_chromo])<< "bp" << endl;
        
#pragma mark START LOOP OVER LOCI

        for(int locus_a=first_locus; locus_a<=last_locus; locus_a++){
                
            if((last_locus-locus_a)%100==0){   
                cout << "Chr " << chromosomes[index_chromo] <<" " << locus_a-first_locus << " out of " << last_locus-first_locus << endl;  
            }
            
            for(int locus_b=first_locus; locus_b<=locus_a; locus_b++){

                if(locus_b==locus_a){
                    pi=pi_out; 
                    os_pi << pi << " " ;   
                    continue;
                }
                   
                
#pragma mark COUNT N2,N1,N0 AND COMPUTE THE RATIO

                count_n2_n1_n0(n2,n1a,n1b,n0, locus_a, locus_b, n_columns, matrix);

                if(n1a+n2<=under || n1b+n2<=under || n1a+n2>over || n1b+n2>over){
                    pi=pi_out;
                    os_pi << pi << " ";                   
                    continue;
                }
                n1=n1a+n1b; 
                ratio=(double)(n2)/((double)(n2)+(double)(n1));
                               
#pragma mark FIND U0
                
                g_dist=((locus_a-locus_b)>max_dist_u0)? max_dist_u0: (locus_a-locus_b);
                u0_value=u0[g_dist-1];
               
#pragma mark ESTIMATE THE INT. PROBABILITY
                
                estimate_pi_min_ratio(pi, diff_exp_theo, ratio, u0_value, eps_value, v0_value, n_p, tol_pi);
                                
#pragma mark WRITE RESULTS IN OUTPUT

                os_pi << pi << " ";
                   
            }
            
            os_pi << endl;

        }
  
        os_pi.close();

        temp_file_name_2=file_pi_out;
        temp_file_name_2=temp_file_name_2+"_"+chromosomes[index_chromo]+".out";
        cout << "Interaction probabilities written in " <<  temp_file_name_2 << " ." << endl;

        cout << "************" << endl;
        
        free(avg_m2_dist);
        free(avg_m1a_dist);
        free(avg_m1b_dist);
        free(avg_m0_dist);
        free(avg_ratio_couples);
        
        free(var_m2_dist);
        free(var_m1a_dist);
        free(var_m1b_dist);
        free(var_m0_dist);
        free(var_ratio_couples);
        
        free(n_dist);
        free(n_tot_dist);
        free(n_chromo_dist);
        
        u0.clear();
        distances_u0.clear();
        quantiles.clear();
        pi_thr.clear();
   
    }     

    free(eps);   

    return(0);
    
}

int main(){
    return run_slice();
}
