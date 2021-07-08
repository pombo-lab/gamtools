#include<vector>
#include<string>
#include<iostream>
#include<fstream>
#include<sstream>
#include<cmath>
#include<gsl/gsl_math.h>
#include<gsl/gsl_rng.h>
#include<gsl/gsl_randist.h>
#include<gsl/gsl_statistics_double.h>
#include<gsl/gsl_sort.h>
#include<time.h>



using namespace std;

const long double tol_ext=5e-2;
const double quant_value=0.95;
const long double pi_out=sqrt(-2);

const int under=0;
const int over=100;
const int thr_u0=100;
const long int N_draws=1e4;

const int single_chr=-1;

long double compute_v0(long double avg_m1s_m, int n_p, const double eps_fix ){
    
    double coef=(double)(1.)/(double)(n_p);
    
    return(  1 - ( 1-sqrt( pow((1-avg_m1s_m),coef) ) )/( eps_fix )  );
      
}


long double compute_eps(long double avg_m1s_m, const long double heff_R, int n_p ){
    
    double coef=(double)(1.)/(double)(n_p);
    
    return(  ((2+heff_R)/(heff_R)) * (1-sqrt(  pow((1-avg_m1s_m),coef)  ))  );
        
}


long double compute_theo_ratio(long double u0, long double eps, long double v0_value, const unsigned int np_value){
        
    const long double c0=u0;
    const long double epsno=1-eps;
    
    long double n00, nab, n2a2b, na0, n2a0, n2ab, neps00, nepsa0, neps2a0;
    long double m0,m1,m2;
        
    n00=c0*c0;
    nab=2*( (1 - 2*v0_value + c0)*c0 + (v0_value - c0)*(v0_value-c0) );
    n2a2b= (1 - 2*v0_value + c0)*(1 - 2*v0_value + c0);
    na0=2*(v0_value - c0)*c0;
    n2a0= (v0_value - c0)*(v0_value-c0);
    n2ab=2*(1 - 2*v0_value + c0)*(v0_value - c0);
    
    neps00=	n00 + epsno*(2*na0) + epsno*epsno*(2*n2a0+nab) + epsno*epsno*epsno*(2*n2ab) + epsno*epsno*epsno*epsno*n2a2b;
    nepsa0=eps*( na0 + epsno*(nab + 2*n2a0) + epsno*epsno*(3*n2ab) + 2*epsno*epsno*epsno*n2a2b);
    neps2a0=eps*eps*( n2a0 + epsno*n2ab + epsno*epsno*n2a2b);
            
    m0=pow((double)(neps00),(double)(np_value));
    m1=(2*(pow((double)(neps2a0+nepsa0+neps00),(double)(np_value)) - pow((double)(neps00),(double)(np_value))));
    m2=1-m1-m0;
            
    return ((m2/(m1+m2)));

}

long double compute_theo_ratio_pi(long double pi, long double u0, long double eps, long double v0_value, const unsigned int np_value){
        
    const long double c0=pi*v0_value+(1-pi)*u0;
    const long double epsno=1-eps;
    long double n00, nab, n2a2b, na0, n2a0, n2ab, neps00, nepsa0, neps2a0;
    long double m0,m1,m2;
        
    n00=c0*c0;
    nab=2*( (1 - 2*v0_value + c0)*c0 + (v0_value - c0)*(v0_value-c0) );
    n2a2b= (1 - 2*v0_value + c0)*(1 - 2*v0_value + c0);
    na0=2*(v0_value - c0)*c0;
    n2a0= (v0_value - c0)*(v0_value-c0);
    n2ab=2*(1 - 2*v0_value + c0)*(v0_value - c0);
    
    neps00=	n00 + epsno*(2*na0) + epsno*epsno*(2*n2a0+nab) + epsno*epsno*epsno*(2*n2ab) + epsno*epsno*epsno*epsno*n2a2b;
    nepsa0=eps*( na0 + epsno*(nab + 2*n2a0) + epsno*epsno*(3*n2ab) + 2*epsno*epsno*epsno*n2a2b);
    neps2a0=eps*eps*( n2a0 + epsno*n2ab + epsno*epsno*n2a2b);
            
    m0=pow((double)(neps00),(double)(np_value));
    m1=(2*(pow((double)(neps2a0+nepsa0+neps00),(double)(np_value)) - pow((double)(neps00),(double)(np_value))));
    m2=1-m1-m0;
            
    return ((m2/(m1+m2)));
   
}


void compute_probabilities(double * probs, long double pi, long double u0, long double eps, long double v0_value, const unsigned int np_value){
            
    const long double c0=u0;
    const long double epsno=1-eps;
    
    long double n00, nab, n2a2b, na0, n2a0, n2ab, neps00, nepsa0, neps2a0;
    long double m0,m1,m2;
        
    n00=c0*c0;
    nab=2*( (1 - 2*v0_value + c0)*c0 + (v0_value - c0)*(v0_value-c0) );
    n2a2b= (1 - 2*v0_value + c0)*(1 - 2*v0_value + c0);
    na0=2*(v0_value - c0)*c0;
    n2a0= (v0_value - c0)*(v0_value-c0);
    n2ab=2*(1 - 2*v0_value + c0)*(v0_value - c0);
    
    neps00=	n00 + epsno*(2*na0) + epsno*epsno*(2*n2a0+nab) + epsno*epsno*epsno*(2*n2ab) + epsno*epsno*epsno*epsno*n2a2b;
    nepsa0=eps*( na0 + epsno*(nab + 2*n2a0) + epsno*epsno*(3*n2ab) + 2*epsno*epsno*epsno*n2a2b);
    neps2a0=eps*eps*( n2a0 + epsno*n2ab + epsno*epsno*n2a2b);
            
    m0=pow((double)(neps00),(double)(np_value));
    m1=(2*(pow((double)(neps2a0+nepsa0+neps00),(double)(np_value)) - pow((double)(neps00),(double)(np_value))));
    m2=1-m1-m0;
    
    probs[0]=m0;
    probs[1]=m1;
    probs[2]=m2;
            
}



long double estimate_pi_ratio(long double ratio, long double u0, long double eps, long double v0_value, unsigned int n_p, long double tol){
    
    long double min=0.0, max=1.0, current=0.5, test;
    double temp_double_1, temp_double_2;
    
    temp_double_1=ratio-compute_theo_ratio_pi(min,u0, eps, v0_value,n_p);
    temp_double_2=ratio-compute_theo_ratio_pi(max, u0,eps, v0_value, n_p);
    
    test=ratio-compute_theo_ratio_pi(current, u0, eps, v0_value,n_p);
    
    while(fabs(test)>tol && current>0.01 && current<0.99){
        
        if(test<0){            
            max=current;
            current=(max+min)/2.;
            test=ratio-compute_theo_ratio_pi(current, u0, eps, v0_value, n_p );            
        }
        else{            
            min=current;
            current=(max+min)/2.;
            test=ratio-compute_theo_ratio_pi(current, u0, eps, v0_value, n_p );                        
        }                
    }

    current=(current<=0.01)? 0: current;
    current=(current>=0.99)? 1: current;
    
    return(current);

}

void estimate_pi_min_ratio(long double &pi, long double &diff_exp_theo, long double ratio, long double u0, long double eps, long double v0_value, unsigned int n_p, long double tol){
    
    long double min=0.0, max=1.0, current=0.5, test;
    double temp_double_1, temp_double_2;
    
    temp_double_1=ratio-compute_theo_ratio_pi(min,u0, eps, v0_value,n_p);
    temp_double_2=ratio-compute_theo_ratio_pi(max, u0,eps, v0_value, n_p);
    
    test=ratio-compute_theo_ratio_pi(current, u0, eps, v0_value,n_p);
    
    while(fabs(test)>tol && current>0.01 && current<0.99){
        if(test<0){            
            max=current;
            current=(max+min)/2.;
            test=ratio-compute_theo_ratio_pi(current, u0, eps, v0_value, n_p );                        
        }
        else{            
            min=current;
            current=(max+min)/2.;
            test=ratio-compute_theo_ratio_pi(current, u0, eps, v0_value, n_p );                        
        }                
    }
    
    current=(current<=0.01)? 0: current;
    current=(current>=0.99)? 1: current;
        
    pi=current;
    diff_exp_theo=test;
     
}


void count_n2_n1_n0(int &n2, int &n1a, int &n1b, int &n0, int locus_a, int locus_b, int n_tubes, vector< vector<bool> > &matrix){
    
    bool flag_a, flag_b;
    
    n2=0;
    n0=0;
    n1a=0;
    n1b=0;
    
    for(int k=0; k<n_tubes; k++){
        flag_a=matrix[locus_a][k];
        flag_b=matrix[locus_b][k];
        if(flag_a==1 && flag_b==1)  n2++;
        else if(flag_a==0 && flag_b==0) n0++;
        else if(flag_a==1 && flag_b==0) n1a++;
        else n1b++;        
    }
 
}


vector<vector<bool> > load_matrix_tubes(string file_tube, int &n_rows, int &n_columns, const int &m){
    
    ifstream is;
    bool flag;
    vector< vector<bool> > matrix;
    string temp_string;
    vector<bool> row_bool;
    const char* file_tube_char= file_tube.c_str();
    
    is.open(file_tube_char, ios::in);
    if(is.fail()){
        cerr << "File " << file_tube << " not found!" << endl;
        exit(1);
        
    }
    
    getline(is,temp_string);
    n_columns=1+temp_string.size()/2; 
    n_rows=1; 
    while(getline(is, temp_string)){
        ++n_rows;
    } 
       
    if(n_columns != m){        
        cerr << "The number of tubes calculated from the input file is not consistent with the expected value!" << endl;
        exit(1);
    }

    is.clear();
    is.seekg(0,ios::beg);
    
    for(int i=0; i<n_rows; i++){
        matrix.push_back(row_bool);        
        for(int j=0; j<n_columns; j++){
            is >> flag;
            matrix[i].push_back(flag);            
        }       
    }
    
    is.close();
  
    return(matrix);

}

void load_chromo_info(vector<string> & chromosomes,  vector<string> &chromo_to_analyse, vector<unsigned int> & start_ind, vector<unsigned int> & end_ind, vector<unsigned int> & chromo_length, vector<int> & ind_chromo_to_analyse,  int & n_chromo_to_analyse,  string file_chr_names, string file_chr_indices, int single_chr){
        
    ifstream is;
    stringstream ss;
    string chromo_current;
    int temp_int, n_chromosomes, n_loci_temp;
    bool flag;
    const char * file_chr_names_char=file_chr_names.c_str();
    const char * file_chr_indices_char=file_chr_indices.c_str();

    is.open(file_chr_names_char,ios::in); 
    if(is.fail()){
        cerr << "File with chromosomes' names not found! (' " << file_chr_names << " ')" << endl;
        exit(1);
        
    }
    while(is.eof()==0){
        is >> chromo_current;
        chromosomes.push_back(chromo_current);
    }
    is.close();
    
    is.open(file_chr_indices_char,ios::in); 
    if(is.fail()){
        cerr << "File with chromosomes' indices not found!" << endl;
        exit(1);
    }
    
    while(is.eof()==0){
        is >> temp_int;
        start_ind.push_back(temp_int);
        is >> temp_int;
        end_ind.push_back(temp_int);
        
    }
    is.close();
    
    n_chromosomes=chromosomes.size();
    for(int i=0; i<n_chromosomes; i++){
        n_loci_temp=end_ind[i]-start_ind[i]+1;
        chromo_length.push_back(n_loci_temp);
    }
    
    if(single_chr==-1){
        for(int i=1; i<20; i++){
            ss << "chr" << i;
            chromo_to_analyse.push_back(ss.str());                        
            ss.str("");             
        }
          
    }else{
        for(int i=single_chr; i<single_chr+1; i++){
            ss << "chr" << i;
            chromo_to_analyse.push_back(ss.str());
            ss.str("");        
        }
   
    }
       
    n_chromo_to_analyse=chromo_to_analyse.size();

    for(int i=0; i<n_chromo_to_analyse; i++){
        flag=0;
        for(int j=0; j<chromosomes.size(); j++){            
            if(chromosomes[j]==chromo_to_analyse[i]){                
                flag=1;
                ind_chromo_to_analyse.push_back(j);
                break;
            }          
        }
        if(flag==0){            
            cerr << "Chromosome " << chromo_to_analyse[i] << " not found!" << endl;
            exit(1);            
        }        
    }
    
}


void compute_m1s_eps(long double *avg_m1s, long double * var_m1s, long int *nused,  long double *eps, long double &avg_m1s_tot, long double &var_m1s_tot, long int &nused_tot, long int &nloci_tot, int n_chromo_to_analyse, vector<int> ind_chromo_to_analyse,  vector<unsigned int>  start_ind, vector<unsigned int> chromo_length, int m, vector<vector<bool> > &matrix, const long double heff_R, int n_p){
        
    int temp_int, m1s, index_chromo, ind_locus;
        
    avg_m1s_tot=0.;
    var_m1s_tot=0.;
    nused_tot=0;
    nloci_tot=0;
        
    for(int c=0; c<n_chromo_to_analyse; c++){
                        
        index_chromo=ind_chromo_to_analyse[c];
       
        temp_int=chromo_length[index_chromo];
        for(int locus=0; locus<temp_int; locus++){
            nloci_tot=nloci_tot+1;                        
            ind_locus=start_ind[index_chromo]+locus;
            
            m1s=0;
            for(int t=0; t<m; t++){
                if(matrix[ind_locus][t]==1){
                    m1s++;
                }
            }
            
            if(m1s!=0){
                avg_m1s[c]=avg_m1s[c]+(long double)(m1s);
                var_m1s[c]=var_m1s[c]+(long double)(m1s)*(long double)(m1s);
                nused[c]=nused[c]+1;
                
                avg_m1s_tot=avg_m1s_tot+(long double)(m1s);
                var_m1s_tot=var_m1s_tot+(long double)(m1s)*(long double)(m1s);
                nused_tot=nused_tot+1;
                
                
            }
        }
                
        avg_m1s[c]=avg_m1s[c]/nused[c];
        
        if(nused[c]>1){
            var_m1s[c]=(nused[c]/(nused[c]-1))* ((var_m1s[c]/nused[c]) - (avg_m1s[c]*avg_m1s[c]) );
        }
                                
        eps[c]=compute_eps(avg_m1s[c]/m, heff_R, n_p);        
        
    }
        
    avg_m1s_tot=avg_m1s_tot/nused_tot;
    var_m1s_tot=(nused_tot/(nused_tot-1))* ((var_m1s_tot/nused_tot) - (avg_m1s_tot*avg_m1s_tot) );
    
    eps[n_chromo_to_analyse]=compute_eps(avg_m1s_tot/m, heff_R, n_p);
        
}


void output_pi_thr(string temp_file_name, int m, int n_p, double b, double R, double eps_value, double h, vector<unsigned int> distances_u0, vector<long double> pi_thr){
        
    ofstream os;

    os.open(temp_file_name, ios::out);
    if(os.fail()){
        cerr << "Output file " << temp_file_name << " impossible to open!" << endl;
        exit(1);
    }
    
    os <<"### Number of tubes used: " << m <<endl;
    os << "### Nuclear profiles per tube: " << n_p << endl;
    os <<"### Resolution (bp): "<< b << endl;
    os <<"### Slice thickness (microns): " << h << endl;
    os <<"### Nuclear radius (microns): " << R << endl;
    os << "### Detection efficiency: " << eps_value << endl;
    
    
    os << "distance\t pi_threshold" << endl;

    for(int i=0; i<distances_u0.size(); i++){
               
        os << distances_u0[i] << "\t" << pi_thr[i] <<  endl;

    }
    
    os.close();

}


void compute_stat_pairs_single(long double *avg_m2_dist,long double *avg_m1a_dist,long double *avg_m1b_dist,long double *avg_m0_dist, long double *avg_ratio_couples, long double *var_m2_dist,long double *var_m1a_dist,long double *var_m1b_dist,long double *var_m1_dist, long double *var_m0_dist, long double*var_ratio_couples, long int *n_dist,long int *n_tot_dist,long int *n_chromo_dist,  int ind_chromo_to_analyse,  vector<unsigned int>  start_ind, vector<unsigned int>  end_ind, unsigned int chromo_length, int m, vector<vector<bool> > &matrix, int under, int over, string chromosome,long int max_dist){
 
    int locus1, locus2, n0, n1a, n1b, n2, index, index_chromo, temp_int;
    bool flag1, flag2;

    for(int d=1; d<=max_dist; d++){
        
        if(d%100==0){
            cout << "Chromosome " << chromosome <<  ", calculating statistics of pairs for distance: " << d << "/" << max_dist << endl;
        }
       
        index_chromo=ind_chromo_to_analyse;
        
        if(end_ind[index_chromo]-start_ind[index_chromo]<d){
            cerr << "Problem with chromosome length!!" << endl;
            exit(1);
        }
            
        n_chromo_dist[d-1]=n_chromo_dist[d-1]+1;
                        
        temp_int=chromo_length-d;
        for(int pair=0; pair<temp_int; pair++){
            n_tot_dist[d-1]=n_tot_dist[d-1]+1;
                
            locus1=start_ind[index_chromo]+pair;
            locus2=locus1+d;
                
            n1a=0;
            n1b=0;
            n2=0;
            n0=0;
                
            for(int t=0; t<m; t++){
                flag1=matrix[locus1][t];
                flag2=matrix[locus2][t];
                    
                if(flag1==1 && flag2==1) n2++;
                else if(flag1==0 && flag2==0) n0++;
                else if(flag1==1 && flag2==0) n1a++;
                else n1b++;
                    
            }
                
            if(n2+n1a+n1b+n0!=m){
                    
                cerr << "Problem distance " << d << "chromosome " << chromosome << "loci " << locus1 << " " << locus2 << endl;
                exit(1);
                    
            }
                
            if(n1a+n2<=under || n1b+n2<=under || n1a+n2>over || n1b+n2>over){
                    
                continue;
                    
            }
                
            else{
                    
                n_dist[d-1]=n_dist[d-1]+1;
                avg_m2_dist[d-1]=avg_m2_dist[d-1]+n2;
                avg_m1a_dist[d-1]=avg_m1a_dist[d-1]+n1a;
                avg_m1b_dist[d-1]=avg_m1b_dist[d-1]+n1b;
                avg_m0_dist[d-1]=avg_m0_dist[d-1]+n0;
                avg_ratio_couples[d-1]=avg_ratio_couples[d-1]+((long double)(n2)/(long double)(n1a+n1b+n2));
                    
                var_m2_dist[d-1]=var_m2_dist[d-1]+(n2*n2);
                var_m1a_dist[d-1]=var_m1a_dist[d-1]+(n1a*n1a);
                var_m1b_dist[d-1]=var_m1b_dist[d-1]+(n1b*n1b);
                var_m1_dist[d-1]=var_m1_dist[d-1]+((n1a+n1b)*(n1a+n1b));
                var_m0_dist[d-1]=var_m0_dist[d-1]+(n0*n0);
                var_ratio_couples[d-1]=var_ratio_couples[d-1]+ ((long double)(n2)/(long double)(n1a+n1b+n2)) * ((long double)(n2)/(long double)(n1a+n1b+n2));
            }
    
        }

        avg_m2_dist[d-1]=avg_m2_dist[d-1]/n_dist[d-1];
        avg_m1a_dist[d-1]=avg_m1a_dist[d-1]/n_dist[d-1];
        avg_m1b_dist[d-1]=avg_m1b_dist[d-1]/n_dist[d-1];
        avg_m0_dist[d-1]=avg_m0_dist[d-1]/n_dist[d-1];
        avg_ratio_couples[d-1]=avg_ratio_couples[d-1]/n_dist[d-1];
                
        if(n_dist[d-1]>1){
            var_m2_dist[d-1]=(n_dist[d-1]/(n_dist[d-1]-1))*( (var_m2_dist[d-1]/n_dist[d-1]) - (avg_m2_dist[d-1]*avg_m2_dist[d-1])  );
            var_m1a_dist[d-1]=(n_dist[d-1]/(n_dist[d-1]-1))*( (var_m1a_dist[d-1]/n_dist[d-1]) - (avg_m1a_dist[d-1]*avg_m1a_dist[d-1])  );
            var_m1b_dist[d-1]=(n_dist[d-1]/(n_dist[d-1]-1))*( (var_m1b_dist[d-1]/n_dist[d-1]) - (avg_m1b_dist[d-1]*avg_m1b_dist[d-1])  );
            var_m1_dist[d-1]=(n_dist[d-1]/(n_dist[d-1]-1))*( (var_m1_dist[d-1]/n_dist[d-1]) - (  (avg_m1a_dist[d-1]+avg_m1b_dist[d-1])*(avg_m1a_dist[d-1]+avg_m1b_dist[d-1]) ));
            var_m0_dist[d-1]=(n_dist[d-1]/(n_dist[d-1]-1))*( (var_m0_dist[d-1]/n_dist[d-1]) - (avg_m0_dist[d-1]*avg_m0_dist[d-1])  );
            var_ratio_couples[d-1]=(n_dist[d-1]/(n_dist[d-1]-1))*( (var_ratio_couples[d-1]/n_dist[d-1]) - (avg_ratio_couples[d-1]*avg_ratio_couples[d-1])  );
        }
        else{
        var_m2_dist[d-1]=0;
        var_m1a_dist[d-1]=0;
        var_m1b_dist[d-1]=0;
        var_m1_dist[d-1]=0;
        var_m0_dist[d-1]=0;
        var_ratio_couples[d-1]=0;                        
        }
    }
    
}

void compute_u0(vector<long double> &u0, vector<unsigned int> &distances_u0,long double eps, const long double v0_value, const unsigned int n_p, long double max_dist, long double *avg_m2, long double *avg_m1a, long double *avg_m1b, long int * n_dist, const long double thr_u0, const long double tol, const long double tol_ext, string chromosome){
    
    long double initial_max=v0_value, initial_min=v0_value*v0_value;
    long double min, max;
    long double m2_value, m1_value, ratio_value, current_best, current_diff;
    long double temp_double_1, temp_double_2;

    for(int d=1; d<=max_dist; d++){
        
        if(d%100==0){
            cout << "Chromosome " << chromosome <<  ", calculating u0 for distance: " << d << "/" << max_dist << endl;
        }

        if(n_dist[d-1]<thr_u0){
            continue;
        }
        
        m2_value=avg_m2[d-1];
        m1_value=avg_m1a[d-1]+avg_m1b[d-1];
        
        ratio_value=(m2_value)/(m1_value+m2_value);
        
        min=initial_min;
        max=initial_max;

        temp_double_1=ratio_value-compute_theo_ratio(min, eps, v0_value,n_p);
        temp_double_2=ratio_value-compute_theo_ratio(max, eps, v0_value, n_p);
       
        if(temp_double_1*temp_double_2>0 & temp_double_2>tol_ext){
        
            cerr << "**********" << endl;
            cerr << "Problem estimating u0 for chromosome " << chromosome << ", distance " << d << endl;
            cerr << "Min&max u0: " << min << " " << max << endl;
            cerr << "Ratio computed: " << ratio_value << endl;
            cerr << "Min&Max expected ratio: " << compute_theo_ratio(min, eps, v0_value,n_p) << " " << compute_theo_ratio(max, eps, v0_value, n_p) << endl; 
            cerr << "np: " << n_p << endl;           
            cerr << "Diff. between expected and theo at min&max: " << temp_double_1 << " " << temp_double_2 << endl;
            cerr << "*********" << endl;
            exit(1);
        }
       
        else if(temp_double_1*temp_double_2<0){
        
            current_best=(min+max)/2;
            current_diff=ratio_value-compute_theo_ratio(current_best, eps, v0_value, n_p);
        
            while(fabs(current_diff)>tol){
            
                if(current_diff<0){
                
                    max=current_best;
                    current_best=(max+min)/2.;
                    current_diff=ratio_value-compute_theo_ratio(current_best, eps, v0_value,n_p);
                
                }
                else{
                    min=current_best;
                    current_best=(max+min)/2.;
                    current_diff=ratio_value-compute_theo_ratio(current_best, eps, v0_value, n_p);
                }
            }
            u0.push_back(current_best);
            distances_u0.push_back(d);
        }
        
        else if(temp_double_1*temp_double_2>0 & temp_double_2>0 & temp_double_2<tol_ext){
            current_best=max;
            u0.push_back(current_best);
            distances_u0.push_back(d);
        
        }
        else if(temp_double_1*temp_double_2>0 & temp_double_2<0 & (-1)*(temp_double_1)<tol_ext){
            current_best=min;
            u0.push_back(current_best);
            distances_u0.push_back(d);
        }
        
        else if(temp_double_1*temp_double_2>0 & temp_double_2<0 & (-1)*(temp_double_1)>tol_ext){
            cerr << "**********" << endl;
            cerr << "Problem estimating u0 for chromosome " << chromosome << ", distance " << d << endl;
            cerr << "Min&max " << min << " " << max << endl;
            cerr << "Diff. between expected and theo at min&max: " << temp_double_1 << " " << temp_double_2 << endl;
            cerr << "*********" << endl;
            
            exit(1);
            
        }
        
    }
    
}

void compute_quantiles(vector<long double> & quantiles, double quant_value, vector<long double> &u0, vector<unsigned int> &distances_u0,long double eps, const long double v0_value, const unsigned int n_p,  long double *avg_m2, long double *avg_m1a, long double *avg_m1b, long int * n_dist, string chromosome, long int N_draws, int m){
    
    double probs[3];
    size_t cat=3;
    double ratios[N_draws];
    unsigned int result[3];
    long int max_dist=u0.size();

    gsl_rng *r;
    const gsl_rng_type * T;
   
    gsl_rng_env_setup();
    T = gsl_rng_default;
    r = gsl_rng_alloc (T);

    gsl_rng_set (r, time(NULL));

    for(int d=1; d<=max_dist; d++){
  
        if(d%100==0){
            cout << "Chromosome " << chromosome <<  ", computing the quantiles for distance: " << d << "/" << max_dist << endl;
        }
        
        compute_probabilities(probs, 0.0, u0[d-1], eps, v0_value, n_p);

        for(int i=0; i<N_draws; i++){
            
            gsl_ran_multinomial(r, cat, m, probs, result);
            ratios[i]=((double)(result[2])/((double)(result[2])+(double)(result[1])));
        
        }
        
        gsl_sort(ratios, 1, N_draws);
        
        quantiles.push_back((double)(gsl_stats_quantile_from_sorted_data(ratios, 1, N_draws, quant_value)));
        
    }

}


void compute_pi_thresholds(vector<long double> & pi_thr, vector<long double>  quantiles, vector<long double> &u0, long double eps, const long double v0_value, const unsigned int n_p, int m, long double tol){
    
    int max_dist=quantiles.size();
    
    for(int i=0; i<max_dist; i++){
        
        pi_thr.push_back(estimate_pi_ratio(quantiles[i], u0[i], eps, v0_value, n_p, tol));
        
    }
    
}

long double compute_heff_R(long double h, long double R, long double bL){

    return ((long double)(h)/(long double)(R))+2*pow((long double)(bL),(long double)(1./3.));

}

int run_slice(string file_tube, string file_pi_out, string file_out_pi_thr, string file_chr_names, string file_chr_indices, int m, long L, int b, double h, double R, int n_p){

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

    string temp_file_name;

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
    const long double heff_R=compute_heff_R(h, R, bL);
    const long double v0=2./(2.+(heff_R));
 
    const double v0_value=2./(2.+(heff_R));

    cout << "heff/R: " << heff_R << endl;

#pragma mark LOAD THE MATRIX OF TUBES 
    
    cout << "Loading matrix: " << file_tube << endl;    
    matrix=load_matrix_tubes(file_tube, n_rows, n_columns, m);    
    cout << "Rows (loci) and columns (tubes) counted (" << n_rows << " " << n_columns << ")" << endl;
    cout << "Matrix loaded!" << endl;
    cout << "************" << endl;
    
#pragma mark LOAD CHROMOSOMES INFO 
    
    cout << "Loading chromosomes info..." << endl;
    load_chromo_info(chromosomes, chromo_to_analyse, start_ind, end_ind, chromo_length, ind_chromo_to_analyse, n_chromo_to_analyse, file_chr_names, file_chr_indices, single_chr);
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
                
        temp_file_name=file_out_pi_thr;
        temp_file_name=temp_file_name+"_"+chromosomes[index]+".out";
        
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
        if(os_pi.fail()){
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

int compilation_test(){
    return(4);
}

