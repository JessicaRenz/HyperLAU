#include <cmath>
#include <cstdlib>
#include <cstring>
#include <vector>
#include <armadillo>
#include <stdio.h>
#include <numeric>
#include <fstream>
#include <unordered_map>
#include <chrono>

using namespace std;
using namespace arma;
using namespace std::chrono;
using std::ofstream;
using Matrix = std::vector<std::vector<double>>;

int L;
string str;
vector<string> R;

/* (just taken from HyperHMM)
A function for converting a number into a binary string.
Input variables:
- int n: The number you want to convert from int to binary string
- int L: The total length of the binary string. This number needs to be greater then log_2(n)
Output:
- A string of length L representing the integer n as a binary string
*/
string number2binary(int n, int L){
	string binary = "";

	//Create the binary number (without added 0's to the left)
	for (int v = pow(2,L-1); v>=1; v /=2){
		if (n >= v){
			n -= v;
			binary = binary + "1";
		}else{
			binary = binary + "0";
		}
	}
	//Add zeros to the left such that it has length L
	while (binary.length() < L){
		binary = "0" + binary;
	}
	return binary;
}

/*(just taken from HyperHMM)
A function for converting a binary string into a number.
Input variables:
- string bin: The binary string you want to convert to a number
- int L: The total length of the binary string.
Output:
- An integer form of the binary string
*/
int binary2int(string bin, int L){
	int number = 0;

	//Loop through the binary number, checking at each position if the entry is a 1, if yes, add the corresponding value to the output
	for (int i = 0; i<L;i++){
		if (bin[i] == '1'){
			number +=pow(2,L-i-1);
		}
	}
	return number;
}

//just taken from HyperHMM
/*
A function for creating numerous lists with information about the possible states to go to from a given state.
Input variables (All of the vectors needs to be empty as input):
- vector<int> n_partners: The number of possible states to go to from the state represented by the index.
- vector<int> cumulative_partners: A vector where the element at index i represents the index where t=i starts in the vector possible
- vector<string> partners: A vector containing all the indices (in CRS format) for states in correct order.
- int L: The total length of the binary string.
Output:
- Updated versions of n_partners, cumulative_partners, and partners.
*/
void possible_transitions(vector<int>& n_partners, vector<int>& cumulative_partners, vector<int>& partners, int L){
	int index = 0;
	int end_node_int;
	string end_node;
	//Loop through all the states
	for (int i = 0; i < pow(2,L); i++){
		string vertex = number2binary(i,L);
		n_partners.push_back(0);
		for (int j = 0; j < L; j ++){
			if (vertex[j] == '0'){
				//if we find a 0 in the state it is possible to add a 1. Hence it is possible to go to the state with all equal elements except this place.
				end_node = vertex;
				end_node[j] = '1';
				end_node_int = binary2int(end_node, L);
				partners.push_back(end_node_int);
				n_partners[i]++;
			}
		}
		//Add the correct number to the cumulative vector
		if ( i == 0){
			cumulative_partners.push_back(0);
		}else{
			int c = cumulative_partners[i-1] + n_partners[i -1];
			cumulative_partners.push_back(c);
		}
	}
}

// a function that removes duplications in the input data and counts how often they occured in the original dataset
vector<string> red_data(vector<string> original_data, vector<int>& frequencies) {
    	vector<string> unique_data;
    	unordered_map<string, int> frequency;

    	for (const auto& str : original_data) {
        	if (frequency.find(str) == frequency.end()) {
            	frequency[str] = 1; // First occurrence of the number
            	unique_data.push_back(str);
        	} else {
            	++frequency[str]; // Increment the frequency
        	}
    	}
    
    	for (string str : unique_data) {
        	frequencies.push_back(frequency[str]);
    	}

    	return unique_data;
}


//creates a uniform rate matrix, as an initial input, with a base rate of 1 for every feature and 0 for all other pairwise interactions
void uniform_rate_matrix(vector<double>& rate_matrix, int model, int L){

	if (model == 1){
		for (int i = 0; i < L; i++){
			rate_matrix[i] = 1;
		}
	}else if (model == 2){
		for (int i = 0; i < L; i++){
			rate_matrix[i*L] = 1;
		}
	}else if (model == 3){
		for (int i = 0; i < L; i++){
			rate_matrix[i*L*L] = 1;
		}
	}else if (model == 4){
		for (int i = 0; i < L; i++){
			rate_matrix[i*L*L*L] = 1;
		}
	}else if (model == -1){
		for (int i = 0; i < L; i++){
			rate_matrix[i*pow(2,L)] = 1;
		}
	}else{
		printf("Invalid value for model \n");
      		exit(1);
	}

}

//just taken from HyperHMM

/* A  function that calculates the uniformized values >0 for the transition matrix, with the corresponding information about their position in the matrix.
Input variables: (All of the vectors need to be empty as input):
- vector A_val: contains the values (>0) that occur in the uniform transition matrix
- vector A_row_ptr: stores the information about the corresponding row in the transition matrix
- vector A_col_idx: stores the information about the corresponding column in the transition matrix
- int L: the total length of the binary string
Output:
- Updated versions of A_val, A_row_ptr and A_col_idx
*/
void uniform_transition_matrix(arma::vec& A_val, arma::vec& A_row_ptr, arma::vec& A_col_idx, int L){
	vector<int> n_partners;
	vector<int> c_partners;
	vector<int> partners;
	possible_transitions(n_partners, c_partners, partners, L);
	int k = 0;
	int c = 0;
	//Loop through all the transition probabilities (>0) and update the specific values to be uniform
	for (int i = 0; i<pow(2,L); i++){
		int n_end_vertices = n_partners[i];
		int r = A_row_ptr(i);
		A_row_ptr(i+1) = n_end_vertices +r;
		c = 0;
		for (int j=0; j<n_end_vertices; j++){
			int j2 = partners[k];
			A_val(r+c) = 1./n_end_vertices;
			A_col_idx(r+c) = j2;
			k++;
			c++;
		}
	}
}

//calculates the rate for an transition from j to i
double rate(int i, int j, vector<double> x_current, int L, int model){

	double rate = 0;
	string i_str = number2binary(i,L);
	string j_str = number2binary(j,L);
	int n_diff = 0;
	for (int k = 0; k < L; k++){
		if (i_str[k] != j_str[k]){
			n_diff = n_diff + 1;
		}
	}

	int locus;
	if (n_diff != 1){
		
		rate = 0;
	}else {
		for (int k = 0; k< L; k++){
			if (i_str[k] != j_str[k]){
				locus = k;
	
				if (model == 1){
					rate = x_current[locus];
				}else if (model == 2){
					rate = x_current[locus*L + locus];
					for (int l= 0; l<L; l++){
						int digit = j_str[l]- '0';
						rate = rate + digit*x_current[l*L + locus];
					}
				}else if (model == 3){
					rate = x_current[locus*L*L+ locus*L + locus];
					for (int l = 0; l<L; l++){
						for (int m = l; m<L; m++){
							int digit_l = j_str[l] - '0';
							int digit_m = j_str[m] - '0';
							rate = rate + digit_l*digit_m*x_current[l*L*L + m*L + locus];
						}
					}
				}else if (model == 4){
					rate = x_current[locus*L*L*L+locus*L*L + locus*L + locus];
					for (int l = 0; l<L; l++){
						for (int m = l; m<L;m++){
							for (int n = m; n<L; n++){
								int digit_l = j_str[l] - '0';
								int digit_m = j_str[m] - '0';
								int digit_n = j_str[m] - '0';
								rate = rate + digit_l*digit_m*digit_n*x_current[l*L*L*L + m*L*L + n*L + locus];
							}
						}
					}
				}else{
					exit(1);
				}
			}
		}

		rate =  exp(rate);
		
	}
	
	return rate;
}

//transforms an vector with transition rates into the form of a matrix including transition probabilities
void build_trans_matrix(int L, vector<double> x_current, int model, mat& mat_temp){
	for (int j = 0; j< pow(2,L); j++){
		for (int i = 0; i < pow(2,L); i++){
			if (j< i){
				mat_temp(i,j) = rate(i,j,x_current,L,model);
				if (mat_temp(i,j) < 0){
					mat_temp(i,j) = 0;
				}
			}
		}	
		double sum = 0;
		for (int i = 0; i<pow(2,L); i++){
			sum = sum + mat_temp(i,j);
		}
		for (int i = 0; i<pow(2,L); i++){
			mat_temp(i,j) = mat_temp(i,j)/sum;
		}
	}
	for (int k = 0; k < pow(2,L); k++){
		mat_temp(k,pow(2,L)-1) = 0;
	}
}

// calculates the log-likelihood to see a the data given a rate vector
double loglh(vector<string> before, vector<string> after, vector<int> frequ, int L, vector<double> x_current, int model){

	vec loglh(before.size(), fill::zeros);
	
	// 2.) Use the transition matrix state vector approach to simulate the system over timesteps t
        mat P_desc(pow(2,L),L+1, fill::zeros);
	P_desc(0,0) = 1;
	mat mat_temp(pow(2,L),pow(2,L), fill::zeros);
	build_trans_matrix(L,x_current,model,mat_temp);
	for (int t = 1; t <= L; t++){
		P_desc.col(t) = mat_temp * P_desc.col(t-1);
	}

	// 1.) Write down the sets of all the compatible states for the before and the after datapoint
		for (int l = 0; l< before.size(); l++){


			string DP_b = before[l];

			vec c_before(pow(2,L), fill::zeros);

			for (int i = 0; i<pow(2,L); i++){
				for (int k = 0; k<L;){
					string node = R[i];
					if((DP_b[k] == '?')&&(k<L-1)){
						k++;
					}else if((DP_b[k] == '?')&&(k=L-1)){
						c_before(i) = 1;
						break;
					}else if((DP_b[k] == node[k])&&(k<L-1)){
						k++;
					}else if((DP_b[k]==node[k])&&(k=L-1)){
						c_before(i) = 1;
						break;
					}else{
						break;
					}

				}
			}

			string DP_a = after[l];

			vec c_after(pow(2,L), fill::zeros);

			for (int i = 0; i<pow(2,L); i++){
				for (int k = 0; k<L;){
					string node = R[i];
					if((DP_a[k] == '?')&&(k<L-1)){
						k++;
					}else if((DP_a[k] == '?')&&(k=L-1)){
						c_after(i) = 1;
						break;
					}else if((DP_a[k] == node[k])&&(k<L-1)){
						k++;
					}else if((DP_a[k]==node[k])&&(k=L-1)){
						c_after(i) = 1;
						break;
					}else{
						break;
					}
				}
			}

			
                        mat P = P_desc;

			// 3.) Loop through t from 0 to L

			vector<double> comp_probs;
			mat Q(pow(2,L),L+1, fill::zeros);

			for (int t = 0; t<=L; t++){
				vec P_0 = P.col(t);

				vec P_0_comp(pow(2,L), fill::zeros);
				for(int i = 0; i< c_before.size(); i++){
					P_0_comp[i] = c_before[i]* P_0[i];
				}

				if (sum(P_0_comp) != 0){
					vec Q_0_comp(pow(2,L), fill::zeros);
					vec Q_col_0 = P_0_comp;

					mat Q(pow(2,L),L+1, fill::zeros);
					Q.col(t) = Q_col_0;
					for (int s = t+1; s<=L; s++){
						Q.col(s) = mat_temp * Q.col(s-1);
					}

					for(int s = 0; s<=L; s++){
						vec Q_0 = Q.col(s);
						for(int i = 0; i< c_after.size(); i++){
							Q_0_comp[i] = c_after[i]*Q_0[i];
						}

						double sum_0 = sum(Q_0_comp);
						if (sum_0 != 0){
							comp_probs.push_back(sum_0);
						}

					}
				}
			}

			double lh = 0;
			for (int i = 0; i< comp_probs.size(); i++){
				lh = lh + comp_probs[i];
			}

			double factor = (L+1)*(L+2);
			double normlh = 2/factor*lh;
			loglh[l] = log(normlh)*frequ[l];

		}

		return sum(loglh);

}


// calculates the log-likelihood to see a the data given a transition matrix 
double loglh_minus_1(vector<string> before, vector<string> after, vector<int> frequ, int L, mat x_current){

//Calculate the new likelihood given the new matrix
		vec loglh(before.size(), fill::zeros);
		
	        // 2.) Use the transition matrix state vector approach to simulate the system over timesteps t
                mat P_desc(pow(2,L),L+1, fill::zeros);
		P_desc(0,0) = 1;
		for (int t = 1; t<=L; t++){
			P_desc.col(t) = x_current * P_desc.col(t-1);
		}

		// 1.) Write down the sets of all the compatible states for the before and the after datapoint
		for (int l = 0; l< before.size(); l++){


			string DP_b = before[l];

			vec c_before(pow(2,L), fill::zeros);

			for (int i = 0; i<pow(2,L); i++){
				for (int k = 0; k<L;){
					string node = R[i];
					if((DP_b[k] == '?')&&(k<L-1)){
						k++;
					}else if((DP_b[k] == '?')&&(k=L-1)){
						c_before(i) = 1;
						break;
					}else if((DP_b[k] == node[k])&&(k<L-1)){
						k++;
					}else if((DP_b[k]==node[k])&&(k=L-1)){
						c_before(i) = 1;
						break;
					}else{
						break;
					}

				}
			}


			string DP_a = after[l];

			vec c_after(pow(2,L), fill::zeros);

			for (int i = 0; i<pow(2,L); i++){
				for (int k = 0; k<L;){
					string node = R[i];
					if((DP_a[k] == '?')&&(k<L-1)){
						k++;
					}else if((DP_a[k] == '?')&&(k=L-1)){
						c_after(i) = 1;
						break;
					}else if((DP_a[k] == node[k])&&(k<L-1)){
						k++;
					}else if((DP_a[k]==node[k])&&(k=L-1)){
						c_after(i) = 1;
						break;
					}else{
						break;
					}
				}
			}


			mat P = P_desc;


			// 3.) Loop through t from 0 to L

			vector<double> comp_probs;
			mat Q(pow(2,L),L+1, fill::zeros);

			for (int t = 0; t<=L; t++){
				vec P_0 = P.col(t);

				vec P_0_comp(pow(2,L), fill::zeros);
				for(int i = 0; i< c_before.size(); i++){
					P_0_comp[i] = c_before[i]* P_0[i];
				}

				if (sum(P_0_comp) != 0){
					vec Q_0_comp(pow(2,L), fill::zeros);
					vec Q_col_0 = P_0_comp;

					mat Q(pow(2,L),L+1, fill::zeros);
					Q.col(t) = Q_col_0;
					for (int s = t+1; s<=L; s++){
					Q.col(s) = x_current * Q.col(s-1);
					}

					for(int s = 0; s<=L; s++){
						vec Q_0 = Q.col(s);
						for(int i = 0; i< c_after.size(); i++){
							Q_0_comp[i] = c_after[i]*Q_0[i];
						}

						double sum_0 = sum(Q_0_comp);
						if (sum_0 != 0){
							comp_probs.push_back(sum_0);
						}

					}
				}
			}

			double lh = 0;
			for (int i = 0; i< comp_probs.size(); i++){
				lh = lh + comp_probs[i];
			}

			double factor = (L+1)*(L+2);
			double normlh = 2/factor*lh;
			loglh[l] = log(normlh)*frequ[l];
			if(normlh <=0){
				printf("Formal error in dataset\n");
      				exit(1);
			}

		}
		return sum(loglh);


}

//performs the simulated annealing for model 1 - 4 (based on a rate vector)
void simulated_annealing(vector<double> x_initial, vector<double>& best_mat, int L, vector<string> before, vector<string> after, vector<int> frequ, int model, vector<double>& prog_best_lik, double denom){
	double temp = 1;
	vector<double> x_old = x_initial;
	vector<double> x_current = x_initial;
	vector<double> x_current_temp = x_initial;
	
	double lik_initial = loglh(before, after, frequ, L, x_initial, model);
	
	double best_lik = lik_initial;
  	double new_lik = lik_initial;
  	double old_lik = lik_initial;
  	double thresh = 0.000001;
	double up = log(thresh/temp);
	double down = log(1/denom);
	double num_it = up/down;
  	while (temp > thresh){
  		
  		
  		auto start = high_resolution_clock::now();
  		
  		
  		//Small perturbations to get a new matrix
  		for (int i = 0; i < x_old.size(); i++){
  			x_current[i] = x_old[i] + (drand48() - 0.5)*0.05;
  		}

  		x_current_temp = x_current;

  		new_lik = loglh(before, after, frequ, L, x_current, model);

  		if (new_lik > old_lik || exp(-(old_lik -new_lik)/temp) > drand48()){
			old_lik = new_lik;
	 		x_old = x_current;
		}

		if (new_lik > best_lik){
			best_lik = new_lik;
			best_mat = x_current;
		}

		prog_best_lik.push_back(best_lik);
		
		if (temp == 1){
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<milliseconds>(stop - start);
			cout << "Estimated time: " << num_it*duration.count()/60000 << " minutes" << endl;
			
		}

		temp = temp/denom;
  	}
}

//performs simulated annealing for model -1 (based on a transition matrix)
void simulated_annealing_minus_1(mat x_initial, mat& best_mat, int L, vector<string> before, vector<string> after, vector<int> frequ, int model, vector<double>& prog_best_lik, double denom){
	double temp = 1;
	mat x_old = x_initial;
	mat x_current = x_initial;
	mat x_current_temp = x_initial;

	double lik_initial = loglh_minus_1(before, after, frequ, L, x_initial);

	double best_lik = lik_initial;
  	double new_lik = lik_initial;
  	double old_lik = lik_initial;

	double thresh = 0.000001;
	double up = log(thresh/temp);
	double down = log(1/denom);
	double num_it = up/down;
	
  	while (temp > thresh){
  	
  		auto start = high_resolution_clock::now();
  		
  		//Small perturbations to get a new matrix
  		for(int i = 0; i < pow(2,L); i++){
  			for( int j = 0; j < pow(2,L); j++){
  				if(x_old(i,j) == 0){
  					x_current(i,j) = x_old(i,j);
  				}else{
  					x_current(i,j) = x_old(i,j) + (drand48()-0.5)*0.05;
  					if( x_current(i,j) < 0){
  						x_current(i,j) = 0.0001;
  					}
  				}
  			}
		}


  		x_current_temp = x_current;

  		for(int i = 0; i < pow(2,L); i++){
			for( int j = 0; j < pow(2,L)-1; j++){
				x_current(i,j) /= sum(x_current_temp.col(j));
			}
		}


  		new_lik = loglh_minus_1(before, after, frequ, L, x_current);

  		if (new_lik > old_lik || exp(-(old_lik -new_lik)/temp) > drand48()){
			old_lik = new_lik;
	 		x_old = x_current;
		}

		if (new_lik > best_lik){
			best_lik = new_lik;
			best_mat = x_current;
		}

		prog_best_lik.push_back(best_lik);
		
		if (temp == 1){
			auto stop = high_resolution_clock::now();
			auto duration = duration_cast<milliseconds>(stop - start);
			cout << "Estimated time: " << (num_it*duration.count())/60000 << " minutes" <<  endl;
			
		}

		temp = temp/denom;
  	}
}


int main(int argc, char**argv){
// arguments: [name of the input file] [name of the output file] [number of bootstraps] [number for a random seed] [model] [denom]
        const int expectedArgs = 7;
        
        if(argc != expectedArgs){
                std::cerr << "Non-valid number of arguments.\nUsage: " << argv[0] << "[name of input file] [label for output files] [number of bootstrap resamples] [random seed] [model structure] [annealing rate]\n\n";

		printf("[name of input file]\n\tName of the file that contains the input data, including possible extensions like .txt. HyperLAU expects as an input a textfile containing a list of ancestor and descendant states separated by a blank space, for example 01? 011. Both states are encoded by binary strings, but can contain one or more ? to mark missing or uncertain data. Every line is considered as a sample independent of the others. For using cross-sectional data, just set all ancestor states to the zero-string.\n\n[label for the output files]\n\tHyperLAU will output several text-files. With this input parameter, you can specify the basis of the names of all these outputs. For every run, you will get the two output files best_likelihood_[name of the output file].txt and transitions_[name of the output file].txt. If the number of bootstrap resamples is specified as > 0 , you will also get two additional files called mean_[name of the output file].txtand sd_[name of the output file].txt.\n\n[number of bootstrap resamples]\n\tNumber of resamples to simulated. If specified as 0, no bootstrapping will be done.\n\n[random seed]\n\tSpecifies the integer random seed to use for simulations.\n\n[model]\n\tHas to be one of the following integers: -1 , 1 , 2 , 3 , 4 . With this input parameter you choose which model, i.e. what degree of allowed interaction between the features should be used. Model -1 corresponds to the model of arbitrary dependencies, where all combinations of features can influence each other. In the article introducing HyperLAU it is labeled as F for full. In model 1, every feature occurs with a fixed rate, independent of other features already obtained. In model 2 there are already pairwise interactions allowed, and in model 3 and 4, pairs or triples can influence the probability of features to occur next, respectively.\n\n[annealing rate]\n\tThis parameter specifies how fast the temperature in the Simulated Annealing Process should be decreased. After every optimization loop, the current temperature is devided by this parameter: temp = temp/denom. This parameter should be a double > 1 (e.g. 1.001).\n\n");
                return 1;
        }
	
	string file_name = argv[1];
	string out_name = argv[2];
	int bootstrap = atoi(argv[3]);
	int seed = atoi(argv[4]);
	srand48(seed);
	int model = atoi(argv[5]);
	double denom = atof(argv[6]);
	
	
	
	printf("Reading in data \n");
	//read in the data
	ifstream in(file_name);
    	vector<string> data;
    	while(getline(in, str)){
        	if(str.size()>0){
            		data.push_back(str);
        	}
    	}
    	in.close();
    	
    	
    	if(data.empty()){
    		printf("Formal error in dataset\n");
      		exit(1);
    	}
    	
    	vector<int> frequ;
    	vector<string> reduced_data = red_data(data,frequ);
    	
    	
    	string first_row = reduced_data[0];
    	L = first_row.find(' ');
    	cout << "L: " << L << endl;
    		
	for (int i = 0; i<pow(2,L); i++){
		string x = number2binary(i,L);
		R.push_back(x);
  	}
  	
  	
    	string S;
	vector<string> all;

    	for(int i =0; i < reduced_data.size(); i++){
    		S = reduced_data[i];
    		stringstream ss(S);
    		string word;
        	while (ss >> word) { // Extract word from the stream.
        		all.push_back(word);
    		}
  	}

  	string next;
  	vector<string> before;
  	vector<string> after;
  	for(int j = 0; j < all.size(); j++){
  		next = all[j];
		if (j % 2 == 0){
			before.push_back(next);
		} else{
			after.push_back(next);
		}
  	}
	
	cout << "Inference started with dataset " << file_name << " under model " << model << " with " << bootstrap << " bootstrap resamples" << endl;

	//initialising the rate matrix
	int num_param = pow(L,model);
	vector<double> rate_vec(num_param, 0.0);
	arma::mat rate_matrix(pow(2,L),pow(2,L), arma::fill::zeros);
	if (model != -1){
		uniform_rate_matrix(rate_vec, model, L);
	}else{
		arma::vec A_val(pow(2,L-1)*L, arma::fill::zeros);
		arma::vec A_row_ptr(pow(2,L)+1, arma::fill::zeros);
		arma::vec A_col_idx(pow(2,L-1)*L, arma::fill::zeros);

		//Calculation of the values in the transition matrix and corresponding informations about their placement
		uniform_transition_matrix(A_val, A_row_ptr, A_col_idx, L);
		

		//Fitting these values into the shape of a matrix
		int c = 0;
		int j = 1;
		for (int i = 0; i < pow(2,L-1)*L; i++){
			int r = A_col_idx(i);
			if (i >= A_row_ptr(j)){
				c++;
				j++;
			}
			rate_matrix(r,c) = A_val(i);
		}
	}
	

	vector<double> x_initial_vec;
	vector<double> best_vec;
	mat x_initial;
	mat best_mat;
	if (model != -1){
		x_initial_vec = rate_vec;
		best_vec = rate_vec;
	}else{
		x_initial = rate_matrix;
		best_mat = rate_matrix;
	}
	
	
  	
  	//Simulated annealing 
	vector<double> prog_best_lik;
	printf("Simulated annealing started \n");
	if (model != -1){
		simulated_annealing(x_initial_vec, best_vec, L, before, after, frequ, model, prog_best_lik, denom);
	}else{
		simulated_annealing_minus_1(x_initial, best_mat, L, before, after, frequ, model, prog_best_lik, denom);
	}

	mat final_trans_matrix(pow(2,L), pow(2,L), fill::zeros);
	if (model !=-1){
		build_trans_matrix(L, best_vec, model, final_trans_matrix);
	}else{
		final_trans_matrix = best_mat;
	}
	
	//final probability distribution
	mat P_fin(pow(2,L),L+1, fill::zeros);
	P_fin(0,0) = 1;
	for (int t = 1; t<=L; t++){
		P_fin.col(t) = final_trans_matrix * P_fin.col(t-1);
	}

  	//Calculating the fluxes
  	mat fluxes(pow(2,L),pow(2,L), fill::zeros);
  	for(int i = 0; i<pow(2,L); i++){
  		for(int j = 0; j < pow(2,L); j++){
  			double flux = sum(P_fin.row(j))* final_trans_matrix(i,j);
  			fluxes(i,j) = flux;
  		}
  	}

	//creating outputs
  	ofstream outdata;
  	outdata.open("transitions_" + out_name + ".txt");
  	outdata << "From " << "To " << "Probability " << "Flux"<< endl;
  	for(int i = 0; i<pow(2,L); i++){
  		for(int j = 0; j < pow(2,L); j++){
  			if (final_trans_matrix(j,i) != 0){
  				outdata << i << " " << j << " " << final_trans_matrix(j,i) << " " << fluxes(j,i) << endl;
  			}
  		}
  	}
  	outdata.close();

  	ofstream liklh;
  	liklh.open("best_likelihood_" + out_name + ".txt");
  	for(int i = 0; i < prog_best_lik.size();i++){
  		liklh << prog_best_lik[i] << endl;
  	}
  	liklh.close();

	
	// Bootstrapping
	int l = data.size();
	std::random_device rd;
    	std::mt19937 gen(rd());
  	uniform_int_distribution<int> distribution(0, l-1);
  	

  	std::vector<Matrix> all_matrices(bootstrap, Matrix(pow(2,L), std::vector<double>(pow(2,L))));
  	std::vector<Matrix> all_fluxes(bootstrap, Matrix(pow(2,L), std::vector<double>(pow(2,L))));
  	
  	
	
	if (bootstrap != 0){
  		for (int i = 1; i<=bootstrap; i++){
  			cout << "Bootstrap " << i << endl;
  			string d_bt;
  			vector<string> data_bt;
  			for (int j = 0; j<l; j++){
  			 	int ri = distribution(gen);
  			 	d_bt = data[ri];
  				data_bt.push_back(d_bt);
  			}
  			
  			vector<int> frequ_bt;
    			vector<string> reduced_data_bt = red_data(data_bt,frequ_bt);
    			
    			string S_bt;
			vector<string> all_bt;

    			for(int i =0; i < reduced_data_bt.size(); i++){
    				S_bt = reduced_data_bt[i];
    				stringstream ss(S_bt);
    				string word_bt;
        			while (ss >> word_bt) { // Extract word from the stream.
        				all_bt.push_back(word_bt);
    				}
  			}

  			string next_bt;
  			vector<string> before_bt;
  			vector<string> after_bt;
  			for(int j = 0; j < all_bt.size(); j++){
  				next_bt = all_bt[j];
				if (j % 2 == 0){
					before_bt.push_back(next_bt);
				} else{
					after_bt.push_back(next_bt);
				}
  			}

  			vector<double> prog_best_lik_bt;
  			vector<double> best_vec_bt;
  			mat best_mat_bt;
  			mat best_trans_mat_bt(pow(2,L),pow(2,L), fill::zeros);
  			if (model != -1){
  				best_vec_bt = x_initial_vec;
  				simulated_annealing(x_initial_vec, best_vec_bt, L, before_bt, after_bt, frequ_bt, model, prog_best_lik_bt, denom);
  				build_trans_matrix(L,best_vec_bt,model, best_trans_mat_bt);
  			}else{
  				best_mat_bt = x_initial;
  				simulated_annealing_minus_1(x_initial, best_mat_bt, L, before_bt, after_bt, frequ_bt, model, prog_best_lik_bt, denom);
  				best_trans_mat_bt = best_mat_bt;
  			}
			for (int j = 0; j < pow(2,L); j++) {
        	    		for (int k = 0; k < pow(2,L); k++) {
        	        		all_matrices[i-1][j][k] = best_trans_mat_bt(j,k);
        	    		}
        		}
        		
        		//final probability distribution
			mat P_fin(pow(2,L),L+1, fill::zeros);
			P_fin(0,0) = 1;
			for (int t = 1; t<=L; t++){
				P_fin.col(t) = best_trans_mat_bt * P_fin.col(t-1);
			}

  			//Calculating the fluxes
  			mat fluxes(pow(2,L),pow(2,L), fill::zeros);
  			for(int j = 0; j<pow(2,L); j++){
  				for(int k = 0; k < pow(2,L); k++){
  					double flux = sum(P_fin.row(k))* best_trans_mat_bt(j,k);
  					fluxes(j,k) = flux;
  					all_fluxes[i-1][j][k] = fluxes(j,k);
  				}
  			}

				
			//creating outputs
  			ofstream bt;
  			bt.open("bootstrap_"+ to_string(i) + "_" + out_name + ".txt");
  			bt << "From " << "To " << "Probability " << "Flux"<< endl;
  			for(int i = 0; i<pow(2,L); i++){
  				for(int j = 0; j < pow(2,L); j++){
  					if (final_trans_matrix(j,i) != 0){
  						bt << i << " " << j << " " << best_trans_mat_bt(j,i) << " " << fluxes(j,i) << endl;
  					}
  				}
  			}
  			bt.close();
		}

		mat mean(pow(2,L),pow(2,L),fill::zeros);
    		//calculating the mean
    		for (int j = 0; j<pow(2,L); j++){
    			for (int k = 0; k<pow(2,L); k++){
    				double sum = 0;
    				for (int i = 0; i< bootstrap; i++){
    					sum = sum + all_matrices[i][j][k];
    				}
    				double m = sum/bootstrap;
    				mean(j,k) = m;
    			}
    		}
    		
    		mat mean_flux(pow(2,L),pow(2,L),fill::zeros);
    		//calculating the mean
    		for (int j = 0; j<pow(2,L); j++){
    			for (int k = 0; k<pow(2,L); k++){
    				double sum = 0;
    				for (int i = 0; i< bootstrap; i++){
    					sum = sum + all_fluxes[i][j][k];
    				}
    				double m = sum/bootstrap;
    				mean_flux(j,k) = m;
    			}
    		}
	
		ofstream trans_mean;
  		trans_mean.open("mean_" + out_name + ".txt");
  		trans_mean << "From " << "To " << "Probability " << "Flux" << endl;
  		for(int i = 0; i<pow(2,L); i++){
  			for(int j = 0; j < pow(2,L); j++){
  				if((mean(j,i) != 0) || (mean_flux(j,i) != 0)){
  					trans_mean << i << " " << j << " " << mean(j,i) << " " << mean_flux(j,i) << endl;
  				}
  			}
  		}
  		trans_mean.close();
	
  		mat sd(pow(2,L),pow(2,L),fill::zeros);
  		//calculating the standard deviation
  		for (int j = 0; j<pow(2,L); j++){
    			for (int k = 0; k<pow(2,L); k++){
    				double sum = 0;
    				for (int i = 0; i< bootstrap; i++){
    					sum = sum + pow(all_matrices[i][j][k] - mean(j,k),2);
    				}
    				double s = sqrt(sum/(bootstrap - 1));
    				sd(j,k) = s;
    			}
    		}
    		
    		mat sd_flux(pow(2,L),pow(2,L),fill::zeros);
  		//calculating the standard deviation
  		for (int j = 0; j<pow(2,L); j++){
    			for (int k = 0; k<pow(2,L); k++){
    				double sum = 0;
    				for (int i = 0; i< bootstrap; i++){
    					sum = sum + pow(all_fluxes[i][j][k] - mean_flux(j,k),2);
    				}
    				double s = sqrt(sum/(bootstrap - 1));
    				sd_flux(j,k) = s;
    			}
    		}

    		ofstream trans_sd;
  		trans_sd.open("sd_" + out_name + ".txt");
  		trans_sd << "From " << "To " << "Probability " << "Flux " << endl;
  		for(int i = 0; i<pow(2,L); i++){
  			for(int j = 0; j < pow(2,L); j++){
  				if((sd(j,i) != 0) || (sd(j,i) !=0)){
  					trans_sd << i << " " << j << " " << sd(j,i) << " " << sd_flux(j,i) << endl;
  				}
  			}
  		}
  		trans_sd.close();
  		
	}
	return 0;
	
}







