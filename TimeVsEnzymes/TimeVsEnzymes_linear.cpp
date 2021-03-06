//This is to simulate the working of Topoisomerase enzyme which
//removes catenations from DNA which is necessary while replication
//of DNA. 

//Enzymes randomly surround medium of DNA, we pick latice points
//equal to number of points and check if the catenation exists
//and if it does whether it can be removed or not

//This is to see the relation of time taken to solve all catenations
//vs the number of enzymes in the system in a 
//chromosome like scenario

//With a linearly decreasing Force Profile along the chromosome

#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<fstream>

using namespace std;

#define N 10000	//SSize of latice
#define FORCE 100.0	//Force applied at kinetochore (pN)
#define V_MAX	3.38  //Max number of cycles per second at a given force (s)
#define Km	270 //Michaelis Constant for the enzyme (uM)
#define Kb	1.38	//Boltzmann Constant (10^-23 units)
#define T 300	//Temperature (K)
#define E_PRODUCT	Kb*T //To increase speed
#define TIME_STEP 0.01		//Time step
#define DELTA 0.241	//angie's parameter (pN), 0.835 is parameter of motion
#define T_MAX 200000		//Seconds
#define FILE_NAME "timeVsATP_linear.dat"		//To see where to output data
#define MAX_CATS 1000		//Number of catenations to insert initially
#define ATP_MAX 2000	//Maximum ATP till which readings are taken
#define N_RUNS 20		//Number of runs for averaging
#define P_FIND 0.02	//Probabilty of landing on a catenation
#define ENZYMES_MAX 100 //Maximum number of enzymes in the system


int n_cats = N;	//No. of catenations in the latice

int track[N];

double prob;	//Probabability of resolving a catenation
double f_each; 	//Force experienced by each catenation
float ATP_conc = ATP_MAX;		//Concentration of ATP

float MAX_ENZYMES;

struct Cats_list{


	int a[MAX_CATS];
	int q;
	int cats_left;		//Required for the linear force profile
	int cats_right;		//Required for the linear force profile
	long int left_sum;
	long int right_sum;


	Cats_list(){

		q = 0;
		cats_left = 0;
		cats_right = 0;
		left_sum = 0;
		right_sum = 0;

	}

	void initialize(){

		q = 0;
		cats_left = 0;
		cats_right = 0;
		left_sum = 0;
		right_sum = 0;

	}

	void insert(int pos){

		a[q] = pos;

		if(pos <= N/2){
			
			cats_left++;
			left_sum += pos;
		
		}
		
		else{
			
			cats_right++;
			right_sum += pos;
		
		}

		q++;

	}

	void remove(int pos){

		int i;

		if( pos <= N/2){
			
			cats_left--;
			left_sum -= pos;

		}
		else{

			cats_right--;
			right_sum -= pos;
		}

		for(i = 0; i<q; i++){

			if(a[i] == pos){

				break;
			}

		}

		track[pos] = 0;

		n_cats--;

		q--;

		for(int j = i; j<q; j++){

			a[j] = a[j+1];

		}

	}

	int givePos(){

		if(q ==0 ){

			cout<<"Here"<<endl;

		}

		int list_place = rand()%q;

		return a[list_place];

	}

	void print(){

		for(int i = 0; i<q; i++){
			cout<<a[i]<<" ";
		}
		cout<<endl;

		cout<<"Left: "<<cats_left<<endl;
		cout<<"Right: "<<cats_right<<endl;

		cout<<"Left Sum: "<<left_sum<<endl;
		cout<<"Right Sum: "<<right_sum<<endl;

	}

};

Cats_list list;		//To maintain the positions of catenations

void calcProb(int pos){	//To update the probability

	if( pos <= N/2){

		f_each = FORCE*list.cats_left*pos/(n_cats*list.left_sum);

	}

	else{
			
		f_each = FORCE*list.cats_right*(N - pos)/(n_cats*(N*list.cats_right - list.right_sum));
	}
	

	float v = 2.937*(exp(-0.505*f_each)+0.2754)*((float)ATP_conc/(ATP_conc + Km));		//Adjust 100 if order changes

	prob = v*TIME_STEP;
}


void fill_cat(){	//To add randomly placed catenations

	for(int i = 0; i< MAX_CATS; i++){
		
		int pos = rand()%N;

		if(track[pos] == 1){

			i--;
			continue;
		}
		
		else{

			track[pos] = 1;
			list.insert(pos);

			//cout<<"here2"<<endl;
		}

	}	

	n_cats = MAX_CATS;

}

float rand1(){		//To generate values b/w 0,1

	float val = (float)rand()/(RAND_MAX);

	return val;
}

void initialize(int a[], int n){	//To make all values of an array zero

	for(int i = 0; i<n; i++){

		a[i] = 0;

	}

	list.initialize();


}

void work(float time1){	//To simulate a single time step 

	for(int i = 0; i< MAX_ENZYMES; i++){

		if(n_cats == 0){
			break;
		}

		int pos = list.givePos();

		calcProb(pos);	//To calculate the probability now
	

		//cout<<prob<<endl;

		float prob_find = rand1();

		if(prob_find < P_FIND){
			
			if( rand1() <prob){

				list.remove(pos);

			}


		}

	}

}


void print(int a[], int n){		//To display an array

	for(int i = 0; i<n; i++){

		cout<<a[i]<<" ";

	}

	cout<<endl;

}


int main(){

	ofstream f1;	//To output data

	f1.open(FILE_NAME);	

	//print(track, N);

	//cout<<"here1";

	for(MAX_ENZYMES = 1 ; MAX_ENZYMES <= ENZYMES_MAX ; MAX_ENZYMES++){	//Never start from 0

		//cout<<ATP_conc;

		float avg_time = 0;

		for(int runs = 0; runs< N_RUNS; runs++){

			//cout<<"Runs"<<" "<<runs<<endl;

			float time1; //To see the time elapsed

			initialize(track, N);	//To make all entries 0

			fill_cat();	//To insert catenations

			for(time1 = TIME_STEP; ;time1 += TIME_STEP){

/*
				if(runs == 0){

					cout<<"Time: "<<time1<<" Cats: "<<n_cats<<endl;
				}
*/

				if(n_cats == 0){
					break;
				}

				work(time1);

			}

			avg_time += time1;

		}

		avg_time /= N_RUNS;

		f1<<MAX_ENZYMES<<" "<<avg_time<<endl;
		cout<<"ATP: "<<ATP_conc<<" Time: "<<avg_time<<" Force: "<<FORCE<<" Enzymes: "<<MAX_ENZYMES<<endl;

	}	


	f1.close();

	return 0;
}