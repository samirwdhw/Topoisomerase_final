#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<fstream>

using namespace std;

#define N 10000	//SSize of latice
#define MAX_ENZYMES 100		//No. Of enzymes around the latice
#define FORCE 1.0	//Force applied at kinetochore (pN)
#define V_MAX	3.38  //Max number of cycles per second at a given force (s)
#define Km	270 //Michaelis Constant for the enzyme (uM)
#define Kb	1.38	//Boltzmann Constant (10^-23 units)
#define T 300	//Temperature (K)
#define E_PRODUCT	Kb*T //To increase speed
#define TIME_STEP 0.1		//Time step
#define DELTA 0.241	//angie's parameter (pN), 0.835 is parameter of motion
#define T_MAX 200000		//Seconds
#define FILE_NAME "Force_EXPO.dat"		//To see where to output data
#define MAX_CATS 1000		//Number of catenations to insert initially
#define ATP_MAX 5000	//Maximum ATP till which readings are taken
#define N_RUNS 200		//Number of runs for averaging
#define NORMAL 1000.0	//Normalising factor to keep exponent within range

int n_cats = N;	//No. of catenations in the latice

int track[N];

float prob;	//Probabability of resolving a catenation
float f_each; 	//Force experienced by each catenation
float ATP_conc = ATP_MAX;		//Concentration of ATP




struct Cats_list{


	int a[MAX_CATS];
	int q;
	int cats_left;		//Required for the linear force profile
	int cats_right;		//Required for the linear force profile
	double left_sum;
	double right_sum;


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
			left_sum += exp((float)pos/NORMAL);
		
		}
		
		else{
			
			cats_right++;
			right_sum += exp((float)(N-pos)/NORMAL);
		
		}

		q++;

	}

	void remove(int pos){

		int i;

		if( pos <= N/2){
			
			cats_left--;
			left_sum -= exp((float)pos/NORMAL);

		}
		else{

			cats_right--;
			right_sum -= exp((float)(N-pos)/NORMAL);
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

Cats_list list;


void calcProb(int pos){	//To update the probability

	if( pos <= N/2){

		f_each = FORCE*list.cats_left*exp(pos/NORMAL)/(n_cats*list.left_sum);

	}

	else{
			
		f_each = FORCE*list.cats_right*exp((N-pos)/NORMAL)/(n_cats*(list.right_sum));
	}
	

	float v = 2.937*(exp(-0.505*f_each)+0.2754)*((float)ATP_conc/(ATP_conc + Km));		//Adjust 100 if order changes

	prob = v*TIME_STEP;
}

int main(){

	ofstream f1;
	f1.open(FILE_NAME);

	n_cats = MAX_CATS;

	for(int i = 0; i<MAX_CATS; i++){

		list.insert(i*10);

	}

	//list.print();

	for(int i = 0; i<N; i++){

		calcProb(i);

		f1<<i<<" "<<f_each<<endl;

	}

	f1.close();

	return 0;
}