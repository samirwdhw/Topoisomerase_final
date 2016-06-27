#include<iostream>
#include<stdlib.h>
#include<math.h>
#include<fstream>

using namespace std;

#define N 5	//SSize of latice
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
#define FILE_NAME "time_non-uniform_ATP1.dat"		//To see where to output data
#define MAX_CATS 100		//Number of catenations to insert initially
#define ATP_MAX 5000	//Maximum ATP till which readings are taken
#define N_RUNS 200		//Number of runs for averaging

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

int main(){

	Cats_list trial;

	trial.insert(1);
	trial.insert(3);
	
	trial.print();

	trial.insert(5);
	trial.insert(2);
	
	trial.print();

	trial.remove(3);
	trial.print();
	trial.remove(2);
	trial.print();


	return 0;
}