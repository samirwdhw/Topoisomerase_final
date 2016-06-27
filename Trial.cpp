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
#define FILE_NAME "time_non-uniform_ATP1.dat"		//To see where to output data
#define MAX_CATS 100		//Number of catenations to insert initially
#define ATP_MAX 5000	//Maximum ATP till which readings are taken
#define N_RUNS 200		//Number of runs for averaging

struct Cats_list{


	int a[MAX_CATS];
	int q;

	Cats_list(){

		q = 0;

	}

	void insert(int pos){

		a[q] = pos;
		q++;

	}

	void remove(int pos){

		int i;

		for(i = 0; i<q; i++){

			if(a[i] == pos){

				break;
			}

		}

		q--;

		for(int j = i; j<q; j++){

			a[j] = a[j+1];

		}

	}

	void print(){

		for(int i = 0; i<q; i++){
			cout<<a[i]<<" ";
		}

		cout<<endl;

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