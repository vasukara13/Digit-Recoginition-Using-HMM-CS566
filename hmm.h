// hmm test.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
// 244101064_HMM3.cpp : Defines the entry point for the console application.
//
//#include "stdafx.h"
#include<stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
using namespace std;

#include <stdio.h>
#include <stdlib.h>
#include "testing.h"

//const int N = 5;     //------------------------ Number of states-------------------------
//const int T = 100; // ----------------------Length of observation sequence-----------



//double A[N][N]={0.0};
//double B[N][32]={0.0};
//double PI[N]={0.0};

double A_avg[N][N]={0.0};
double B_avg[N][32]={0.0};
double PI_avg[N]={0.0};

//int observationMatrix[T]={0};

int state_seq[T];
double beta_values[T][N]={0.0};
//double alpha_values[T][N]={0.0}; // to store alpha matrix
double gamma[T][N]={0.0};
double zai[T][N][N]={0.0};


void output_A_B_PI_Matrix(int i) {
    char name[120];
    sprintf(name, "output/Model/average_model_%d.txt", i);
    FILE *write_file = fopen(name, "w");

    if (!write_file) {
        printf("Cannot open file for writing\n");
        exit(0);
    }

    // Write A Matrix
    fprintf(write_file, "A Matrix\n");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            fprintf(write_file, "%e ", A_avg[i][j]); // Scientific notation
        }
        fprintf(write_file, "\n");
    }
    fprintf(write_file, "\n");

    // Write B Matrix
    fprintf(write_file, "B Matrix\n");
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 32; ++j) {
            fprintf(write_file, "%e ", B_avg[i][j]); // Scientific notation
        }
        fprintf(write_file, "\n");
    }
    fprintf(write_file, "\n");

    // Write PI Matrix
    fprintf(write_file, "PI Matrix\n");
    for (int j = 0; j < N; ++j) {
        fprintf(write_file, "%lf ", PI[j]); // Scientific notation
    }
    fprintf(write_file, "\n");

    fclose(write_file);
    // printf("Matrices written to file successfully.\n");
}


double viterbie_algo() {

    // Ccreating delta and psi arrays
    double delta[T][N]={0.0} ;
    int psi[T][N]={0}; //stores arg max index

    // setting first row at t=1
    for (int i = 0; i < N; ++i) {
        delta[0][i] = PI[i] * B[i][observationMatrix[0]-1];
    }

    for (int t = 1; t < T; ++t) {
        for (int j = 0; j < N; ++j) {
            double max_probability = -9999999;
            int best_state_sequence = 0;

            //max function instead
            for (int i = 0; i < N; ++i) {
                double prob = delta[t - 1][i] * A[i][j];
                if (prob > max_probability) {
                    max_probability = prob;
                    best_state_sequence = i;
                }
            }
            delta[t][j] = max_probability * B[j][observationMatrix[t]-1];
            psi[t][j] = best_state_sequence;
        }
    }

    // backtracking to find the most probable state sequence
   
    double max_probability = -9999999999; // -99999999999
    int last_state = 0;

    for (int i = 0; i < N; ++i) {
        if (delta[T - 1][i] > max_probability) {
            max_probability = delta[T - 1][i];
            last_state = i;
        }
    }
    return max_probability;

    
}

// calculating alpha values using forward algorithm

void calculate_alpha() {
    // initilize alpha_values[0][i] = PI[i] * B[i][O[0]]
    for (int i = 0; i < N; ++i) {
        alpha_values[0][i] = PI[i] * B[i][observationMatrix[0]-1];
    }

    for (int t = 1; t < T; ++t) {
        for (int j = 0; j < N; ++j) {
            double sum = 0.0;
			alpha_values[t][j]=0;
            for (int i = 0; i < N; ++i) {
                alpha_values[t][j]+= alpha_values[t - 1][i] * A[i][j];
            }
            int temp_obs=observationMatrix[t];
            alpha_values[t][j] *= B[j][temp_obs-1];
        }
    }
    // summation of last row to get probability 
    double prob = 0.0;
    for (int i = 0; i < N; ++i) {
        prob += alpha_values[T - 1][i];
    }

/*	cout<<endl;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            cout<<setw(10)<<A[i][j];
        }
        cout<<endl;
    }*/

}

/*
void input_obsMatrix(const string& filename) {
    ifstream read_file(filename);
    if (!read_file) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }
    string line;
    int count = 0;
    while (getline(read_file, line) && count < T) {
        istringstream iss(line);
        while (iss >> observationMatrix[count] && count < T) {
            ++count;
        }
    }
    read_file.close();
}
*/
// Calcualting beta values using backward algorithm
void computeBeta() {

    //setting last row to 1
    for (int i = 0; i < N; ++i) {
        beta_values[T - 1][i] = 1.0;
    }
    // claculating rest of the beta_values array
    for (int t = T - 2; t >= 0; --t) {
        for (int i = 0; i < N; ++i) {
            beta_values[t][i] = 0.0;
            for (int j = 0; j < N; ++j) {
                int temp_obs=observationMatrix[t + 1];
                beta_values[t][i] += A[i][j] * B[j][temp_obs-1] * beta_values[t + 1][j];
            }
        }
    }
   /* 
    cout<<endl;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            cout<<B[i][j]<<"    ";
        }
        cout<<endl;
    }*/  
}


void reset_avg(){   //--------------------------------- TO RESET THE AVERAGE MATRIX--------------------------------------

    
    // Copy values to A
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A_avg[i][j] = 0.0;
        }
    }
    
    
    // Copy values to B
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 32; ++j) {
            B_avg[i][j] = 0.0;
        }
    }
    double initPI[N]={1.0,0,0,0,0};

    for(int i=0;i<N;i++){
        PI_avg[i]=initPI[i];
    }

}

void reset_A_B(){           //---------------------------------to reset the a and b and GET AVG MODLES VALUES------------------------------------------------------//
    // Initializing A matrix

/*
	for (int i = 0; i < N; ++i) {
		PI_avg[i]+=PI[i];
	}*/


    double initA[N][N] = {
        {0.8, 0.2, 0.0, 0.0, 0.0},
        {0.0, 0.8, 0.2, 0.0, 0.0},
        {0.0, 0.0, 0.8, 0.2, 0.0},
        {0.0, 0.0, 0.0, 0.8, 0.2},
        {0.0, 0.0, 0.0, 0.0, 1.0}
    };
    
    // Copy values to A
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = initA[i][j];
        }
    }
    
    // Initializing B matrix
    double initB=0.03125;
    
    // Copy values to B
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 32; ++j) {
            B[i][j] = initB;
        }
    }
    double initPI[5]={1.0,0,0,0,0};

    for(int i=0;i<N;i++){
        PI[i]=initPI[i];
    }
}



void gamma_calculater() {

    //-------------------------------------ZAI matrix-------------------------------------------

    for (int t = 0; t < T-1; ++t) {
        double denomenator=0;
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                denomenator += B[j][observationMatrix[t + 1]-1] * beta_values[t + 1][j]*alpha_values[t][i] * A[i][j] ;
            }
        }
        for (int i = 0; i < N; i++) {
            gamma[t][i]=0;
            for (int j = 0; j < N; j++) {
                zai[t][i][j] = (alpha_values[t][i] * A[i][j] * B[j][observationMatrix[t + 1]-1] * beta_values[t + 1][j]) / denomenator;           
                gamma[t][i]+=zai[t][i][j];
            }
        }  
    }
    //--------------------------------------GAMMA matrix------------------------------------
    
  /*  double summm=0;
    for(int t=0;t<T;t++){
        for(int i=0;i<N;i++){
            for(int j=0;j<N;j++){
                summm+=zai[t][i][j];
            }
            gamma[t][i]=summm;
            summm=0;
            
        }
    }*/
    
}
void best_model(){

    //calculate the best a and B

}

void avg_add(){

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A_avg[i][j]+=A[i][j];
        }
    }

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < 32; ++j) {
            B_avg[i][j]+=B[i][j];
        }
    }

}

void re_estimation_algortihm(){
 
    // ----------------------------------------A_Bar-----------------------------------------
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
                double n=0;
                double d=0;
            for(int t =0;t<T-1;t++){
                n+=zai[t][i][j];
                d+=gamma[t][i];
            }
            A[i][j]=n/d;
            
            
        }
        //A[i][i] += (1 - sum_adjust);
        
    }

    // ----------------------------------------B_Bar-----------------------------------------
  
    for(int i=0;i<N;i++){
        
        for(int j=0;j<32;j++){
            long double n=0;
            long double d=0;
            for(int t =0;t<T;t++){
                if(observationMatrix[t]-1==j){
                    n+=gamma[t][i];
				}
                d+=gamma[t][i];
            }
            B[i][j]=n/d;
               
        }
    }
    int count =0;
    int max_i=0;
    int max_j=0;
    double max=0;
    for(int i=0;i<N;i++)
    {
        for(int j=0;j<32;j++)
        {
            if(B[i][j]==0)
            {
                B[i][j] = 1.000000e-030;
                count = count + 1;
            }
            if(max<B[i][j])
            {
                max_i=i;
                max_j=j;
            }
        }
    }
    B[max_i][max_j] = B[max_i][max_j] - (count *1.000000e-030);

}

void calulate_avg_model(int digit){

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            A_avg[i][j]/=30;
        }
    }


    for(int i=0;i<N;i++){
        for(int j=0;j<32;j++){
            B_avg[i][j]/=30;
        }
    }

    for(int i=0;i<N;i++){
        PI_avg[i]/=30;
    }

    output_A_B_PI_Matrix(digit);
    reset_avg();


}


int hmm(){ 

    char digits[] = {'0','1','2','3','4','5','6','7','8','9'};

//================================================Training OF DIGITS==================================//
    
    
    for(int digit=0;digit<10;digit++) // should be 5 for digits
    {
		reset_A_B();
        printf("Training of digit ====>  %c \n",digits[digit]);
        cout<<"<";
        for(int j=1;j<=30;j++)   // should be 30 for training
        {
            cout<<"==";
            char output_name[150]="";
            char observation_input_name[150]="";

            sprintf(observation_input_name,"output/Observation Sequence/obs_%c_%d.txt",digits[digit],j);
            sprintf(output_name,"Model/Model%c_%d_ai.txt",digits[digit],j);

            input_obsMatrix(observation_input_name);
            
            int i=1;
        	double prev=2.2250738585072014e-308;
        	double current=2.2250738585072014e-308;
            //Scout<<"starting";
            while(i<100){
                prev=current;
                calculate_alpha() ;
                computeBeta();
                gamma_calculater();
                re_estimation_algortihm();
        		current=viterbie_algo();
				i++;
                //cout << "\nProbability of the P * : " << prev;
				if(prev>current){
                    break;
                }
                
                
            }
			
			avg_add();
			reset_A_B();       // caluclates plus avg and resets array also
        }
        cout<<">"<<endl;
        calulate_avg_model(digit);
        reset_avg(); //resst best also

    }  

    reset_all();
    prerecored_test();
    return 0;

}


