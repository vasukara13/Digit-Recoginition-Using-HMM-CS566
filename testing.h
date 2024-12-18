
#include<stdio.h>
#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <limits>
using namespace std;


const int N = 5;     //------------------------ Number of states-------------------------
const int T = 100; // ----------------------Length of observation sequence-----------

double A[N][N]={0.0};
double B[N][32]={0.0};
double PI[N]={0.0};
int observationMatrix[T]={0};
double alpha_values[T][N]={0.0}; // to store alpha matrix

int accuracy=0;


double probabilty[10]={0.0};

void reset_all() {
    //fill(&A[0][0], &A[0][0] + N * N, 0.0);
    //std::fill(&B[0][0], &B[0][0] + N * 32, 0.0);
    //std::fill(PI, PI + N, 0.0);
    double initPI[N]={1.0,0,0,0,0};

    for(int i=0;i<<N;i++){
        PI[i]=initPI[i];
    }
    fill(observationMatrix, observationMatrix + T, 0);
    fill(&alpha_values[0][0], &alpha_values[0][0] + T * N, 0.0);
}

// inputing A,B,PI matrices
void input_A_B_PI_Matrix(const char* name){
    ifstream read_file(name);
    string line;

    if (!read_file) {
        cout << "Cant open file" <<name<< endl;
        exit(0);
    }

    int matrixIndex = 0;  // 1 for A, 2 for B, 3 for PI
    int row = 0;

  
    while (getline(read_file, line)) {
        if (line == "A Matrix") {
            matrixIndex = 1;  // for A matrix 
            row = 0;
        } else if (line == "B Matrix") {
            matrixIndex = 2;  // for B matrix
            row = 0;
        } else if (line == "PI Matrix") {
            matrixIndex = 3;  // for PI matrix
            row = 0;
        } else if (!line.empty()) {
            istringstream iss(line);
            if (matrixIndex == 1) {  
                for (int j = 0; j < 5; ++j) {
                    iss >> A[row][j];
                }
                row++;
            } else if (matrixIndex == 2) {  
                for (int j = 0; j < 32; ++j) {
                    iss >> B[row][j];
                }
                row++;
            } else if (matrixIndex == 3) {  
                for (int j = 0; j < 5; ++j) {
                    iss >> PI[j];
                }
            }
        }
    }

    read_file.close();
}

void load_model(int digit){

    for(int i=0;i<10;i++){
            char name[100];
            sprintf(name,"output/Model/average_model_%d.txt",digit);
            input_A_B_PI_Matrix(name);
        }


}


int predict(){

    double max_probability=-99999999;
    int index=-1;
    for(int i=0;i<10;i++){
        if (probabilty[i]>max_probability){
            max_probability=probabilty[i];
            index=i;
        }
    }
    reset_all();
    return index;


}



void calculate_alpha(int dig_index) {
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
    probabilty[dig_index]=prob;

/*

	cout<<endl;
    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            cout<<setw(10)<<A[i][j];
        }
        cout<<endl;
    }*/

}

void input_obsMatrix(const string& filename) {
    ifstream read_file(filename);
    if (!read_file) {
        cerr << "Error opening obseq file: " << filename << endl;
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


void prerecored_test(){

/*    read one file of 0_30 
    take 1 avg model and get probablity 
    the ncompare the max probability
    and the npredict
    */
    char digits[] = {'0','1','2','3','4','5','6','7','8','9'};

    
    for(int i=0;i<10;i++){
        for(int j=31;j<40;j++){
            char name[100];
            sprintf(name,"output/Observation Sequence/obs_%c_%d.txt",digits[i],j);
            input_obsMatrix(name);

            for(int dig=0;dig<10;dig++){
                load_model(dig);
                calculate_alpha(dig);
            }
            int pred=predict();
            if(i==pred)
                accuracy++;
            cout<<i<<" is PREDICTED AS ====> "<<pred<<endl;
           
        }
    }
    cout<<endl<<"Accuracy of your model is = "<<accuracy<<"%"<<endl;

}
