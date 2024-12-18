// 244101064_digitrecoginition.cpp : Defines the entry point for the console application.
//

#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <iostream>
#include <vector>
#include <limits>
#include <fstream>
#include <sstream>
#include <string>

#include "lbg.h"
#include "hmm.h"

using namespace std;

#define max_length 32000 //max signal length
#define FRAME_LENGTH 320  // Frame size for analysis
#define LPC_ORDER 12      // LPC analysis order
#define pi 3.14159265358979323846
const int num_frames = 100;//100  78

void get_hamming_window(double *window, int N) {
        //for(int j=0;j<5;j++)
        for (int n = 0; n < N; ++n) {   //  N should be 320
                window[n] = 0.54 - 0.46 * cos((2 * pi * n) / (N - 1));
        }
}

void apply_hamming(double *frames, double window[320], int N) {// hamming windiow
    for (int n = 0; n < N; ++n) {
        frames[n] *= window[n];
    }
}

void write_data_universe(const char* filename,double Ci_values[num_frames][13]){//-----------------------------write data to universe-------------------------
    FILE* file = fopen(filename, "a+");
    if (!file) {
        printf("Unable to open file output %s\n", filename);
        exit(1);
    }
    for(int j=0;j<num_frames;j++){ //101

		//fprintf(file, "Frame [%d] =  ", (j+1));
        for (int i = 1; i < 13; i++) {
            //fprintf(file, "%lf \t",Ci_values[j][i]);
            fprintf(file, "%lf,",Ci_values[j][i]);
        
        }

            fprintf(file, "\n");
    }
    fclose(file);
}




void read_signal(const char* filename, double* signal_data, int* signal_size) {
    FILE* file = fopen(filename, "r");
    if (!file) {
        printf("Unable to open file %s\n", filename);
        exit(1);
    }

    *signal_size = 0;
    while (fscanf(file, "%lf", &signal_data[*signal_size]) != EOF && *signal_size<40000) {
        (*signal_size)++;
    }
	
    fclose(file);
}

void write_data(const char* filename, double* data_arr, int order,int frame_position) {
    FILE* file = fopen(filename, "a+");
    if (!file) {
        printf("Unable to open file %s\n", filename);
        exit(1);
    }
	fprintf(file, "Frame %d:", frame_position);
	fprintf(file, "\n");
    for (int i = 1; i <= order; i++) {
        fprintf(file, "a[%d] = %lf \n", i,data_arr[i]);
    }
    fprintf(file, "\n");

    fclose(file);
}

void removeDCShift(double* signal, int length) { //-------------------------remove dc shift-------------------------
    double sum = 0.0;
    for (int i = 0; i < length; i++) {
        sum += signal[i];
    }
    double mean = sum / length;
    //printf("DC offset is : %f \t",mean);
    for (int i = 0; i < length; i++) {
        signal[i] -= mean;
    }
}

void normalizeSignal(double* signal, int length){  //-------------------------normalize signal------------------------


    double min_val = signal[0];
    double max_val = signal[0];

    for (int i = 1; i < length; i++) {
        if (signal[i] < min_val) {
            min_val = signal[i];
        }
        if (signal[i] > max_val) {
            max_val = signal[i];
        }
    }
    for (int i = 0; i < length; i++) {
        signal[i] = -5000 + ((signal[i] - (min_val)) / (max_val - min_val)) *10000;
    }
}



//---------------------------------Calculate autocorrealtion matrix-------------------------------------//

void Ri_calc(double frame_data[320], int frame_length, double* Ri_values, int order) {
    for (int i = 0; i <= order; i++) {
        Ri_values[i] = 0.0;
        for (int j = 0; j < frame_length - i; j++) {
            Ri_values[i] += frame_data[j] * frame_data[j + i];
        }
        //cout<<Ri_values[i]<<endl;
    }
}
//---------------------------------------Calculate AI matrix-------------------------------------------//

void levinson_durbin_recursion(double* Ri_values, double* Ai_values, int order) {
    double E[LPC_ORDER + 1]={0.0};
    double Alpha[LPC_ORDER + 1][LPC_ORDER + 1]={0.0};
    double k[LPC_ORDER + 1]={0.0};

    E[0] = Ri_values[0];
    for (int i = 1; i <= order; i++) {
        double sum = 0.0;

        for (int j = 1; j < i; j++) {      // i-1
            sum += Alpha[i - 1][j] * Ri_values[i - j];
        }
        if(E[i-1]!=0)
            k[i] = (Ri_values[i] - sum) / E[i - 1];
        else
            k[i]=0 ;  //(Ri_values[i] - sum); // or zero check this


        Alpha[i][i] = k[i];
        for (int j = 1; j <= i - 1; j++) {
            Alpha[i][j] = Alpha[i - 1][j] - k[i] * Alpha[i - 1][i - j];
        }
        E[i] = (1 - k[i] * k[i]) * E[i - 1];
    }
    for (int i = 1; i <= order; i++) {  //i =1
        Ai_values[i] = Alpha[order][i];  //alpha[0][i]
        
    }
}
//---------------------------------Calculate cepstral matrix-----------------------------------------//
void cepstral_coefficients(double *Ci_values,double* Ri_values,double* Ai_values, int p) {
    //Ci_values[0] += log(Ri_values[0]); // First cepstral coefficient
    

    //----------log of r[0] try this also
    //Ci_values[0] = Ai_values[0];
    for (int i = 1; i <= p; i++) {
        //double sum = 0.0;
        Ci_values[i] = Ai_values[i];
        for (int k = 1; k< i; k++) {// maybe <=
            //sum+= ((double)k / (double)i) * Ci_values[k] * Ai_values[i - k];
            Ci_values[i]+= ((double)k / (double)i) * Ci_values[k] * Ai_values[i - k];
        }  
    }
    for (int i = 1; i <= p; i++) {
        Ci_values[i] *= (1 + (p / 2.0) * sin((pi * i) / p));
    }


  }

void write2DArrayTo1DFile(int array[num_frames][320], const string& filename) {
    ofstream file(filename);

    if (!file.is_open()) {
        cerr << "Error: Could not open file for writing" << endl;
        return;
    }

    for (int i=0;i<num_frames;i++) {
        for (int j=0;j<320;j++) {
            file << array[i][j] << "\n";  // Write values in 1D form separated by space
        }
    }

    file.close();
    cout << "Array written to " << filename << " successfully." << endl;
	exit(0);
}


/*void get_ci(){

    DCshift(data_r);
    normalize(data_r);
    frame_select(data_r,frames_r);
    h_window(frames_r);
    tapered_window();

        // Calculate Ri, Ai, and Ci for each frame
    for (int i = 0; i < FRAMES; i++) {
        autocorrelation(frames_r, Ri_r, i);
        levinson(i,Ri_r,Ai_r);
    }

        // Compute Ci values using the new function
    calculate_ci(frames_r,Ci_r,Ai_r);
}*/


int main(){

    char digits[] = {'0','1','2','3','4','5','6','7','8','9'};
    
//=========================================PREPROCESSING OF DIGITS==================================//

    for(int digit=0;digit<10;digit++) // should be 5 for digits
    {
        printf("Pre_Processing for digit ====>  %c \n",digits[digit]);
        cout<<"<";
        for(int j=1;j<=40;j++)// should be 40
        {
            cout<<"=";
            double Ci_values[num_frames][13]={0.0};
            double signal[max_length]={0};

			int signal_size=0;
            char input_name[num_frames]="";
            sprintf(input_name,"txt/244101064_E_%c_%d.txt",digits[digit],j);
            char output_name[num_frames]="";
            sprintf(output_name,"output/Ai's/244101064_%c_%d_ai.txt",digits[digit],j);

            //-----------------------------------PROCESSING STARTING-------------------------

            read_signal(input_name, signal, &signal_size); // Read the speech signal from a file     
            removeDCShift(signal, signal_size); // Remove DC offset
            normalizeSignal(signal, signal_size); // Normalize the signal
            
            double window[num_frames][FRAME_LENGTH];
            // ---------------------------------------FRAME MAKING-------------------------

            double frames[num_frames][320]={0};
            int marker = 0;
            for (int j = 0; j < num_frames; j++) {
                for (int f = 0; f < FRAME_LENGTH; f++) {
                    frames[j][f] = signal[marker];
                    marker++;
                }
            }
            //write2DArrayTo1DFile(frames, "output.txt");

            //-----------------------------------FRAMES PROCESSING of 1 FILE-----------------------------
            
            for (int i = 0; i <num_frames; i++) {
                
				double Ri_values[13]={0.0};
				double Ai_values[13]={0.0};
                
                get_hamming_window(window[i],FRAME_LENGTH);  // CALCULATION OF HAMMING WINDOW
                apply_hamming(frames[i],window[i],FRAME_LENGTH); // APPLYING IT TO EACH FRAME
                
                Ri_calc(frames[i], 320, Ri_values, LPC_ORDER); // CALCULATING AUTOCORRELATION COEFFICIENTS
                levinson_durbin_recursion(Ri_values, Ai_values, LPC_ORDER); //CALCULATING LPC COEFFICIENTS
                cepstral_coefficients(Ci_values[i],Ri_values,Ai_values,LPC_ORDER);
                

            }
            char ci_name[num_frames]="";
            sprintf(ci_name,"output/Ci's/244101064_%c_%d_ci.txt",digits[digit],j);  
            //write_data_universe(ci_name,Ci_values);
            write_data_universe("universe.csv",Ci_values);

        }

        cout<<">"<<endl<<endl;
    }
	//removeZeroRows("universe.csv");

	cout<<"Starting LBG kMeans...............\n\n";
    read_csv_to_float_array("universe.csv");
    lbg_kmeans();
    cout<<"\nWriting codeBook........\n";
    write_csv();
	cout<<endl;

	hmm();

	exit(0);
}




