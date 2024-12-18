#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <limits>

using namespace std;

#define SPLITTING_FACTOR 0.03
#define CSV_COLUMN 12
#define CSV_ROWS 40000 //31200
#define K 32
#define num_frame 100 // 78
vector<vector<double>> universe(CSV_ROWS, vector<double>(CSV_COLUMN));
vector<vector<double>> centroid(K, vector<double>(CSV_COLUMN));
vector<int> cluster(CSV_ROWS, 0);
vector<int> cluster_size(K, 0);
vector<double> lbg_centroid(CSV_COLUMN, 0.0);
int obs[CSV_ROWS] = {0};

void write_obs_data() {
    int chunkSize = num_frame;
    int index = 0;
    char digits[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

    for (int x = 0; x < 10 && index < CSV_ROWS; ++x) {
        for (int y = 1; y <= 40 && index < CSV_ROWS; ++y) {
            char obs_name[150];
            sprintf(obs_name, "./output/Observation Sequence/obs_%c_%d.txt", digits[x], y);
            ofstream outfile(obs_name);
            if (!outfile) {
                cerr << "Error opening file " << obs_name << endl;
                return;
            }

            for (int i = 0; i < chunkSize && index < CSV_ROWS; ++i) {
                outfile << obs[index++] << "\n";
            }

            outfile.close();
        }
    }

    cout << "Files generated successfully." << std::endl;
}

double tokura_distance(int index, const vector<double>& point) {
    double dist = 0;
    double tokura_weights[CSV_COLUMN] = {1.0, 3.0, 7.0, 13.0, 19.0, 22.0, 25.0, 33.0, 42.0, 50.0, 56.0, 61.0};
    for (int i = 0; i < CSV_COLUMN; i++) {
        double diff = point[i] - centroid[index][i];
        dist += tokura_weights[i] * (diff * diff);
    }
    return dist;
}

void update_centroids(int x) {
    vector<vector<double>> updated_centroids(x, vector<double>(CSV_COLUMN, 0.0));
    vector<int> count(x, 0);

    if (x == 32) {
        for (int i = 0; i < cluster.size(); i++) {
            obs[i] = cluster[i];
        }
    }

    for (int j = 0; j < CSV_ROWS; j++) {
        if (cluster[j] > 0) {
            int cent = cluster[j] - 1;
            for (int i = 0; i < CSV_COLUMN; i++) {
                updated_centroids[cent][i] += universe[j][i];
            }
            count[cent]++;
        }
    }

    for (int cent = 0; cent < x; cent++) {
        if (count[cent] > 0) {
            for (int i = 0; i < CSV_COLUMN; i++) {
                centroid[cent][i] = updated_centroids[cent][i] / count[cent];
            }
        }
    }

    fill(cluster.begin(), cluster.end(), 0);
}

double distortion_calculator(int x) {
    double distortion = 0.0;
    for (int c = 0; c < x; c++) {
        double dist = 0.0;
        for (int u = 0; u < CSV_ROWS; u++) {
            if (cluster[u] == (c + 1)) {
                dist += tokura_distance(c, universe[u]);
            }
        }
        distortion += dist;
    }
    return distortion / CSV_ROWS;
}

void read_csv_to_float_array(const char* filename) {
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Could not open file " << filename << endl;
        return;
    }

    string line;
    int row_index = 0;
    while (getline(file, line) && row_index < CSV_ROWS) {
        istringstream ss(line);
        string token;
        int col_index = 0;
        while (getline(ss, token, ',') && col_index < CSV_COLUMN) {
            universe[row_index][col_index] = stod(token);
            col_index++;
        }
        row_index++;
    }

    file.close();
}

void split_centroid(int x) {
    for (int cent = 0; cent < x; cent++) {
        for (int i = 0; i < CSV_COLUMN; i++) {
            centroid[cent + x][i] = centroid[cent][i] * (1 - SPLITTING_FACTOR);
            centroid[cent][i] = centroid[cent][i] * (1 + SPLITTING_FACTOR);
        }
    }

    cout << "For K= " << x * 2 << endl;
}

void initial_centroid() {
    for (int i = 0; i < CSV_COLUMN; i++) {
        for (int j = 0; j < CSV_ROWS; j++) {
            lbg_centroid[i] += universe[j][i];
        }
        lbg_centroid[i] /= CSV_ROWS;
        centroid[0][i] = lbg_centroid[i];
    }
}

void lbg_kmeans() {
    int current_codebook_size = 1;

    initial_centroid();

    while (current_codebook_size < K) {
        split_centroid(current_codebook_size);
        current_codebook_size *= 2;

        double prev_distortion = 0.0;
        int iterations = 0;

        double distortion = numeric_limits<double>::max();
        while (fabs(distortion - prev_distortion) > 0.00001) {
            if (iterations != 0)
                prev_distortion = distortion;

            for (int i = 0; i < CSV_ROWS; i++) {
                double min_dist = numeric_limits<double>::max();
                int cluster_index = -1;

                for (int j = 0; j < current_codebook_size; j++) {
                    double dist = tokura_distance(j, universe[i]);
                    if (dist < min_dist) {
                        cluster_index = j + 1;
                        min_dist = dist;
                    }
                }
                cluster[i] = cluster_index;
            }

            for (int index = 0; index < current_codebook_size; index++) {
                cluster_size[index] = 0;
                for (int size = 0; size < CSV_ROWS; size++) {
                    if (cluster[size] == (index + 1)) {
                        cluster_size[index]++;
                    }
                }
            }

            distortion = distortion_calculator(current_codebook_size);
            iterations++;
            cout << "For codebook size " << current_codebook_size << " ====> iteration : " << iterations << " distortion is : " << distortion << endl;
            update_centroids(current_codebook_size);
        }

        cout << endl;
    }
    write_obs_data();
}

void write_csv() {
    FILE *file = fopen("Final_codebook.csv", "w");

    if (file == NULL) {
        fprintf(stderr, "Error: Could not open file for writing\n");
        return;
    }



    for (int i = 0; i < K; i++) {
        for (int j = 0; j < CSV_COLUMN; j++) {
            fprintf(file, "%lf", centroid[i][j]);
            if (j < CSV_COLUMN - 1) {
                fprintf(file, ",");
            }
        }
        fprintf(file, "\n");
    }

    fclose(file);
    printf("Codebook written successfully.\n");
}
