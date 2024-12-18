# Digit-Recoginition-Using-HMM-CS566
Below is a sample README for a digit recognition project using Hidden Markov Models (HMM) and LBG k-means clustering, structured similarly to the original README but tailored to the new specifications.

---


## Overview

This project implements digit recognition (0 through 9) using Hidden Markov Models (HMM) and LBG k-means clustering for codebook generation. The system uses Cepstral Coefficients (CI coefficients) for feature extraction and applies the Forward Algorithm for probability calculation and digit prediction. Unlike live speech recognition, this system operates on pre-recorded speech data stored in text files.

The project has two main components:

1. **Training the Model**:
   - Input: A folder of text files (`txt` folder) containing training data (audio features).
   - Output: An observation sequence stored in `output/Observation_sequence` and trained HMM models stored in `output/model`.
   - Codebook generation: LBG k-means clustering is used to create the codebook.
   
2. **Testing the Model**:
   - Input: Pre-recorded speech data (in `txt` files) for testing.
   - Process: Uses the trained models and observation sequences to calculate the probability of the digit being recognized.
   - Accuracy: The system currently achieves about 78% accuracy.

## Directory Structure

```
project_root/
│
├── training/
│   └── 0_1.txt, 0_2.txt, 1_1.txt, 1_2.txt, ...  (Training files for each digit)
│
├── output/
│   ├── Observation_sequence/    (Contains the generated observation sequences)
│   └── model/                   (Contains average HMM models for each digit 0-9)
│       ├── 244101064_0_1.txt
│       ├── 244101064.txt
│       ├── ...
│       └── 244101064_9_0.txt
│
├── testing/
│   └── testing.h                (Contains functions for probability calculation and digit prediction)
│
├── universe.csv                 (Generated feature data - CI coefficients)
│
└── README.md                    (Project documentation)
```

## Setup and Installation

1. **Prerequisites**:
   - C++ Development Environment (e.g., GCC, Clang)
   - Suitable OS (Windows/Linux)
   - Standard C++ libraries (`iostream`, `cmath`, `vector`)
   
2. **No Live Recording**:
   - Since the system does not use live recording, no special audio APIs are required. All inputs are provided via text files.

## Training the Model

### Input
- **Training Data**: Text files (e.g., `0_1.txt`, `0_2.txt`, ..., `9_1.txt`, `9_2.txt`, ...) located in the `training/` folder.
- Each file contains pre-extracted or preprocessed audio features representing spoken digits.

### Process
1. **Preprocessing**:
   - Ensure that the audio features are properly normalized and DC-shift corrected before use.
   
2. **Feature Extraction**:
   - CI coefficients are computed from the feature data.
   - These coefficients represent the short-term power spectrum of the speech frames and are critical for capturing the audio signal characteristics.

3. **Codebook Generation (LBG K-means)**:
   - The Linde-Buzo-Gray (LBG) algorithm is used to perform k-means clustering on the feature vectors.
   - A codebook of a suitable size (e.g., 16 clusters) is generated.
   
4. **Training the HMM**:
   - Using the observation sequences derived from the codebook, HMM parameters (transition, observation, and initial probabilities) are trained.
   - An average model for each digit (0 through 9) is computed and stored in the `output/model` directory.

### Output
- **universe.csv**: Contains the extracted feature vectors (CI coefficients) for all training samples.
- **output/Observation_sequence/**: Contains the observation sequences for each digit.
- **output/model/**: Contains the final average HMM model for each digit (0 to 9).

## Testing the Model

### Input
- Test data (in `.txt` files) similar to the training format, stored in the same `txt` folder or a separate test folder.

### Process
1. **Observation Sequence Generation**:
   - The test audio features are processed to generate an observation sequence using the previously created codebook.

2. **Forward Algorithm**:
   - The Forward Algorithm calculates the probability of the observed sequence given each digit's HMM model.
   - The digit corresponding to the highest probability is chosen as the recognized digit.

3. **Accuracy**:
   - The current system achieves approximately 78% accuracy in digit recognition.

### Output
- The predicted digits for the test cases are recorded or printed as output.
- The `testing.h` file contains the functions for loading the models, computing probabilities, and generating predictions.

## Key Features

- **HMM-based Recognition**: Employs Hidden Markov Models for robust statistical modeling of digit utterances.
- **LBG K-means Clustering**: Generates a concise and representative codebook of feature vectors.
- **Cepstral Coefficients**: CI coefficients provide effective feature representation for spoken digits.
- **Model Storage**: Average HMM models for each digit are easily accessible in `output/model/`.
- **Observation Sequence Utilization**: The observation sequences are a critical link between the raw features and the HMM-based probability calculations.

## Limitations

- **Accuracy**: Current accuracy is around 78%, which may require further optimization, more training data, or enhanced feature extraction techniques.
- **Noise Sensitivity**: Models might not perform well in highly noisy or reverberant environments due to the limited robustness of the training data and features.
- **No Live Processing**: The system only works on pre-recorded data and does not perform real-time recognition.

## To Run the System

1. **For Training the Model**:
   - Place the digit `.txt` files (for each digit 0-9) in the `training/` folder.
   - Run the training script to generate `universe.csv`, observation sequences, and average models for each digit.

2. **For Testing the Model**:
   - Place the test `.txt` files in the appropriate folder.
   - Use the `testing.h` functions to compute the probability for each test sequence and output the predicted digits.

## License

This project is licensed under the IIT Guwahati License - By [Vasu Kara (244101064)].
