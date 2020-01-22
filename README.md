# A-Deep-Learning-Framework-for-Assessing-Physical-Rehabilitation-Exercises

The codes in this repository are based on the eponymous research project <a href="https://arxiv.org/abs/1901.10435">A Deep Learning Framework for Assessing Physical Rehabilitation Exercises</a>. The proposed framework for automated quality assessment of physical rehabilitation exercises encompasses metrics for quantifying movement performance, scoring functions for mapping the performance metrics into numerical scores of movement quality, techniques for dimensionality reduction, and deep neural network models for regressing quality scores of input movements via supervised learning. 

# Data
<a href="https://www.webpages.uidaho.edu/ui-prmd/">UI-PRMD dataset</a> of rehabilitation movements is used. It contains full-body skeletal joint displacements for 10 movements performed by 10 healthy subjects. The codes employ 117-dimensional skeletal angles acquired with a Vicon optical tracker for the deep squat exercise. 

# Neural Network Codes
The codes are developed using the Keras library.
* SpatioTemporal NN - the proposed deep spatio-temporal model in the paper
* CNN - a basic convolutional neural network
* RNN - a basic recurrent neural network
* Augmented data models - are based on data augmentation with random noise. Data augmentation was applied in an initial version of the project, but it is not included in the validation results published in the final version of the paper.

# Distance Functions
The codes are developed using MATLAB.
* Maximum Variance - distance functions on reduced-dimensionality data using the maximum variance approach.
* PCA - distance functions on reduced-dimensionality data using PCA.
* Autoencoder - distance functions on reduced-dimensionality data using an autoencoder neural network.
* No Dimensionality Reduction - distance functions on full-body data (117 dimensions).

Please see the List of Files and Functions document for a complete list and brief descriptions of all files in the repository. 

# License
MIT License

# Acknowledgments
This work was supported by the <a href="https://imci.uidaho.edu/get-involved/about-cmci/">Institute for Modeling Collaboration and Innovation (IMCI)</a> at the University of Idaho through NIH Award #P20GM104420.
