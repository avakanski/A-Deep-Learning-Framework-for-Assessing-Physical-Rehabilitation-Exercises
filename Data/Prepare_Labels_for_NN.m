% Trains a Gaussian Mixture Model (GMM) on the data with reduced dimensionality
% obtained by the autoencoder neural network, and uses the loglikelihood 
% as a performance indicator
% A scaling function is applied to obtain movement quality scores

clear; close; clc;

addpath('../Utility Functions')

%% Load the data

% Correct repetitions
Data_NN = csvread('Autoencoder_Output_Correct.csv');

% Incorrect repetitions
Data_NN_inc = csvread('Autoencoder_Output_Incorrect.csv');

%% Reshape the data for GMM 

% All individual sequences will be concatenated into one large matrix

% The output dimension of the autoencoder is 4
nDim = 4;

% Number of frames in each movement repetition
L = size(Data_NN,2)/nDim;

% Number of correct sequences
n_seq_corr = size(Data_NN,1);

% Create a row for the time indices of correct data
Data=repmat([1:L],1,n_seq_corr);

% Concatenate the data
Data_position=[];
for i=1:n_seq_corr
    temp = [];
    for j=1:nDim
        temp = [temp; Data_NN(i,j:nDim:nDim*L)];
    end
    Data_position=[Data_position,temp];
end
Data=[Data;Data_position];

% Incorrect sequences

% Number of incorrect sequences
n_seq_inc = size(Data_NN_inc,1);

% Create a row for the time indices of incorrect data
Data_inc = repmat([1:L],1,n_seq_inc);

% Concatenate the data
Data_inc_position=[];
for i=1:n_seq_corr
    temp = [];
    for j=1:nDim
        temp = [temp; Data_NN_inc(i,j:nDim:nDim*L)];
    end
    Data_inc_position=[Data_inc_position,temp];
end
Data_inc = [Data_inc;Data_inc_position];

%% Train GMM

% Define the number of states for GMM
nbStates = 6;

% Number of rows in the data matrix
nbVar = size(Data,1);

% Training by EM algorithm, initialized by k-means clustering
[Priors, Mu, Sigma] = EM_init_regularTiming(Data, nbStates);
[Priors, Mu, Sigma] = EM_boundingCov(Data, Priors, Mu, Sigma);

%% Plot the GMM encoding results

h = figure('position',[20,120,700,700],'name','GMM Autoencoder');
for n=1:nbVar-1
  subplot(nbVar-1,1,n); hold on;
  plotGMM1(Mu([1,n+1],:), Sigma([1,n+1],[1,n+1],:), [1 0.4 0], 1);
  for j=1:n_seq_corr
     plot(Data(n+1,(j-1)*L+1:j*L)','color', [0, 0, 0, 0.1], 'LineWidth', 0.25), hold on
  end
  axis([min(Data(1,:)) max(Data(1,:)) min(Data(n+1,:))-0.01 max(Data(n+1,:))+0.01]);
  set(findobj('type','axes'),'fontsize',10,'box','off')
  xticks(0:20:229)
  yticks(-200:50:200)
  xlabel('Time Frame', 'fontsize',12)
  ylabel('Angle (Degrees)', 'fontsize',12)
end
print(h,'../Results/GMM_Encoded_Movements','-dpng','-r300');

%% Calculate data loglikelihood 

% Correct sequences loglikelihood
for j=1:n_seq_corr
    loglikelihood_corr(j) = -loglik(Data(:,(j-1)*L+1:j*L), nbStates, Priors, Mu, Sigma);
end
 
% Incorrect sequences loglikelihood
for j=1:n_seq_inc
    loglikelihood_inc(j) = -loglik(Data_inc(:,(j-1)*L+1:j*L), nbStates, Priors, Mu, Sigma);
end

%% Plot the scaled loglikelihood for correct and incorrect sequences

MAX = max(max(loglikelihood_corr),max(loglikelihood_inc));
MIN = min(min(loglikelihood_corr),min(loglikelihood_inc));
for s=1:n_seq_corr
    Loglikelihood_sc_corr(s) = 19*(loglikelihood_corr(s)-MIN)/(MAX-MIN)+1;
    Loglikelihood_sc_inc(s) = 19*(loglikelihood_inc(s)-MIN)/(MAX-MIN)+1;
end

% Plot the loglikelihood 
h = figure;
plot(Loglikelihood_sc_corr,'go','LineWidth',2), 
hold on, plot(Loglikelihood_sc_inc,'rs','LineWidth',2)
xlabel('Sequence Number', 'fontsize',18);
ylabel('Loglikelihood', 'fontsize',18);
title('GMM Loglikelihood Scores','fontsize',18);
legend({'Correct Sequences','Incorrect Sequences'}, 'fontsize',16,'location','NW')
set(gca,'box','off','fontweight','bold','LineWidth',2);
set(gcf,'Units','inches','position',[0 0 5.5 4.5]);
print(h,'../Results/GMM_Loglikelihood_Scores','-dpng','-r300');

%% Plot the movement quality scores

% Mean and standard deviation
mean_abs_corr = mean(abs(loglikelihood_corr));
std_abs_corr = std(abs(loglikelihood_corr));

% Scaling function
for r =1:n_seq_corr
        corr_Y(r) = loglikelihood_corr(r)/(mean_abs_corr+3*std_abs_corr);
        inc_Y(r) = loglikelihood_inc(r)/(mean_abs_corr+3*std_abs_corr);
end

for r=1:n_seq_corr
        Correct_Y(r)=1/(1+exp(-3.2+corr_Y(r)));
        Incorrect_Y(r)=1/(1+exp(-3.2+corr_Y(r))+(inc_Y(r)-corr_Y(r))/10);
end

% Plot the movement quality scores
h=figure; 
plot(Correct_Y,'ko', 'MarkerFaceColor','g'), hold on, 
plot(Incorrect_Y,'ks', 'MarkerFaceColor','r')
xlabel('Sequence Number', 'fontsize',18);
ylabel('Quality Score', 'fontsize',18);
axis([0 n_seq_corr 0 1]);
title('GMM Movement Quality Scores','fontsize',18);
legend({'Correct Sequences','Incorrect Sequences'}, 'fontsize',16,'location','SW')
set(gca,'box','off','fontweight','bold','LineWidth',2);
set(gcf,'Units','inches','position',[0 0 5.5 4.5]);
print(h,'../Results/GMM_Movement_Quality_Scores','-dpng','-r300');

%% Save the labels
csvwrite('Labels_Correct.csv',Correct_Y')
csvwrite('Labels_Incorrect.csv',Incorrect_Y')
