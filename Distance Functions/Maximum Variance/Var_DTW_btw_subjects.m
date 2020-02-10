% Var_DTW_btw_subjects: uses maximum variance to reduce the dimensionality of raw data, and 
% afterward calculates DTW distance on the low-dimensional data for the between-subject case

clear; clc; close;

%% Load the data

% Correct repetitions
Data_NN = csvread('../../Data/Data_Correct.csv');

% Incorrect repetitions
Data_NN_inc = csvread('../../Data/Data_Incorrect.csv');

% Data dimensions
n_dim = 117;
% Number of timesteps
T1 = size(Data_NN,2);
% Number of repetitions
T2 = size(Data_NN,1)/n_dim;

% Tranform the data into cells
% Correct repetitions
Train_Data = cell(1,90);
for i=1:T2
    Train_Data{1,i} = Data_NN((i-1)*n_dim+1:i*n_dim,:)';
end

% Incorrect repetitions
Test_Data = cell(1,90);
for i=1:T2
    Test_Data{1,i} = Data_NN_inc((i-1)*n_dim+1:i*n_dim,:)';
end

%% Extract 3 dimensions with greatest variations

D = n_dim;
n_dim = 3;

% Find the variances for each correct sequence
dim_var=zeros(T2,D);
for i = 1:T2
    dim_var(i,:) = var(Train_Data{i},1);
end

% mean variance 
mean_var = mean(dim_var);

% Sort in descending order
bb = sort(mean_var,'descend');

% Extract the indices of the most varying dimensions
ind =zeros(1,n_dim);
i=1;
while i<=n_dim
    cc = find(mean_var == bb(i));
    % in case that two or above have the same variance
    for h=1:length(cc) 
         ind(1,i) = cc(h); 
    i = i+1;
    end
end

% cancell the redundant data that 
% when i=n_dim and length of cc is two or above
if length(ind)>n_dim
    ind = ind(1,1:n_dim);
end

% reduce correct sequences
Train_Var=cell(1,T2);
for i = 1:T2
    Train_Var{i} = Train_Data{i}(:,ind);
end

% reduce incorrect sequences
Test_Var=cell(1,T2);
for i = 1:T2
    Test_Var{i} = Test_Data{i}(:,ind);
end

%% Calculate the DTW between subjects

% Correct sequences
dtw_train = zeros(1,T2);
for i=1:T2
    for j=1:T2
        dtw_train(i) =dtw_train(i)+dtw(Train_Var{i},Train_Var{j});
    end
end

% Incorrect sequences
dtw_test = zeros(1,T2);  
for i=1:T2
    for j=1:T2
        dtw_test(i) =dtw_test(i)+dtw(Test_Var{i},Train_Var{j});
    end
end

%% Scale and plot data

% Scale data
dtw_train = dtw_train/T1/T2;
dtw_test = dtw_test/T1/T2;

% Scale data in the [1,20] range
MAX = max(max(dtw_train),max(dtw_test));
MIN = min(min(dtw_train),min(dtw_test));
for s = 1: T2
    dtw_train(s) = 19*(dtw_train(s)-MIN)/(MAX-MIN)+1;
    dtw_test(s) = 19*(dtw_test(s)-MIN)/(MAX-MIN)+1;
end

% Plot
h = figure; 
plot(dtw_train,'go','LineWidth',2); hold on, ...
plot(dtw_test,'rs','LineWidth',2);
xlabel('Sequence Number', 'fontsize',18);
ylabel('DTW Distance', 'fontsize',18);
% Title: DTW Distance Between Subjects - Maximum Variance Dimensionality Reduction
title('DTW Distance BS MV','fontsize',18);
legend({'Correct Sequences','Incorrect Sequences'}, 'fontsize',16,'location','NW')
set(gca,'box','off','fontweight','bold','LineWidth',2);
set(gcf,'Units','inches','position',[0 0 5.5 4.5]);
% print(h,'../../Results/DTW_Distance_BS_MV','-dpng','-r300');

%% Separation degree betwween the correct and incorrect sequences
SD = zeros(1,1);
for i=1:T2
    for j=1:T2
        SD=SD+(dtw_test(i)-dtw_train(j))/(abs(dtw_test(i))+abs(dtw_train(j)));
    end
end
SD=SD/T2/T2;
