% Mahalanobis_wthn_subjects: calculates Mahalanobis distance on the raw data for the within-subject case

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

% rt is the sequence of repetion for each person
rt=[9 9 9 9 9 10 8 9 8 10];

%% Mahalanobis distance based on variances of each dimension

% Find the variances for each correct sequence
dim_var=zeros(T2,n_dim);
for i = 1:T2
    dim_var(i,:) = var(Train_Data{i},1);
end

% Calculate the mean variance for all correct sequences
% to be used for calculating Mahalanobis distance
Mhlnbs_dist = mean(dim_var);
Mhlnbs_cov=diag(Mhlnbs_dist);

%% M distance for the correct sequences

rms_train = zeros(1,T2);
M=0;
for s=1:10
    for i=1:rt(1,s)
        for j=1:rt(1,s)
        for t = 1:T1
            rms_train(M+i) = rms_train(M+i)+sqrt((...
            Train_Data{M+i}(t,:)-Train_Data{M+j}(t,:))*Mhlnbs_cov*...
           (Train_Data{M+i}(t,:)-Train_Data{M+j}(t,:))');
        end
        end
        % scale data
        if rt(1,s)>0
           rms_train(M+i)=rms_train(M+i)/rt(1,s);
        else
           rms_train(M+i)=rms_train(M+i);
        end
    end
    M=M+rt(1,s);
end


%% Find the variances for each incorrect sequence 

dim_var_test=zeros(T2,n_dim);
for i = 1:T2
    dim_var_test(i,:) = var(Test_Data{i},1);
end

%% M distance for the incorrect sequences

rms_test = zeros(1,T2);
N=0;
for s=1:10
    for i=1:rt(1,s)
        Mhlnbs_cov_test = diag(dim_var_test(i,:));
        for j=1:rt(1,s)
            for t = 1:T1
                rms_test(N+i) = rms_test(N+i)+sqrt((...
                Test_Data{N+i}(t,:)-Train_Data{N+j}(t,:))*...
                Mhlnbs_cov_test*(Test_Data{N+i}(t,:)-...
                Train_Data{N+j}(t,:))');
            end
        end
        % scale data
        if rt(1,s)>0
           rms_test(N+i)=rms_test(N+i)/rt(1,s);
        else
           rms_test(N+i)=rms_test(N+i);
        end
    end
    N=N+rt(1,s);
end

%% Scale and plot data

% Scale data
rms_train = rms_train/T1;
rms_test = rms_test/T1;

% Scale data in the [1,20] range
MAX = max(max(rms_train),max(rms_test));
MIN = min(min(rms_train),min(rms_test));
for s = 1: T2
    rms_train(s) = 19*(rms_train(s)-MIN)/(MAX-MIN)+1;
    rms_test(s) = 19*(rms_test(s)-MIN)/(MAX-MIN)+1;
end

% Plot
h = figure; 
plot(rms_train,'go','LineWidth',2); hold on, ...
plot(rms_test,'rs','LineWidth',2);
xlabel('Sequence Number', 'fontsize',18);
ylabel('Mahalanobis Distance', 'fontsize',18);
% Title: Mahalanobis Distance Within Subjects - No Dimensionality Reduction
title('Mahalanobis Distance WS NDR','fontsize',18);
legend({'Correct Sequences','Incorrect Sequences'}, 'fontsize',16,'location','NW')
set(gca,'box','off','fontweight','bold','LineWidth',2);
set(gcf,'Units','inches','position',[0 0 5.5 4.5]);
% print(h,'../../Results/Mahalanobis_Distance_WS_NDR','-dpng','-r300');

%% Separation degree betwween the correct and incorrect sequences
SD = zeros(1,1);
for i=1:T2
    for j=1:T2
        SD=SD+(rms_test(i)-rms_train(j))/(abs(rms_test(i))+abs(rms_train(j)));
    end
end
SD=SD/T2/T2;

