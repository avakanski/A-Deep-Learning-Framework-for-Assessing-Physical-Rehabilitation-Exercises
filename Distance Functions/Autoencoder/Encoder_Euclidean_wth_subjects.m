% Encoeder_Euclidean_wth_subjects: uses autoencoder neural network to reduce the dimensionality of raw data, 
% and afterward calculates Euclidean distance on the low-dimension data for the within-subject case

clear; clc; close;

%% Load the data

% Correct repetitions
Data_NN = csvread('../../Data/Autoencoder_Output_Correct.csv');

% Incorrect repetitions
Data_NN_inc = csvread('../../Data/Autoencoder_Output_Incorrect.csv');

% Data dimensions
nDim = 4;
% Number of timesteps
T1 = size(Data_NN,2)/nDim;
% Number of repetitions
T2 = size(Data_NN,1);

% Tranform the data into cells
% Correct repetitions
Train_Data_Reduced = cell(1,90);
for i=1:T2
    temp = [];
    for j=1:nDim
        temp = [temp; Data_NN(i,j:nDim:nDim*T1)];
    end
    Train_Data_Reduced{1,i} = temp';
end

% Incorrect repetitions
Test_Data_Reduced = cell(1,90);
for i=1:T2
    temp = [];
    for j=1:nDim
        temp = [temp; Data_NN_inc(i,j:nDim:nDim*T1)];
    end
    Test_Data_Reduced{1,i} = temp';
end

%repetitions for each subject
rt=[9 9 9 9 9 10 8 9 8 10];

%% Calculate RMS within subjects

% Correct sequences
M=0;
rms_train=zeros(1,T2);
for n=1:10
   for i=1:rt(1,n)
       for j=1:rt(1,n)
           for t = 1:T1
               rms_train(M+i) = rms_train(M+i)+norm...
               (Train_Data_Reduced{M+i}(t,:)-Train_Data_Reduced{M+j}(t,:));
           end
       end
        % scale data
        if rt(1,n)>0
           rms_train(M+i)=rms_train(M+i)/rt(1,n);
        else
           rms_train(M+i)=rms_train(M+i);
        end
    end
    M=M+rt(1,n);
end

% Incorrect sequences
N=0;
rms_test=zeros(1,T2);
for n=1:10
    for i=1:rt(1,n)
        for j=1:rt(1,n)
            for t = 1:T1
                rms_test(N+i)= rms_test(N+i)+norm...
                (Test_Data_Reduced{N+i}(t,:)-Train_Data_Reduced{N+j}(t,:));
            end   
        end
        % scale data
        if rt(1,n)>0
           rms_test(N+i)=rms_test(N+i)/rt(1,n);
        else
           rms_test(N+i)=rms_test(N+i);
        end
    end
    N=N+rt(1,n);
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
ylabel('Euclidean Distance', 'fontsize',18);
% Title: Euclidean Distance Within Subjects - Autoencoder Dimensionality Reduction
title('Euclidean Distance WS ENC','fontsize',18);
legend({'Correct Sequences','Incorrect Sequences'}, 'fontsize',16,'location','NW')
set(gca,'box','off','fontweight','bold','LineWidth',2);
set(gcf,'Units','inches','position',[0 0 5.5 4.5]);
% print(h,'../../Results/Euclidean_Distance_WS_ENC','-dpng','-r300');

%% Separation degree betwween the correct and incorrect sequences

SD = zeros(1,1);
for i=1:T2
    for j=1:T2
        SD=SD+(rms_test(i)-rms_train(j))/(abs(rms_test(i))+abs(rms_train(j)));
    end
end
SD=SD/T2/T2;

