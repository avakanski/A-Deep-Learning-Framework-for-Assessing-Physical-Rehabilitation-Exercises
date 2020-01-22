% M1_Augmentation: generates new instances by adding random noise to the correct instances. 
% The input data is â€œM1-DeepSquatâ€? in the folder â€œData for Distance Functionsâ€?.

close; clear; clc;

% load data 
load('../Data for Distance Functions/M1-DeepSquat-Correct.mat');
load('../Data for Distance Functions/M1-DeepSquat-Incorrect.mat');


%% T1 is the sequence of timesteps and T2 is the length of each sequence
% timesteps
T1 = size(Train_Data{1},1);
% repetitions
T2=length(Train_Data);

%% make noise  [0.01; 0.03; 0.05; 0.07]
mu=0;
sigma =[10^-2; 3*10^-2; 5*10^-2; 7*10^-2]; 
N = cell(4,T2);
for a = 1:4
    for r=1:T2
        N{a,r} = sigma(a) * rand(T1,117) + mu;
        N{a,r} = smoothdata(N{a,r},'SmoothingFactor',1);
    end 
end
h=plot(N{1,1}(:,1));


%% Data Augmentation
M1_Aug = cell(1,6);
M1_Aug{1} = Train_Data;  % correct movement
M1_Aug{6} = Test_Data; % incorrect movment
for r=1:T2
    for a = 1:4
        M1_Aug{a}{r} = Train_Data{r} + N{a,r};
    end
end

% Plot the original data with blue and the augmented data with red lines
for r=1:2
    h = figure;
    plot(M1_Aug{1}{r},'b'), hold on,
    plot(M1_Aug{2}{r},'r');
end
