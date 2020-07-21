

% analysis for HG-Direct-LFO
% analysis for HG-Avg-LFO
% adding brain-area based analysis
% adding analysis using sig channels
% cross validation of A matrix for LFO and HG

clear all
close all
clc
% load signals:
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\LFO_signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Direct_LFO_Signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Avg_LFO_Signals.mat')

% in case of using indexing for flexsion/extension
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\FingersKinIndexes.mat')

%% 1- performing general ERP for trials and signals based on max vel 
% window for analysis
win = 250;

for Fi = 1:5
    HG_Direct = [];
    HG_Avg = [];
    LFO_trials = [];
    HG_Direct_trials = [];
    HG_Avg_trials = [];
    
    HG_Direct = abs(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    HG_Avg = abs(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    
    for j = 1:length(FingersKinIndexes.Finger(Fi).Fs_MaxVel)
        
        index = FingersKinIndexes.Finger(Fi).Fs_MaxVel(j);
        LFO_trials(:,:,j) = LFO_signals(Fi).Hilbert(index-win+1:index+win,:);
        HG_Direct_trials(:,:,j) = HG_Direct(index-win+1:index+win,:);
        HG_Avg_trials(:,:,j) = HG_Avg(index-win+1:index+win,:);
          
    end
    Single_Trials(Fi).Fs(Fi).LFO = LFO_trials;
    Single_Trials(Fi).Fs(Fi).HG_Direct_LFO = HG_Direct_trials;
    Single_Trials(Fi).Fs(Fi).HG_Avg_LFO = HG_Avg_trials;
    
    ERPs.Fs(Fi).LFO = squeeze(mean(LFO_trials,3));
    ERPs.Fs(Fi).HG_Direct_LFO = squeeze(mean(HG_Direct_trials,3));
    ERPs.Fs(Fi).HG_Avg_LFO = squeeze(mean(HG_Avg_trials,3));
    
    LFO_trials = [];
    HG_Direct_trials = [];
    HG_Avg_trials = [];
    
    for j=1:length(FingersKinIndexes.Finger(Fi).Es_MaxVel)
        index = FingersKinIndexes.Finger(Fi).Fs_MaxVel(j);
        LFO_trials(:,:,j) = LFO_signals(Fi).Hilbert(index-win+1:index+win,:);
        HG_Direct_trials(:,:,j) = HG_Direct(index-win+1:index+win,:);
        HG_Avg_trials(:,:,j) = HG_Avg(index-win+1:index+win,:);
        
    end
    
    Single_Trials(Fi).Es(Fi).LFO = LFO_trials;
    Single_Trials(Fi).Es(Fi).HG_Direct_LFO = HG_Direct_trials;
    Single_Trials(Fi).Es(Fi).HG_Avg_LFO = HG_Avg_trials;
    
    ERPs.Es(Fi).LFO = squeeze(mean(LFO_trials,3));
    ERPs.Es(Fi).HG_Direct_LFO = squeeze(mean(HG_Direct_trials,3));
    ERPs.Es(Fi).HG_Avg_LFO = squeeze(mean(HG_Avg_trials,3));
    
end


%% 2- choose the brain area or significant channels for next section
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2\BrainAreas_Logical.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2\SigChs_Logical.mat')

Electrodes = Selected_Chs;
%Electrodes = SigChs_All;
%Electrodes = Selected_HandChs;
%Electrodes = Selected_PMChs;
%Electrodes = Selected_SMChs;
%Electrodes = Selected_PMSMChs;
NumDim=50; % for next PCAs
%% 3- performing general PCA (without cross validation)on determined brain areas electrodes or sig channels

for Fi=1:5
    % LFO Flexion
    Data=ERPs.Fs(Fi).LFO(:,Electrodes);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    PCA.Fs.LFO(Fi).Coeff=F_coeff;
    PCA.Fs.LFO(Fi).Score=F_score;
    PCA.Fs.LFO(Fi).Explained=F_explained;
    PCA.Fs.LFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,NumDim])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    
    % LFO extension
    Data=ERPs.Es(Fi).LFO(:,Electrodes);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    PCA.Es.LFO(Fi).Coeff=F_coeff;
    PCA.Es.LFO(Fi).Score=F_score;
    PCA.Es.LFO(Fi).Explained=F_explained;
    PCA.Es.LFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,NumDim])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    
    % HG_Direct_LFO Flexion
    Data=ERPs.Fs(Fi).HG_Direct_LFO(:,Electrodes);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    PCA.Fs.HG_Direct_LFO(Fi).Coeff=F_coeff;
    PCA.Fs.HG_Direct_LFO(Fi).Score=F_score;
    PCA.Fs.HG_Direct_LFO(Fi).Explained=F_explained;
    PCA.Fs.HG_Direct_LFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,NumDim])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    
    % HG_Direct_LFO extension
    Data=ERPs.Es(Fi).HG_Direct_LFO(:,Electrodes);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    PCA.Es.HG_Direct_LFO(Fi).Coeff=F_coeff;
    PCA.Es.HG_Direct_LFO(Fi).Score=F_score;
    PCA.Es.HG_Direct_LFO(Fi).Explained=F_explained;
    PCA.Es.HG_Direct_LFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,NumDim])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    
    % HG_Avg_LFO Flexion
    Data=ERPs.Fs(Fi).HG_Avg_LFO(:,Electrodes);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    PCA.Fs.HG_Avg_LFO(Fi).Coeff=F_coeff;
    PCA.Fs.HG_Avg_LFO(Fi).Score=F_score;
    PCA.Fs.HG_Avg_LFO(Fi).Explained=F_explained;
    PCA.Fs.HG_Avg_LFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,NumDim])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    
    % HG_Avg_LFO extension
    Data=ERPs.Es(Fi).HG_Avg_LFO(:,Electrodes);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    PCA.Es.HG_Avg_LFO(Fi).Coeff=F_coeff;
    PCA.Es.HG_Avg_LFO(Fi).Score=F_score;
    PCA.Es.HG_Avg_LFO(Fi).Explained=F_explained;
    PCA.Es.HG_Avg_LFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,NumDim])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    
end

% plot the PC weights on grid
close all
%% 4-1: Nik method: rotational dynamics for the first 2PC with cross validation and plots 

% Give the data for cross validating the PCs using A matrix
% input data from single trials:
Fi = 5;
% choose an input: flexion or extention
% 
Trials = Single_Trials(Fi).Fs(Fi).LFO;
ERP_Trials = ERPs.Fs(Fi).LFO(:,Electrodes);
PCs_Trials = PCA.Fs.LFO(Fi).Score;
Coeff_cv = PCA.Fs.LFO(Fi).Coeff;

% Trials = Single_Trials(Fi).Fs(Fi).HG_Direct_LFO;
% ERP_Trials=ERPs.Fs(Fi).HG_Direct_LFO(:,Electrodes);
% PCs_Trials = PCA.Fs.HG_Direct_LFO(Fi).Score;
% Coeff_cv = PCA.Fs.HG_Direct_LFO(Fi).Coeff;

% Trials = Single_Trials(Fi).Fs(Fi).HG_Avg_LFO;
% ERP_Trials=ERPs.Fs(Fi).HG_Avg_LFO(:,Electrodes);
% PCs_Trials = PCA.Fs.HG_Avg_LFO(Fi).Score;
% Coeff_cv = PCA.Fs.HG_Avg_LFO(Fi).Coeff;

% cross validation
%Input_dynamics=PCs_Trials(:,1:2);
Input_dynamics=[PCs_Trials(:,1),PCs_Trials(:,2)];
Output_Dynamics=[[0;diff(PCs_Trials(:,1))], [0;diff(PCs_Trials(:,2))]];

for i = 0:1:2*win-1
    Input_new = circshift(Input_dynamics,i);
    Output_new = circshift(Output_Dynamics,i);
    
    Input_train_OB = Input_new(1:end-1,:);
    Output_Train_OB = Output_new(1:end-1,:);
    
    A=pinv(Input_train_OB'*Input_train_OB)*(Input_train_OB'*Output_Train_OB);
    A_all(:,:,i+1) = A;
    
    Output_train_Pre=Input_train_OB*A;
    Input_train_Pre = Output_Train_OB*pinv(A);
    
    if i<10
        save('AMatrix.mat','A')
        tspan=[0,10000];
        icond={[0,3]};
        figure;
        PhasePlane(@system1,tspan,icond,'Xlim',[-5 5],'Ylim',[-5 5],...
            'ArrowHeads',true,'ArrowSize',1);
        figure;
        subplot(2,1,1)
        plot(Input_train_OB(:,1)); hold on; plot(Input_train_Pre(:,1),'r-','linewidth',1);
        corr1=corr(Input_train_OB(:,1),Input_train_Pre(:,1));
        title(['Train: PC1: ',num2str(corr1^2)])
        subplot(2,1,2)
        plot(Input_train_OB(:,2)); hold on; plot(Input_train_Pre(:,2),'r-','linewidth',1);
        corr2=corr(Input_train_OB(:,2),Input_train_Pre(:,2));
        title(['Train: PC2: ',num2str(corr2^2)])
    end
    
end

close all
% mean of A matrices and coeffs and observing roatinal dynamics
A_cv = squeeze(mean(A_all,3));
%save('A_LFO_AllChs.mat','A_cv')
%save('A_LFO_SigChs.mat','A_cv')
%save('A_HGDirectLFO_AllChs.mat','A_cv')
%save('A_HGDirectLFO_SigChs.mat','A_cv')
%save('A_HGAvgLFO_AllChs.mat','A_cv')
%save('A_HGAvgLFO_SigChs.mat','A_cv')

% solving for cross validated phase plane
A = A_cv;
save('AMatrix.mat','A')
tspan=[0,10000];
icond={[0,3]};
figure;
PhasePlane(@system1,tspan,icond,'Xlim',[-5 5],'Ylim',[-5 5],...
    'ArrowHeads',true,'ArrowSize',1);

% testing all over single trials
for i = 1:size(Trials,3)
    Trial = squeeze(Trials(:,:,i));
    
    % projection to PC
    PCs_Observed = Trial(:,Electrodes)*Coeff_cv(:,1:2);
    
    % prediction of values
    PCs_Predicted = [0, 0; diff(PCs_Observed)]*pinv(A_cv);
     
    % plotting the final train and test data for all trials
    figure;
    set(gcf, 'Position', [100, 100, 800, 600]);
    subplot(2,1,1)
    plot(PCs_Observed(:,1),'b')
    hold on
    plot(PCs_Predicted(:,1),'r')
    legend('Observed','Predicted')
    Corr_Value= corr(PCs_Observed(:,1),PCs_Predicted(:,1));
    PC1_R2(i)=Corr_Value^2;
    title(['Trial: ',num2str(i),' PC1; Cross Validated; R2:',num2str(Corr_Value^2)])
    
    subplot(2,1,2)
    plot(PCs_Observed(:,2),'b')
    hold on
    plot(PCs_Predicted(:,2),'r')
    legend('Observed','Predicted')
    Corr_Value= corr(PCs_Observed(:,2),PCs_Predicted(:,2));
    PC2_R2(i)=Corr_Value^2;
    title(['Trial: ',num2str(i),' PC2; Cross Validated; R2:',num2str(Corr_Value^2)])
    
end

close all
%plot corrs
figure;
set(gcf, 'Position', [100, 100, 800, 600]);
subplot(2,1,1)
plot(PC1_R2)
xlabel('Trial #')
ylabel('Corr^2')
title (' PC1: Cross validated')
set(gca,'fontsize',14)
subplot(2,1,2)
plot(PC2_R2)
xlabel('Trial #')
ylabel('Corr^2')
title (' PC2: Cross validated')
set(gca,'fontsize',14)

% testing all over just ERP trials
% projection to PC
PCs_Observed = [PCs_Trials(:,1),PCs_Trials(:,2)];

% prediction of values
PCs_Predicted = [0,0;diff(PCs_Observed)]*pinv(A_cv);

% plotting the final train and test data for all trials
figure;
set(gcf, 'Position', [100, 100, 800, 600]);
subplot(2,1,1)
plot(PCs_Observed(:,1),'b')
hold on
plot(PCs_Predicted(:,1),'r')
legend('Observed','Predicted')
Corr_Value= corr(PCs_Observed(:,1),PCs_Predicted(:,1));
title(['ERP-Trials; PC1; Cross Validated; R2:',num2str(Corr_Value^2)])

subplot(2,1,2)
plot(PCs_Observed(:,2),'b')
hold on
plot(PCs_Predicted(:,2),'r')
legend('Observed','Predicted')
Corr_Value= corr(PCs_Observed(:,2),PCs_Predicted(:,2));
title(['ERP-Trials; PC2; Cross Validated; R2:',num2str(Corr_Value^2)])

% an example of plot in phase plane
A = A_cv;
save('AMatrix.mat','A')
%tspan=[0 (1:2*win-1)/508.8];
tspan=[0, 500];
icond={[PCs_Predicted(2,1),PCs_Predicted(2,2)]};
figure;
PhasePlane(@system1,tspan,icond,'Xlim',[-2 2],'Ylim',[-2 2],...
    'ArrowHeads',true,'ArrowSize',1);
hold on;
Traj_PC1=PCs_Predicted(2:400,2);
Traj_PC2=PCs_Predicted(2:400,1);
u=[0;diff(Traj_PC1)];
v=[0;diff(Traj_PC2)];
quiver(Traj_PC1(1:10:end),Traj_PC2(1:10:end),u(1:10:end),v(1:10:end),'MaxHeadSize',1,...
    'LineWidth',1.5,'AutoScale','on','color','m')

%% 4-2: my method1: rotational dynamics for the first 2PC with cross validation and plots 

% Give the data for cross validating the PCs using A matrix
% input data from single trials:
Fi = 5;
% choose an input: flexion or extention
Trials = Single_Trials(Fi).Fs(Fi).LFO;
ERP_Trials = ERPs.Fs(Fi).LFO(:,Electrodes);
PCs_Trials = PCA.Fs.LFO(Fi).Score;
Coeff_cv = PCA.Fs.LFO(Fi).Coeff;

%Trials = Single_Trials(Fi).Fs(Fi).HG_Direct_LFO;
%ERP_Trials=ERPs.Fs(Fi).HG_Direct_LFO(:,Electrodes);

%Trials = Single_Trials(Fi).Fs(Fi).HG_Avg_LFO;
%ERP_Trials=ERPs.Fs(Fi).HG_Avg_LFO(:,Electrodes);

% cross validation
Input_dynamics=PCs_Trials(:,1:2);
Output_Dynamics=[[0;diff(PCs_Trials(:,1))], [0;diff(PCs_Trials(:,2))]];

fold_num = 5;
idx_test{1} = 1:floor((1/fold_num)*(length(Input_dynamics)));
idx_train{1} = setdiff(1:length(Input_dynamics),idx_test{1});

for i=2:fold_num
    idx_test{i} = idx_test{i-1}(end)+1:floor((i)*(1/fold_num)*(length(Input_dynamics)));
    idx_train{i}= setdiff(1:length(Input_dynamics),idx_test{i});
end

figure(1);
suptitle(['Folds for PC1 & PC2']); 
for i=1:fold_num
    
    %Input_dynamics_ones=[ones(size(Input_dynamics,1),1),Input_dynamics];
    Input_dynamics_ones = Input_dynamics;
    Input_train_OB = Input_dynamics_ones(idx_train{i},:);
    Input_test_OB = Input_dynamics_ones(idx_test{i},:);
 
    Output_Train_OB=Output_Dynamics(idx_train{i},:);
    Output_Test_OB=Output_Dynamics(idx_test{i},:);
    
    A=pinv(Input_train_OB'*Input_train_OB)*(Input_train_OB'*Output_Train_OB);
    A_all(:,:,i) = A;
    % solving for phase plane
    save('AMatrix.mat','A')
    tspan=[0 100];
    icond={[0, 3]};
    figure;
    PhasePlane(@system1,tspan,icond,'Xlim',[-5 5],'Ylim',[-5 5],...
        'ArrowHeads',true,'ArrowSize',1);
    
    Output_train_Pre=Input_train_OB*A_all(:,:,i);
    Output_test_Pre=Input_test_OB*A_all(:,:,i);
    
    Input_train_Pre = Output_Train_OB*pinv(A_all(:,:,i));
    Input_test_Pre =  Output_Test_OB*pinv(A_all(:,:,i));
    
    figure(1);
    subplot(fold_num,2,i*2-1)
    plot(Input_train_OB); hold on; plot(Input_train_Pre,'r-','linewidth',1);
    corr1=corr(Input_train_OB(:,1),Input_train_Pre(:,1));
    corr2=corr(Input_train_OB(:,2),Input_train_Pre(:,2));
    title(['Train: PC1: ',num2str(corr1^2),' PC2: ',num2str(corr2^2)])
    
    subplot(fold_num,2,i*2)
    plot(Input_test_OB); hold on; plot(Input_test_Pre,'r-','linewidth',1);
    corr1=corr(Input_test_OB(:,1),Input_test_Pre(:,1));
    corr2=corr(Input_test_OB(:,2),Input_test_Pre(:,2));
    title(['Test: PC1: ',num2str(corr1^2),' PC2: ',num2str(corr2^2)])  
end

% mean of A matrices and coeffs and observing roatinal dynamics
A_cv = squeeze(mean(A_all,3));

% solving for cross validated phase plane
A = A_cv;
save('AMatrix.mat','A')
tspan=[0,10000];
icond={[0,3]};
figure;
PhasePlane(@system1,tspan,icond,'Xlim',[-5 5],'Ylim',[-5 5],...
    'ArrowHeads',true,'ArrowSize',1);

% testing all over single trials
for i = 1:size(Trials,3)
    Trial = squeeze(Trials(:,:,i));
    
    % projection to PC
    PCs_Observed = Trial(:,Electrodes)*Coeff_cv(:,1:2);
    
    % prediction of values
    Predicted = pinv(A_cv)*[0, 0; diff(PCs_Observed)]';
    PCs_Predicted = Predicted';
    
    % plotting the final train and test data for all trials
    figure;
    set(gcf, 'Position', [100, 100, 800, 600]);
    subplot(2,1,1)
    plot(PCs_Observed(:,1),'b')
    hold on
    plot(PCs_Predicted(:,1),'r')
    legend('Observed','Predicted')
    Corr_Value= corr(PCs_Observed(:,1),PCs_Predicted(:,1));
    PC1_R2(i)=Corr_Value^2;
    title(['Trial: ',num2str(i),' PC1; Cross Validated; R2:',num2str(Corr_Value^2)])
    
    subplot(2,1,2)
    plot(PCs_Observed(:,2),'b')
    hold on
    plot(PCs_Predicted(:,2),'r')
    legend('Observed','Predicted')
    Corr_Value= corr(PCs_Observed(:,2),PCs_Predicted(:,2));
    PC2_R2(i)=Corr_Value^2;
    title(['Trial: ',num2str(i),' PC2; Cross Validated; R2:',num2str(Corr_Value^2)])
    
end

%plot corrs
figure;
set(gcf, 'Position', [100, 100, 800, 600]);
subplot(2,1,1)
plot(PC1_R2)
xlabel('Trial #')
ylabel('Corr^2')
title (' PC1: Cross validated')
set(gca,'fontsize',14)
subplot(2,1,2)
plot(PC2_R2)
xlabel('Trial #')
ylabel('Corr^2')
title (' PC2: Cross validated')
set(gca,'fontsize',14)

% testing all over just ERP trials
% projection to PC
PCs_Observed = PCs_Trials(:,1:2);

% prediction of values
Predicted = pinv(A_cv)*[0, 0; diff(PCs_Observed)]';
PCs_Predicted = Predicted';

% plotting the final train and test data for all trials
figure;
set(gcf, 'Position', [100, 100, 800, 600]);
subplot(2,1,1)
plot(PCs_Observed(:,1),'b')
hold on
plot(PCs_Predicted(:,1),'r')
legend('Observed','Predicted')
Corr_Value= corrcoef(PCs_Observed(:,1),PCs_Predicted(:,1));
title(['ERP-Trials; PC1; Cross Validated; R2:',num2str(Corr_Value(1,2).^2)])

subplot(2,1,2)
plot(PCs_Observed(:,2),'b')
hold on
plot(PCs_Predicted(:,2),'r')
legend('Observed','Predicted')
Corr_Value= corrcoef(PCs_Observed(:,2),PCs_Predicted(:,2));
title(['ERP-Trials; PC2; Cross Validated; R2:',num2str(Corr_Value(1,2).^2)])


%% 4-3: my method2: rotational dynamics for the first 2PC with cross validation and plots 

% Give the data for cross validating the PCs using A matrix
% input data from single trials:
Fi = 5;
% flexion or extention
%Trials = Single_Trials(Fi).Fs(Fi).LFO;
%ERP_Trials=ERPs.Fs(Fi).LFO(:,Electrodes);

Trials = Single_Trials(Fi).Fs(Fi).HG_Direct_LFO;
ERP_Trials=ERPs.Fs(Fi).HG_Direct_LFO(:,Electrodes);

%Trials = Single_Trials(Fi).Fs(Fi).HG_Avg_LFO;
%ERP_Trials=ERPs.Fs(Fi).HG_Avg_LFO(:,Electrodes);

% load Fs
load('E:\ECoGLeapMotion\DataPatientTwo\ECoGData\ECoG_data.mat');

for cv = 1:5 % number of folds
    N = size(Trials,3);
    
    trainN = zeros(N,1);
    samples = randperm(N,round(0.75*N));
    trainN(samples,1)=1;
    trainN = logical(trainN);
    
    testN = ones(N,1);
    testN(trainN)=0;
    testN = logical(testN);
    
    ERP_train = squeeze(mean(Trials(:,:,trainN),3)); 
    ERP_test =  squeeze(mean(Trials(:,:,testN),3));
    
    % bring to low dimesion for testing
    Data_1 = ERP_train(:,Electrodes);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data_1);
    PC1_train_Observed = F_score(:,1);
    PC2_train_Observed = F_score(:,2);
    
    Coeff1_train = F_coeff(:,1);
    Coeff2_train = F_coeff(:,2);
    Coeff_all(:,:,cv) = [Coeff1_train, Coeff2_train];
    
    Data_2 = ERP_test(:,Electrodes);
    PC1_test_Observed = Data_2*Coeff1_train;
    PC2_test_Observed = Data_2*Coeff2_train;
    
    % generating, testing, plotting and recroding A
    % Xdot=A*X
    X = [PC1_train_Observed, PC2_train_Observed];
    Xdot = diff(X);
    X = X(1:(end-1),:);
    X = X';
    Xdot = Xdot';
    A = Xdot*X'*pinv(X*X');
    A_all(:,:,cv) = A; 
    % solving for phase plane
    save('AMatrix.mat','A')
    tspan=[0 (1:2*win-1)/Fs];
    icond={[PC1_train_Observed(1), PC2_train_Observed(1)]};
    figure;
    PhasePlane(@system1,tspan,icond,'Xlim',[-5 5],'Ylim',[-5 5],...
        'ArrowHeads',true,'ArrowSize',1);
    
    % prediction of values
    if 0
        [t, Predicted]=ode23(@odefun,[0 (1:2*win-1)/Fs],[PC1_train_Observed(1);PC2_train_Observed(1)]);
        PC1_train_Predicted = Predicted(:,1);
        PC2_train_Predicted = Predicted(:,2);
        
        [t, Predicted]=ode23(@odefun,[0 (1:2*win-1)/Fs],[PC1_test_Observed(1);PC2_test_Observed(1)]);
        PC1_test_Predicted = Predicted(:,1);
        PC2_test_Predicted = Predicted(:,2);
    end
    
    Predicted = pinv(A)*[0, 0; diff(PC1_train_Observed), diff(PC2_train_Observed)]';
    Predicted = Predicted';
    PC1_train_Predicted = Predicted(:,1);
    PC2_train_Predicted = Predicted(:,2);
    
    Predicted = pinv(A)*[0, 0; diff(PC1_test_Observed), diff(PC2_test_Observed)]';
    Predicted = Predicted';
    PC1_test_Predicted = Predicted(:,1);
    PC2_test_Predicted = Predicted(:,2);
    
    % plot the training and testing
    figure;
    set(gcf, 'Position', [100, 100, 1000, 800]);
    subplot(2,2,1)
    plot(PC1_train_Observed,'b')
    hold on
    plot(PC1_train_Predicted,'r')
    legend('Observed','Predicted')
    Corr_Value= corrcoef(PC1_train_Observed,PC1_train_Predicted);
    title(['PC1-train; R2:',num2str(Corr_Value(1,2).^2)])
    
    subplot(2,2,2)
    plot(PC1_test_Observed,'b')
    hold on
    plot(PC1_test_Predicted,'r')
    legend('Observed','Predicted')
    Corr_Value= corrcoef(PC1_test_Observed,PC1_test_Predicted);
    title(['PC1-test; R2:',num2str(Corr_Value(1,2).^2)])
    
    subplot(2,2,3)
    plot(PC2_train_Observed,'b')
    hold on
    plot(PC2_train_Predicted,'r')
    legend('Observed','Predicted')
    Corr_Value= corrcoef(PC2_train_Observed,PC2_train_Predicted);
    title(['PC2-train; R2:',num2str(Corr_Value(1,2).^2)])
        
    subplot(2,2,4)
    plot(PC2_test_Observed,'b')
    hold on
    plot(PC2_test_Predicted,'r')
    legend('Observed','Predicted')
    Corr_Value= corrcoef(PC2_test_Observed,PC2_test_Predicted);
    title(['PC2-test; R2:',num2str(Corr_Value(1,2).^2)])
    
end 

% mean of A matrices and coeffs and observing roatinal dynamics
A_cv = squeeze(mean(A_all,3));
Coeff_cv = squeeze(mean(Coeff_all,3));

% solving for cross validated phase plane
A = A_cv;
save('AMatrix.mat','A')
tspan=[0,10000];
icond={[0,3]};
figure;
PhasePlane(@system1,tspan,icond,'Xlim',[-5 5],'Ylim',[-5 5],...
    'ArrowHeads',true,'ArrowSize',1);

% testing all over single trials
for i = 1:size(Trials,3)
    Trial = squeeze(Trials(:,:,i));
    
    % projection to PC
    PCs_Observed = Trial(:,Electrodes)*Coeff_cv;
    
    % prediction of values
    Predicted = pinv(A_cv)*[0, 0; diff(PCs_Observed)]';
    PCs_Predicted = Predicted';
    
    % plotting the final train and test data for all trials
    figure;
    set(gcf, 'Position', [100, 100, 800, 600]);
    subplot(2,1,1)
    plot(PCs_Observed(:,1),'b')
    hold on
    plot(PCs_Predicted(:,1),'r')
    legend('Observed','Predicted')
    Corr_Value= corr(PCs_Observed(:,1),PCs_Predicted(:,1));
    PC1_R2(i)=Corr_Value^2;
    title(['Trial: ',num2str(i),' PC1; Cross Validated; R2:',num2str(Corr_Value^2)])
    
    subplot(2,1,2)
    plot(PCs_Observed(:,2),'b')
    hold on
    plot(PCs_Predicted(:,2),'r')
    legend('Observed','Predicted')
    Corr_Value= corr(PCs_Observed(:,2),PCs_Predicted(:,2));
    PC2_R2(i)=Corr_Value^2;
    title(['Trial: ',num2str(i),' PC2; Cross Validated; R2:',num2str(Corr_Value^2)])
    
end

%plot corrs
figure;
set(gcf, 'Position', [100, 100, 800, 600]);
subplot(2,1,1)
plot(PC1_R2)
xlabel('Trial #')
ylabel('Corr^2')
title (' PC1: Cross validated')
set(gca,'fontsize',14)
subplot(2,1,2)
plot(PC2_R2)
xlabel('Trial #')
ylabel('Corr^2')
title (' PC2: Cross validated')
set(gca,'fontsize',14)


% testing all over just ERP trials
% projection to PC
PCs_Observed = ERP_Trials*Coeff_cv;

% prediction of values
Predicted = pinv(A_cv)*[0, 0; diff(PCs_Observed)]';
PCs_Predicted = Predicted';

% plotting the final train and test data for all trials
figure;
set(gcf, 'Position', [100, 100, 800, 600]);
subplot(2,1,1)
plot(PCs_Observed(:,1),'b')
hold on
plot(PCs_Predicted(:,1),'r')
legend('Observed','Predicted')
Corr_Value= corrcoef(PCs_Observed(:,1),PCs_Predicted(:,1));
title(['ERP-Trials; PC1; Cross Validated; R2:',num2str(Corr_Value(1,2).^2)])

subplot(2,1,2)
plot(PCs_Observed(:,2),'b')
hold on
plot(PCs_Predicted(:,2),'r')
legend('Observed','Predicted')
Corr_Value= corrcoef(PCs_Observed(:,2),PCs_Predicted(:,2));
title(['ERP-Trials; PC2; Cross Validated; R2:',num2str(Corr_Value(1,2).^2)])

%% checking sequence peacking using ERP real and predicted ERP within a band LFO or HG

% ERP real --->  PCs_Observed ---> PCs_Predicted ---> predicted ERP? 

% PCs_Predicted = (predicted ERP) * Coeff_cv ---> predicted ERP?

% 500*1 = (500* 241) * (241*1)

% for having inverse and better accuracy: 500*2 = (500* 241) * (241*2) 
% projection to PC
PCs_Observed = [PCs_Trials(:,1),PCs_Trials(:,2)];
% prediction of values
PCs_Predicted = [0,0;diff(PCs_Observed)]*pinv(A_cv);

Coeff_ERP = Coeff_cv(:,1:2);
ERP_Predicted = PCs_Predicted*pinv(Coeff_ERP);
ERP_Observed = ERP_Trials;

% plot and compare the sequence peacking between the two ERP
data1=ERP_Predicted;
data2=zeros(size(data1));
for i=1:size(ERP_Observed,2)
    data2(:,i)=normalize(data1(:,i),'range');
end
data3=data2';
Rows=[];
for i=1:2*win
    [m,n]=max(data3(:,i));
    if isempty(find(Rows==n))
        Rows=[Rows;n];
    end
    
end
%save('Rows_LFO_AllChs_Pre.mat','Rows')
OrderedData=data3(Rows,:);
figure;
set(gcf, 'Position', [300, 100, 700,850]);
suptitle(['LFO ERP-Pre Ordered by OB; Finger: ',num2str(Fi)]);
imagesc(OrderedData)
colorbar
yticks('')
yticks(1:3:length(Rows))
ytickangle(0)
yticklabels(num2str(Rows(1:3:end)))


% look at the single trials for any patterns
load('Rows_LFO_AllChs_OB.mat')

for i = 1:size(Trials,3)
    Trial = squeeze(Trials(:,:,i));
    data2 = Trial(:,Electrodes);
    data3 = data2';
    OrderedData=data3(Rows,:);
    figure;
    set(gcf, 'Position', [300, 100, 700,850]);
    suptitle(['LFO SingleTrial Ordered by OB; Finger: ',num2str(Fi)]);
    imagesc(OrderedData)
    colorbar
    yticks('')
    yticks(1:3:length(Rows))
    ytickangle(0)
    yticklabels(num2str(Rows(1:3:end))) 
    
end

%% checking sequence peacking using ERP real and predicted ERP across bands of LFO and HG

% ERP real --->  PCs_Observed ---> PCs_Predicted ---> predicted ERP? 

% PCs_Predicted = (predicted ERP) * Coeff_cv ---> predicted ERP?

% 500*1 = (500* 241) * (241*1)

% for having inverse and better accuracy: 500*2 = (500* 241) * (241*2) 
% projection to PC
PCs_Observed = [PCs_Trials(:,1),PCs_Trials(:,2)];
% prediction of values
load('A_HGDirectLFO_AllChs.mat')
PCs_Predicted = [0,0;diff(PCs_Observed)]*pinv(A_cv);

Coeff_ERP = Coeff_cv(:,1:2);
ERP_Predicted = PCs_Predicted*pinv(Coeff_ERP);
ERP_Observed = ERP_Trials;

% plot and compare the sequence peacking between the two ERP
data1=ERP_Predicted;
data2=zeros(size(data1));
for i=1:size(ERP_Observed,2)
    data2(:,i)=normalize(data1(:,i),'range');
end
data3=data2';
Rows=[];
for i=1:2*win
    [m,n]=max(data3(:,i));
    if isempty(find(Rows==n))
        Rows=[Rows;n];
    end
    
end
OrderedData=data3(Rows,:);
%save('Rows_LFO_AllChs_PreByHG.mat','Rows')
figure;
set(gcf, 'Position', [300, 100, 700,850]);
suptitle(['Ordered Ch; LFO ERP-Pre by HG-LFO; Finger: ',num2str(Fi)]);
imagesc(OrderedData)
colorbar
yticks('')
yticks(1:1:length(Rows))
ytickangle(0)
yticklabels(num2str(Rows(1:1:end)))


% order the prediction of LFO by the order of prediction using HG-LFO
PCs_Observed = [PCs_Trials(:,1),PCs_Trials(:,2)];
% prediction of values
load('A_LFO_AllChs.mat')
PCs_Predicted = [0,0;diff(PCs_Observed)]*pinv(A_cv);

Coeff_ERP = Coeff_cv(:,1:2);
ERP_Predicted = PCs_Predicted*pinv(Coeff_ERP);
ERP_Observed = ERP_Trials;

data3=ERP_Predicted';
load('Rows_LFO_AllChs_PreByHG.mat')
OrderedData=data3(Rows,:);
%save('Rows_LFO_AllChs_PreByHG.mat','Rows')
figure;
set(gcf, 'Position', [300, 100, 700,850]);
suptitle(['LFO ERP-Pre Ordered by (ERP-Pre by HG-LFO); Finger: ',num2str(Fi)]);
imagesc(OrderedData)
colorbar
yticks('')
yticks(1:1:length(Rows))
ytickangle(0)
yticklabels(num2str(Rows(1:1:end)))

