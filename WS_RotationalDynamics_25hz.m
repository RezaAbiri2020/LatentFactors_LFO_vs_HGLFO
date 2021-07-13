
% analysis for whole grid
% cross validation of A matrix for LFO and HG within the band and across
% the bands

clear all
close all
clc

%% load the related data for the subject

% subject 1

% subject 2
load('/media/reza/WindowsDrive/ECoGLeapMotion/ResultsGroupAnalysis/github_Branch_V3/WS2_KinECoG_25hz.mat')

% subject 3



%% performing general ERP for trials  

for Fi = 1:5
    HGLFO_trials = [];
    LFO_trials = [];
    
    for j = 1:length(ECoGKin_hglfo(Fi).FingerECoG)
        
        LFO_trials(:,:,j) = ECoGKin_lfo(Fi).FingerECoG(j).Trials;
        HGLFO_trials(:,:,j) = ECoGKin_hglfo(Fi).FingerECoG(j).Trials;     
    end
    
    ERPs.Finger(Fi).LFO = squeeze(mean(LFO_trials,3));
    ERPs.Finger(Fi).HGLFO = squeeze(mean(HGLFO_trials,3));
        
end


%% performing general ERP PCA (without cross validation)on determined brain areas electrodes or sig channels

for Fi=1
    % LFO 
    Data=ERPs.Finger(Fi).LFO;
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    ERPPCA.LFO(Fi).Coeff=F_coeff;
    ERPPCA.LFO(Fi).Score=F_score;
    ERPPCA.LFO(Fi).Explained=F_explained;
    ERPPCA.LFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,min(size(Data,1),size(Data,2))])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    
    % HGLFO 
    Data=ERPs.Finger(Fi).HGLFO;
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    ERPPCA.HGLFO(Fi).Coeff=F_coeff;
    ERPPCA.HGLFO(Fi).Score=F_score;
    ERPPCA.HGLFO(Fi).Explained=F_explained;
    ERPPCA.HGLFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,min(size(Data,1),size(Data,2))])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    
      
end

%% performing CAT for trials  

for Fi = 1:5
    HGLFO_trials = [];
    LFO_trials = [];
    
    for j = 1:length(ECoGKin_hglfo(Fi).FingerECoG)
        
        LFO_trials = [LFO_trials; ECoGKin_lfo(Fi).FingerECoG(j).Trials];
        HGLFO_trials = [HGLFO_trials; ECoGKin_hglfo(Fi).FingerECoG(j).Trials];     
    end
    
    CAT.Finger(Fi).LFO = LFO_trials;
    CAT.Finger(Fi).HGLFO = HGLFO_trials;
        
end

%% performing CAT PCA (without cross validation)on determined brain areas electrodes or sig channels

for Fi=1
    % LFO 
    Data=CAT.Finger(Fi).LFO;
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    CATPCA.LFO(Fi).Coeff=F_coeff;
    CATPCA.LFO(Fi).Score=F_score;
    CATPCA.LFO(Fi).Explained=F_explained;
    CATPCA.LFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,min(size(Data,1),size(Data,2))])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    
    % HGLFO 
    Data=CAT.Finger(Fi).HGLFO;
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    CATPCA.HGLFO(Fi).Coeff=F_coeff;
    CATPCA.HGLFO(Fi).Score=F_score;
    CATPCA.HGLFO(Fi).Explained=F_explained;
    CATPCA.HGLFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,min(size(Data,1),size(Data,2))])
    ylim([0,100])
    title(['Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    
      
end

%% Nik method: rotational dynamics for the first 2PC with cross validation and plots 

% Give the data for cross validating the PCs using A matrix
% input data from single trials:

for Fi = 1 %:5
    % choose an input: LFO or HGLFO    
    % cross validation
    %Input_dynamics=PCs_Trials(:,1:2);
    Input_dynamics=ERPPCA.LFO.Score(:,1:10);
    Output_Dynamics=[zeros(1,10);diff(ERPPCA.LFO.Score(:,1:10))];
    
    for i = 0:1:75-1
        Input_new = circshift(Input_dynamics,i);
        Output_new = circshift(Output_Dynamics,i);
        
        Input_train_OB = Input_new(1:end-1,:);
        Output_Train_OB = Output_new(1:end-1,:);
        
        A=pinv(Input_train_OB'*Input_train_OB)*(Input_train_OB'*Output_Train_OB);
        A_all(:,:,i+1) = A;
        
        Output_train_Pre=Input_train_OB*A;
        Input_train_Pre = Output_Train_OB*pinv(A);
        
        if i<0
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
    % mean of A matrices and coeffs and observing roatinal dynamics
    A_cv = squeeze(mean(A_all,3));
    %save('A_LFO_AllChs.mat','A_cv')
    %save('A_LFO_SigChs.mat','A_cv')
    %save('A_HGDirectLFO_AllChs.mat','A_cv')
    %save('A_HGDirectLFO_SigChs.mat','A_cv')
    %save('A_HGAvgLFO_AllChs.mat','A_cv')
    %save('A_HGAvgLFO_SigChs.mat','A_cv')
          
end


% solving for cross validated phase plane
A = A_cv;
save('AMatrix.mat','A')
tspan=[0,10000];
icond={[0,3]};
figure;
PhasePlane(@system1,tspan,icond,'Xlim',[-10 10],'Ylim',[-10 10],...
    'ArrowHeads',false,'ArrowSize',0.001);

Fi = 1;
% testing all over single trials
for i = 1:length(ECoGKin_lfo(Fi).FingerECoG)
    Trial = ECoGKin_lfo(Fi).FingerECoG(i).Trials;
    
    % projection to PC
    PCs_Observed = Trial*ERPPCA.LFO.Coeff(:,1:2);
    
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
PCs_Observed = ERPPCA.LFO.Score(:,1:10);

% prediction of values
PCs_Predicted = [zeros(1,10);diff(ERPPCA.LFO.Score(:,1:10))]*pinv(A_cv);

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

%% my method1: rotational dynamics for the first 2PC with cross validation and plots 

% Give the data for cross validating the PCs using A matrix
% input data from single trials:
Fi = 1;

% cross validation
Input_dynamics=[ERPPCA.LFO.Score(:,1), ERPPCA.LFO.Score(:,2)];
Output_Dynamics=[[0;diff(ERPPCA.LFO.Score(:,1))], [0;diff(ERPPCA.LFO.Score(:,2))]];
    
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
for i = 1:length(ECoGKin_lfo(Fi).FingerECoG)
    Trial = ECoGKin_lfo(Fi).FingerECoG(i).Trials;
    
    % projection to PC
    PCs_Observed = Trial*ERPPCA.LFO.Coeff(:,1:2);
    
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
PCs_Observed = ERPPCA.LFO.Score(:,1:2);

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
