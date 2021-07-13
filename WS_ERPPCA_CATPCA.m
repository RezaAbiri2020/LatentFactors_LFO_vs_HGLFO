
% analysis for whole grid
% cross validation of A matrix for LFO and HG

clear all
close all
clc

%% load the related data for the subject

% subject 1

% subject 2
load('/media/reza/WindowsDrive/ECoGLeapMotion/ResultsGroupAnalysis/github_Branch_V3/WS2_KinECoG.mat')

% subject 3


%% Filtering into bands and calculate the power of LFO and HG-LFO

for Fi=1:5
    Signal=ECoG_data(Fi).Finger;
    
    % for lfo
    [b,a]=butter(3,[.5, 4]/(Fs/2));
    LFO(Fi).Finger=abs(hilbert(filtfilt(b,a,Signal)));
    
    % for hg-lfo
    [b,a]=butter(3,[70 150]/(Fs/2));
    F_Filtered=filtfilt(b,a,Signal);
    [F_Envelope,Lower]=envelope(F_Filtered);
    [b,a]=butter(3,[.5, 4]/(Fs/2));
    HG_mean = mean(F_Envelope,1);
    HGLFO(Fi).Finger=filtfilt(b,a,F_Envelope)+repmat(HG_mean,length(F_Envelope),1);
    
end


%% performing general ERP and CAT for trials and signals based on max vel 
% window for analysis
win = 250;

for Fi = 1:5
    
        
    % for flexion
    ERPLFO_trials = [];
    ERPHGLFO_trials = [];
    CATLFO = [];
    CATHGLFO = [];
      
    for j = 1:length(FingersKinInfo.Finger(Fi).Fs_MaxVel)
        
        index = FingersKinInfo.Finger(Fi).Fs_MaxVel(j);
        
        % for ERP
        ERPLFO_trials(:,:,j) = LFO(Fi).Finger(index-win+1:index+win,:);
        ERPHGLFO_trials(:,:,j) = HGLFO(Fi).Finger(index-win+1:index+win,:);
        
        % LFO: randomize the activity within each ch in each trial
        Trial_Data = LFO(Fi).Finger(index-win+1:index+win,:);
        LFO_Shuffled_Trials = [];
        for ch = 1:256
            for shuffle=1:100 % producing shuffled results
                LFO_Shuffled_Trials(j, shuffle,:,ch) = circshift(Trial_Data(:,ch),randi(length(Trial_Data)));
                
            end
            
        end
        
        % HGLFO: randomize the activity within each ch in each trial
        Trial_Data = HGLFO(Fi).Finger(index-win+1:index+win,:);
        HGLFO_Shuffled_Trials = [];
        for ch = 1:256
            for shuffle=1:100 % producing shuffled results
                HGLFO_Shuffled_Trials(j, shuffle,:,ch) = circshift(Trial_Data(:,ch),randi(length(Trial_Data)));
                
            end
            
        end
        
        % for CAT
        CATLFO = [CATLFO ; LFO(Fi).Finger(index-win+1:index+win,:)];
        CATHGLFO = [CATHGLFO ; HGLFO(Fi).Finger(index-win+1:index+win,:)];
        
         
    end
    
    ERP.Fs(Fi).LFO = squeeze(mean(ERPLFO_trials,3));
    ERP.Fs(Fi).HGLFO = squeeze(mean(ERPHGLFO_trials,3));
    ERP_Shuffled.Fs(Fi).LFO = squeeze(mean(LFO_Shuffled_Trials,1)); 
    ERP_Shuffled.Fs(Fi).HGLFO = squeeze(mean(HGLFO_Shuffled_Trials,1)); 
    
    CAT.Fs(Fi).LFO = CATLFO;
    CAT.Fs(Fi).HGLFO = CATHGLFO;
   
    
    % for extension
    ERPLFO_trials = [];
    ERPHGLFO_trials = [];
    CATLFO = [];
    CATHGLFO = [];
    
    for j = 1:length(FingersKinInfo.Finger(Fi).Es_MaxVel)
        
        index = FingersKinInfo.Finger(Fi).Es_MaxVel(j);
        
        % for ERP
        ERPLFO_trials(:,:,j) = LFO(Fi).Finger(index-win+1:index+win,:);
        ERPHGLFO_trials(:,:,j) = HGLFO(Fi).Finger(index-win+1:index+win,:);
        
        % for CAT
        CATLFO = [CATLFO ; LFO(Fi).Finger(index-win+1:index+win,:)];
        CATHGLFO = [CATHGLFO ; HGLFO(Fi).Finger(index-win+1:index+win,:)];
        
        
    end
 
    
    ERP.Es(Fi).LFO = squeeze(mean(ERPLFO_trials,3));
    ERP.Es(Fi).HGLFO = squeeze(mean(ERPHGLFO_trials,3));
    CAT.Es(Fi).LFO = CATLFO;
    CAT.Es(Fi).HGLFO = CATHGLFO;
    
    
end


%% choose the brain area or significant channels for next section

Electrodes = Selected_Chs;
%Electrodes = SigChs_All;
%Electrodes = Selected_HandChs;
%Electrodes = Selected_PMChs;
%Electrodes = Selected_SMChs;
%Electrodes = Selected_PMSMChs;
NumDim=50; % for next PCAs

%% performing general PCA on determined brain areas electrodes or sig channels

for Fi=1 %:5
    % ERP for LFO Flexion
    Data=ERP.Fs(Fi).LFO(:,Electrodes);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    ERPPCA.Fs.LFO(Fi).Coeff=F_coeff;
    ERPPCA.Fs.LFO(Fi).Score=F_score;
    ERPPCA.Fs.LFO(Fi).Explained=F_explained;
    ERPPCA.Fs.LFO(Fi).Variability=F_variability;
    
    figure;
    plot(F_variability,'r');
    xlim([0,NumDim])
    ylim([0,100])
    title(['ERP LFO Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    hold on
    
    % plot the shuffled version
    for shuffle =1:100
        Data = squeeze(ERP_Shuffled.Fs(Fi).LFO(shuffle,:,:));
        Data=Data(:,Electrodes);
        [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
        F_variability=cumsum(F_explained);
        plot(F_variability,'b');
        hold on
        
    end
    
    % CAT for LFO Flexion
    Data=CAT.Fs(Fi).LFO(:,Electrodes);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    CATPCA.Fs.LFO(Fi).Coeff=F_coeff;
    CATPCA.Fs.LFO(Fi).Score=F_score;
    CATPCA.Fs.LFO(Fi).Explained=F_explained;
    CATPCA.Fs.LFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,NumDim])
    ylim([0,100])
    title(['CAT LFO Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
        
    % ERP for HGLFO Flexion
    Data=ERP.Fs(Fi).HGLFO(:,Electrodes);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    ERPPCA.Fs.HGLFO(Fi).Coeff=F_coeff;
    ERPPCA.Fs.HGLFO(Fi).Score=F_score;
    ERPPCA.Fs.HGLFO(Fi).Explained=F_explained;
    ERPPCA.Fs.HGLFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,NumDim])
    ylim([0,100])
    title(['ERP HGLFO Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
    
    % CAT for HGLFO Flexion
    Data=CAT.Fs(Fi).HGLFO(:,Electrodes);
    [F_coeff,F_score,F_latent,F_tsquared,F_explained] = pca(Data);
    F_variability=cumsum(F_explained);
    
    CATPCA.Fs.HGLFO(Fi).Coeff=F_coeff;
    CATPCA.Fs.HGLFO(Fi).Score=F_score;
    CATPCA.Fs.HGLFO(Fi).Explained=F_explained;
    CATPCA.Fs.HGLFO(Fi).Variability=F_variability;
    
    figure;
    stem(F_variability);
    xlim([0,NumDim])
    ylim([0,100])
    title(['CAT HGLFO Finger:',num2str(Fi)]);
    ylabel('% Explained Variation')
    xlabel('Number of PCs')
  
     
end

