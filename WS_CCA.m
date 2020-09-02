
% adding CCA results 

clear all
close all
clc

%% load the related data for the subject
load('E:\ECoGLeapMotion\ResultsGroupAnalysis\github_Branch_V3\Subject2.mat')

%% Filtering into bands of LFO and HG-LFO

for Fi=1:5
    Signal=ECoG_data(Fi).Finger;
    
    % for lfo
    [b,a]=butter(3,[.5, 4]/(Fs/2));
    LFO(Fi).Finger=filtfilt(b,a,Signal);
    
    % for hg-lfo
    [b,a]=butter(3,[70 150]/(Fs/2));
    F_Filtered=filtfilt(b,a,Signal);
    [F_Envelope,Lower]=envelope(F_Filtered);
    [b,a]=butter(3,[.5, 4]/(Fs/2));
    HG_mean = mean(F_Envelope,1);
    HGLFO(Fi).Finger=filtfilt(b,a,F_Envelope)+repmat(HG_mean,length(F_Envelope),1);
    
end

%% performing CCA

%%% LOAD THE LFO AND HG LFO MATRICES AS Xa AND Xb %%%
%%% Rows - time-points, Columns - channels
Fi=1;
%Xa_init = LFO_signals(Fi).Hilbert;
%Xb_init = abs(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));

Xa_init = LFO(Fi).Finger;
Xb_init = HGLFO(Fi).Finger;

%Xa = Xa_init(:,Brain_Area);
%Xb = Xb_init(:,Brain_Area);

% for significant Ch
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2\SigChs.mat')
Xa = Xa_init(:,SigChs_All);
Xb = Xb_init(:,SigChs_All);

N = length(Xa);

for iter=1:10
    
    % de-mean
    Xa=Xa-mean(Xa);
    Xb=Xb-mean(Xb);
    
    % parition into training and testing
    trainN = randperm(N,round(0.75*N));
    Xa_train = Xa(trainN,:);
    Xb_train = Xb(trainN,:);
    heldoutN = ones(N,1);
    heldoutN(trainN)=0;
    heldoutN = logical(heldoutN);
    Xa_heldout = Xa(heldoutN,:);
    Xb_heldout = Xb(heldoutN,:);
    
    %%% null - break temporal relationship between Xa and Xb
    
    %%%% randomize Xb_heldout
    %Xb_heldout = Xb_heldout(randperm(size(Xb_heldout,1)),:);
    
    %%%circularly shuffle Xa_heldout
    %     for i=1:size(Xa_heldout,2)
    %         k = randperm(sum(heldoutN),1);
    %         Xa_heldout(:,i) = circshift(Xa_heldout(:,i),k);
    %     end
    
    % demean again
    Xa_train = Xa_train - mean(Xa_train);
    Xb_train = Xb_train - mean(Xb_train);
    % sample covariance matrices
    Caa = ((size(Xa_train,1)-1)^-1)* (Xa_train'*Xa_train);
    Cbb = ((size(Xa_train,1)-1)^-1)* (Xb_train'*Xb_train);
    Cab = ((size(Xa_train,1)-1)^-1)* (Xa_train'*Xb_train);
    
    % cholesky factorization: numerical stability
    Caa12=chol(Caa);
    Cbb12=chol(Cbb);
    
    % solver
    X = (Caa12')\Cab/(Cbb12);
    [U,S,V]=svd(X,0); % Singular values - uncovered canonical correlations
    Wa = Caa12\U; % Each column of Wa gives a weighting of variables of Xa
    Wb = Cbb12\V; % Each column of Wb gives a weighting of variables of Xb
    S=diag(S);
    
    Za_heldout = Xa_heldout*Wa;
    Zb_heldout = Xb_heldout*Wb;
    Sh(iter,:)=diag(corr(Za_heldout,Zb_heldout));
end


Corr2=mean(Sh,1).^2;

figure;
set(gcf, 'Position', [100, 100, 800, 600]);
plot(Corr2)
xlabel('CCA Dimension')
ylabel('Corr^2')
title (' 114 Sig Chs: CV corr of Chs between LFO & HG-LFO')
set(gca,'fontsize',14)
HighQualityFigs('CCA_Corr2_SigChs')

