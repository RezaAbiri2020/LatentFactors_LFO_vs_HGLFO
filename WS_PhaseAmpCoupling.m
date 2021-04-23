
% performing phase amplitude coupling scenarios

% record the significant channels per finger; 
% add up all sig chs for all fingers and record for all other necessary analysis

% see the section 1 and 2 for final results

clear all
close all
clc

%% observing phase-amplitude coupling analysis
% load the related data for the subject

% subject 1
load('E:\ECoGLeapMotion\ResultsGroupAnalysis\github_Branch_V3\Subject1.mat')
% subject 2
load('E:\ECoGLeapMotion\ResultsGroupAnalysis\github_Branch_V3\Subject2.mat')
% subject 3
load('E:\ECoGLeapMotion\ResultsGroupAnalysis\github_Branch_V3\Subject3.mat')



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

%% 1- phase-amplitude coupling for (whole movement)
% how the phase of LFO is coupled to the power(abs-hilbert) of HG during whole movement 

for Fi=1:5
    fprintf(['Finger:',num2str(Fi),'\n'])
    Phase=angle(hilbert(LFO(Fi).Finger));
    Amp=abs(hilbert(HGLFO(Fi).Finger));
    
    for ch=1:length(Selected_Chs)
        if Selected_Chs(ch)
            fprintf(['Ch:',num2str(ch),'\n'])
            %polarplot(Phase(1:end,ch),Amp(1:end,ch),'-')
            PAC_Value=((sum(Amp(1:end,ch).*exp(1i.*Phase(1:end,ch)))))/length(Phase);
            FingerPAC(Fi).ValueCh(ch).RealComplex = PAC_Value;
            
            PhaseCh=Phase(1:end,ch);
            for j=1:2000 % producing shuffled results
                Shuffled_PhaseCh = circshift(PhaseCh,randi(length(PhaseCh)));
                FingerPAC(Fi).ValueCh(ch).Shuffled(j)=(abs(sum(Amp(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Phase);
                
            end
            % z-score per channel
            FingerPAC(Fi).ValueCh(ch).ZScAll=...
                zscore([abs(PAC_Value), FingerPAC(Fi).ValueCh(ch).Shuffled]);
        end
    end
end

%Calculating p-value for all channels
for Fi=1:5
    for ch=1:length(Selected_Chs)
        if Selected_Chs(ch)
            Pvalue=length(find(FingerPAC(Fi).ValueCh(ch).ZScAll(2:end)>...
                FingerPAC(Fi).ValueCh(ch).ZScAll(1)))...
                /length(FingerPAC(Fi).ValueCh(ch).ZScAll);
            FingerPAC(Fi).PValues.ValueCh(ch)=Pvalue;           
        end
    end
end

%% 1-2 saving
save('E:\ECoGLeapMotion\ResultsGroupAnalysis\github_Branch_V3\WS2_PAC.mat','FingerPAC','-v7.3')

%% 1-3 some example plots 
% example plot of observation
Fi=5;
ch=106;
figure;
set(gcf, 'Position', [100, 100, 800, 600]);
hist(FingerPAC(Fi).ValueCh(ch).ZScAll,100)
hold on 
vline(FingerPAC(Fi).ValueCh(ch).ZScAll(1),'color','r')
xlabel('Zsc values of PAC')
ylabel('Repetation of Observations')
title (' Permutation test of PAC; Fi=5; Ch=106')
set(gca,'fontsize',14)
HighQualityFigs('PAC_PA_Permutation')

% example stem plot of p-values for all channels per finger with hline in 0.05
Fi=1;
PValues=[];
for ch=1:length(Selected_Chs)
    PValues=[PValues; FingerPAC(Fi).PValues.ValueCh(ch)];
end
figure(1);
stem(PValues)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
% HighQualityFigs('WholeMovement_PAC_PA')

% example plot of a significant chanel
Fi=5; ch=27;
Phase=angle(hilbert(LFO(Fi).Finger));
Amp=abs(hilbert(HGLFO(Fi).Finger));
figure;
histogram2(Phase(:,ch),Amp(:,ch),150)
xlabel('Phase LFO')
ylabel('Amplitude HG')
zlabel('Number of repetition')
%HighQualityFigs('hist2D_PAC_PA')

%% 1-4 Subject2: showing significant channels on grid per finger & sig channels in brain areas 

SigChs_Fingers=logical(zeros(length(Selected_Chs),5));
for Fi=1:5
    PValues=[];
    for ch=1:length(Selected_Chs)
        if Selected_Chs(ch)
            PValues=[PValues; FingerPAC(Fi).PValues.ValueCh(ch)];
        end
    end
    SigChs_Fingers(find(PValues<=0.05),Fi)=1; %significant
end

% all significant channels across fingers
SigChs_All=zeros(length(Selected_Chs),1);
SigChs_All=sum(SigChs_Fingers,2);

% showing significant channels on brain map per finger and all finger
% toghether
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')

for Fi=1:6
    if Fi<6
        subplot(2,3,Fi)
        ctmr_gauss_plot(cortex,elecmatrix(65:320,:),zeros(length(Selected_Chs),1),'rh'); % rho is a 256ch vector
        el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
        el_add(elecmatrix(find(SigChs_Fingers(:,Fi)==1)+65,:),'msize',5,'color','r');
        title(['Finger',num2str(Fi)]);
    elseif Fi==6
        subplot(2,3,Fi)
        ctmr_gauss_plot(cortex,elecmatrix(65:320,:),zeros(256,1),'rh'); % rho is a 256ch vector
        el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
        for ch = 1:length(Selected_Chs)
            if SigChs_All(ch)>0  
                el_add(elecmatrix(ch+65,:),'msize',2*SigChs_All(ch),'color','r');
            end
        end
        title(['All Fingers']);
        
    end   
end

%HighQualityFigs('All_SigChs_Cleared')

% number of significant ch per brain areas
SigChs_All = logical(SigChs_All);
SigChs_PM = SigChs_All.*Selected_PMChs;
sum(Selected_PMChs)
sum(SigChs_PM)
SigChs_SM = SigChs_All.*Selected_SMChs;
sum(Selected_SMChs)
sum(SigChs_SM)
SigChs_Hand = SigChs_All.*Selected_HandChs;
sum(Selected_HandChs)
sum(SigChs_Hand)