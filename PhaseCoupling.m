
% for this commit
% performing different phase amplitude coupling scenarios

% record the significant channels per finger; 
% add up all sig chs for all fingers and record for all other necessary analysis

% see the section 1 and 2 for final results


clear all
close all
clc
% observing different coupling analysis for these signals:
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\LFO_signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Direct_LFO_Signals.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\HG_Avg_LFO_Signals.mat')

% using indexing for flexsion/extension
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1\FingersKinIndexes.mat')

%% 1- phase-amplitude coupling for (whole movement)

for Fi=1:5
    fprintf(['Finger:',num2str(Fi),'\n'])
    Phase=angle(hilbert(LFO_signals(Fi).Filtered));
    Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).PureEnv));
    Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).PureEnv));
    
    for ch=1:256
        fprintf(['Ch:',num2str(ch),'\n'])
        %polarplot(Phase(1:end,ch),Amp(1:end,ch),'-')
        PAC_Value1=((sum(Amp_Direct(1:end,ch).*exp(1i.*Phase(1:end,ch)))))/length(Phase);
        PAC_Value2=((sum(Amp_Avg(1:end,ch).*exp(1i.*Phase(1:end,ch)))))/length(Phase);
        FingerPAC(Fi).WholeHG.ValueCh(ch).RealComplex = PAC_Value1;
        FingerPAC(Fi).AvgHG.ValueCh(ch).RealComplex = PAC_Value2;
        
        PhaseCh=Phase(1:end,ch);
        for j=1:2000 % producing shuffled results
            Shuffled_PhaseCh = circshift(PhaseCh,randi(length(PhaseCh)));
            FingerPAC(Fi).WholeHG.ValueCh(ch).Shuffled(j)=(abs(sum(Amp_Direct(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Phase);
            FingerPAC(Fi).AvgHG.ValueCh(ch).Shuffled(j)=(abs(sum(Amp_Avg(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Phase);
            
        end
        % z-score per channel
        FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll=...
            zscore([abs(PAC_Value1), FingerPAC(Fi).WholeHG.ValueCh(ch).Shuffled]);
        FingerPAC(Fi).AvgHG.ValueCh(ch).ZScAll=...
            zscore([abs(PAC_Value2), FingerPAC(Fi).AvgHG.ValueCh(ch).Shuffled]);
    end
end

% plot of observations; example
figure
hist(FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll,100)
set(gca,'fontsize',16)
figure
hist(FingerPAC(Fi).AvgHG.ValueCh(ch).ZScAll,100)
set(gca,'fontsize',16)

%Calculating p-value for all channels
for Fi=1:5
    for ch=1:256
        Pvalue1=length(find(FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll(2:end)>...
            FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll(1)))...
            /length(FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll);
        FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)=Pvalue1;
        
        Pvalue2=length(find(FingerPAC(Fi).AvgHG.ValueCh(ch).ZScAll(2:end)>...
            FingerPAC(Fi).AvgHG.ValueCh(ch).ZScAll(1)))...
            /length(FingerPAC(Fi).AvgHG.ValueCh(ch).ZScAll);
        FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)=Pvalue2;
        
    end
end

%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2/FingerPAC_PA_WholeMove.mat','FingerPAC','-v7.3')

% stem plots of p-values for all channels per finger with hline in 0.05
Fi=1;
% for flexion extension
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)];
    PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)];
end
figure(1);
subplot(2,1,1)
stem(PValuesDirect)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('HG-Direct-LFO')
subplot(2,1,2)
stem(PValuesAvg)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('HG-Avg-LFO')

HighQualityFigs('WholeMovement_PAC_PA')

% plot of example of significant chanel
Fi=5; ch=27;
Phase=angle(hilbert(LFO_signals(Fi).Filtered));
Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).PureEnv));
Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).PureEnv));
figure;
histogram2(Phase(:,ch),Amp_Direct(:,ch),150)
xlabel('Phase LFO')
ylabel('Amplitude HG')
zlabel('Number of repetition')
HighQualityFigs('hist2D_PAC_PA')

% showing significant channels on brain map per finger
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')
Sig_Fis = zeros(256,5);
% for direct HG
for Fi=1:5
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)];
    end
    SigChs_Fingers=find(PValuesDirect<0.05); %significant
    Sig_Fis(SigChs_Fingers,Fi) =1;
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesDirect),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs_Fingers+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end

Sig_Chs = sum(Sig_Fis');
Sig_Chs = Sig_Chs';
subplot(2,3,6)
ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*Sig_Chs),'rh'); % rho is a 256ch vector
el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
for ch = 1:256
    if Sig_Chs(ch)>0
        el_add(elecmatrix(ch+65,:),'msize',Sig_Chs(ch)*1.5,'color','r');
    end
end
title(['All Fingers']);

HighQualityFigs('WholeMove_PAC_PA_Grid_Direct_2')
   

% for avg HG
for Fi=1:5
    PValuesAvg=[];
    for ch=1:256
        PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)];
    end
    SigChs_Fingers=find(PValuesAvg<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesAvg),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs_Fingers+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end
HighQualityFigs('WholeMove_PAC_PA_Grid_Avg')

%% phase-amplitude coupling for (whole movement)
% how the phase of LFO is coupled to the phase of
% HG-Direct-LFO & HG-Avg-LFO during whole movement of one finger
for Fi=1%:5
    fprintf(['Finger:',num2str(Fi),'\n'])
    Phase=angle(hilbert(LFO_signals(Fi).Filtered));
    Phase_Direct=angle(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    Phase_Avg=angle(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    
    for ch=1:256
        fprintf(['Ch:',num2str(ch),'\n'])
        
        PAC_Value1=(sum(exp(1i.*(Phase(1:end,ch)-Phase_Direct(1:end,ch)))))/length(Phase);
        PAC_Value2=(sum(exp(1i.*(Phase_Avg(1:end,ch)-Phase(1:end,ch)))))/length(Phase);
        
        FingerPAC(Fi).WholeHG.ValueCh(ch).RealComplex=PAC_Value1;
        FingerPAC(Fi).AvgHG.ValueCh(ch).RealComplex=PAC_Value2;
        
        PhaseCh=Phase(1:end,ch);
        PhaseCh2=Phase_Direct(1:end,ch);
        for j=1:1000 % producing shuffled results
            Shuffled_PhaseCh = PhaseCh(randperm(numel(PhaseCh)));
            %Shuffled_PhaseCh = circshift(PhaseCh,randi(length(PhaseCh)));
            Shuffled_Phase_Direct = PhaseCh2(randperm(numel(PhaseCh2)));
            %Shuffled_Phase_Direct = circshift(Phase_Direct(1:end,ch),randi(length(Phase_Direct(1:end,ch))));
            FingerPAC(Fi).WholeHG.ValueCh(ch).Shuffled(j)=...
                (sum(exp(1i.*(Phase(1:end,ch)-Shuffled_Phase_Direct(1:end)))))/length(Phase);
            FingerPAC(Fi).AvgHG.ValueCh(ch).Shuffled(j)=...
                (sum(exp(1i.*(Phase_Avg(1:end,ch)-Shuffled_PhaseCh(1:end)))))/length(Phase);
        end
    end
    
end

%Calculating p-value for all channels
for Fi=1%:5
    for ch=1:256
        RealValue1 = FingerPAC(Fi).WholeHG.ValueCh(ch).RealComplex;
        ShuffledValue1 = FingerPAC(Fi).WholeHG.ValueCh(ch).Shuffled;
        
        Pvalue1=length(find(abs((RealValue1))>abs((ShuffledValue1))))/...
            (length(FingerPAC(Fi).WholeHG.ValueCh(ch).Shuffled)+1);
        FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)=Pvalue1;
        
        RealValue2 = FingerPAC(Fi).AvgHG.ValueCh(ch).RealComplex;
        ShuffledValue2 = FingerPAC(Fi).AvgHG.ValueCh(ch).Shuffled;
        
        Pvalue2=length(find(abs((RealValue2))>abs((ShuffledValue2))))/...
            (length(FingerPAC(Fi).AvgHG.ValueCh(ch).Shuffled)+1);
        FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)=Pvalue2;
            
    end
end

save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2/FingerPAC_PP_WholeMove.mat','FingerPAC','-v7.3')

% stem plots of p-values for all channels per finger with hline in 0.05
Fi=1;
% for flexion extension
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)];
    PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)];
end
figure(1);
subplot(2,1,1)
stem(PValuesDirect)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('HG-Direct-LFO')
subplot(2,1,2)
stem(PValuesAvg)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('HG-Avg-LFO')

HighQualityFigs('WholeMovement_PAC_PP')

% plot of example of significant chanel
Fi=5; ch=27;
Phase=angle(hilbert(LFO_signals(Fi).Filtered));
Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).PureEnv));
Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).PureEnv));
figure;
histogram2(Phase(:,ch),Amp_Direct(:,ch),150)
xlabel('Phase LFO')
ylabel('Amplitude HG')
zlabel('Number of repetition')
HighQualityFigs('hist2D_PAC_PP')

% showing significant channels on brain map per finger
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')

% for direct HG
for Fi=1:5
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)];
    end
    SigChs_Fingers=find(PValuesDirect<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesDirect),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs_Fingers+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end
HighQualityFigs('WholeMove_PAC_PP_Grid_Direct')

% for avg HG
for Fi=1:5
    PValuesAvg=[];
    for ch=1:256
        PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)];
    end
    SigChs_Fingers=find(PValuesAvg<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesAvg),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs_Fingers+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end
HighQualityFigs('WholeMove_PAC_PP_Grid_Avg')

%% phase-amplitude coupling for considering a flexion and a extension as one trial
% how the phase of LFO is coupled to the power(abs-hilbert) of HG during flexions or extensions
% it should be the envelope of HG-Direct-LFO or HG-Avg-LFO 

for Fi=1:5
    fprintf(['Finger:',num2str(Fi),'\n'])
    Phase=angle(hilbert(LFO_signals(Fi).Filtered));
    Amp_Direct=abs(hilbert(HG_Direct_LFO_Signals(Fi).PureEnv));
    Amp_Avg=abs(hilbert(HG_Avg_LFO_Signals(Fi).PureEnv));
    
    Index_Fi=sort([1; FingersKinIndexes.Finger(Fi).Fs_MaxPos;...
        FingersKinIndexes.Finger(Fi).Es_MaxPos; length(Phase)]);
    
    % for flexions and extension
    kk=0;
    for k=1:2:length(Index_Fi)-2
        fprintf(['Flexion & Extension:',num2str(k),'\n'])
        Phase_k=Phase(Index_Fi(k):Index_Fi(k+2),:);
        Amp_Direct_k=Amp_Direct(Index_Fi(k):Index_Fi(k+2),:);
        Amp_Avg_k=Amp_Avg(Index_Fi(k):Index_Fi(k+2),:);
        kk=kk+1;
                
        for ch=1:256
            %fprintf(['Ch:',num2str(ch),'\n'])
            %polarplot(Phase(1:end,ch),Amp_Direct(1:end,ch),'-')
            %polarplot(Phase(1:end,ch),Amp_Avg(1:end,ch),'-')
            PAC_Value1=((sum(Amp_Direct_k(1:end,ch).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
            PAC_Value2=((sum(Amp_Avg_k(1:end,ch).*exp(1i.*Phase_k(1:end,ch)))))/length(Phase_k);
            
            FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).RealComplex = PAC_Value1;
            FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).RealComplex = PAC_Value2;
            
            PhaseCh=Phase_k(1:end,ch);
            for j=1:1000 % producing shuffled results
                Shuffled_PhaseCh = circshift(PhaseCh,randi(length(PhaseCh)));
                FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).Shuffled(j)=(abs(sum(Amp_Direct_k(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Phase_k);
                FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).Shuffled(j)=(abs(sum(Amp_Avg_k(1:end,ch).*exp(1i.*Shuffled_PhaseCh(1:end)))))/length(Phase_k);
                
            end
            % z-score per channel
            FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).ZScAll=...
                zscore([abs(PAC_Value1), FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).Shuffled]);
            FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).ZScAll=...
                zscore([abs(PAC_Value2), FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).Shuffled]);     
        end
    end   
end

% plot of observations; example
figure;
histogram2(Phase_k(:,ch),Amp_Direct_k(:,ch))
figure;
histogram2(Phase_k(:,ch),Amp_Avg_k(:,ch))

%Calculating p-value for all channels across trials
for Fi=1:5
    % for flexion and extension
    for ch=1:256
        RealValue1=[];
        ShuffledValue1=[];
        RealValue2=[];
        ShuffledValue2=[];
        for trial=1:length(FingerPAC(Fi).FE)
            RealValue1=[RealValue1; FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).RealComplex];
            ShuffledValue1=[ShuffledValue1; FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).Shuffled];
            RealValue2=[RealValue2; FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).RealComplex];
            ShuffledValue2=[ShuffledValue2; FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).Shuffled];
        end
        Pvalue1=length(find(abs(mean(RealValue1))>abs(mean(ShuffledValue1))))/...
            (length(FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).Shuffled)+1);
        FingerPAC(Fi).PValuesFE.WholeHG.ValueCh(ch)=Pvalue1;
        
        Pvalue2=length(find(abs(mean(RealValue2))>abs(mean(ShuffledValue2))))/...
            (length(FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).Shuffled)+1);
        FingerPAC(Fi).PValuesFE.AvgHG.ValueCh(ch)=Pvalue2;
        
    end
     
end

%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V1/FingerPAC.mat','FingerPAC','-v7.3')
% stem plots of p-values for all channels per finger with hline in 0.05

Fi=5;
% for flexion extension
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValuesFE.WholeHG.ValueCh(ch)];
    PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValuesFE.AvgHG.ValueCh(ch)];
end
figure(1);
subplot(2,1,1)
stem(PValuesDirect)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Flexion-Extension; HG-Direct-LFO')
subplot(2,1,2)
stem(PValuesAvg)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Flexion-Extension; HG-Avg-LFO')

HighQualityFigs('FE_PAC')


% showing significant channels on brain map per finger

load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')

% for flexion-extension; for direct HG
for Fi=1:5
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValuesFE.WholeHG.ValueCh(ch)];
    end
    SigChs_Fingers=find(PValuesDirect<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesDirect),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs_Fingers+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end
HighQualityFigs('FE_PAC_Grid_Direct')

% for flexion-extension; for Avg HG
for Fi=1:5
    PValuesAvg=[];
    for ch=1:256
        PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValuesFE.AvgHG.ValueCh(ch)];
    end
    SigChs_Fingers=find(PValuesAvg<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesAvg),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs_Fingers+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end
HighQualityFigs('FE_PAC_Grid_Avg')


%% phase-amplitude coupling for considering a flexion and a extension as one trial
% how the phase of LFO is coupled to the phase of
% HG-Direct-LFO & HG-Avg-LFO during flexions and extensions

for Fi=1:5
    fprintf(['Finger:',num2str(Fi),'\n'])
    Phase=angle(hilbert(LFO_signals(Fi).Filtered));
    Phase_Direct=angle(hilbert(HG_Direct_LFO_Signals(Fi).DeltaofEnv));
    Phase_Avg=angle(hilbert(HG_Avg_LFO_Signals(Fi).DeltaofEnv));
    
    Index_Fi=sort([1; FingersKinIndexes.Finger(Fi).Fs_MaxPos;...
        FingersKinIndexes.Finger(Fi).Es_MaxPos; length(Phase)]);
    
    % for flexions and extension
    kk=0;
    for k=1:2:length(Index_Fi)-2
        fprintf(['Flexion & Extension:',num2str(k),'\n'])
        Phase_k=Phase(Index_Fi(k):Index_Fi(k+2),:);
        Phase_Direct_k=Phase_Direct(Index_Fi(k):Index_Fi(k+2),:);
        Phase_Avg_k=Phase_Avg(Index_Fi(k):Index_Fi(k+2),:);
        kk=kk+1;
        
        for ch=1:256
            PAC_Value1=(sum(exp(1i.*(Phase_k(1:end,ch)-Phase_Direct_k(1:end,ch)))))/length(Phase_k);
            PAC_Value2=(sum(exp(1i.*(Phase_k(1:end,ch)-Phase_Avg_k(1:end,ch)))))/length(Phase_k);
            
            FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).RealComplex=PAC_Value1;
            FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).RealComplex=PAC_Value2;
            
            PhaseCh=Phase_k(1:end,ch);
            for j=1:1000 % producing shuffled results
                Shuffled_PhaseCh = PhaseCh(randperm(numel(PhaseCh)));
                %Shuffled_PhaseCh = circshift(PhaseCh,randi(length(PhaseCh)));
                
                FingerPAC(Fi).FE(kk).WholeHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Shuffled_PhaseCh(1:end)-Phase_Direct_k(1:end,ch)))))/length(Phase_k);                
                FingerPAC(Fi).FE(kk).AvgHG.ValueCh(ch).Shuffled(j)=(sum(exp(1i.*(Shuffled_PhaseCh(1:end)-Phase_Avg_k(1:end,ch)))))/length(Phase_k);                
                
            end
           
        end
    end
          
end

%Calculating p-value for all channels across trials and fingers

for ch=1:256
    fprintf(['Ch:',num2str(ch),'\n'])
    RealValue1=[];
    ShuffledValue1=[];
    RealValue2=[];
    ShuffledValue2=[];
    for Fi=1:5
        % for flexion and extension
        for trial=1:length(FingerPAC(Fi).FE)
            RealValue1=[RealValue1; FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).RealComplex];
            ShuffledValue1=[ShuffledValue1; FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).Shuffled];
            RealValue2=[RealValue2; FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).RealComplex];
            ShuffledValue2=[ShuffledValue2; FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).Shuffled];
        end
    end
    Pvalue1=length(find(abs(mean(RealValue1))>abs(mean(ShuffledValue1))))/...
        (length(FingerPAC(Fi).FE(trial).WholeHG.ValueCh(ch).Shuffled)+1);
    FingerPAC(1).PValuesFE.WholeHG.ValueCh(ch)=Pvalue1;
    
    Pvalue2=length(find(abs(mean(RealValue2))>abs(mean(ShuffledValue2))))/...
        (length(FingerPAC(Fi).FE(trial).AvgHG.ValueCh(ch).Shuffled)+1);
    FingerPAC(1).PValuesFE.AvgHG.ValueCh(ch)=Pvalue2;
    
end

%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2/FingerPAC.mat','FingerPAC','-v7.3')

% stem plots of p-values for all channels with hline in 0.05
Fi=1;
% for flexion extension
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValuesFE.WholeHG.ValueCh(ch)];
    PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValuesFE.AvgHG.ValueCh(ch)];
end
figure(1);
subplot(2,1,1)
stem(PValuesDirect)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Flexion-Extension; HG-Direct-LFO')
subplot(2,1,2)
stem(PValuesAvg)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('Flexion-Extension; HG-Avg-LFO')

HighQualityFigs('FE_PAC')


if 0
% compass plots of coupled-phase per channel for all finger with 
% significant channel with red lines
Nch_Record=256;

% grid layout
Ch_num_1=1:Nch_Record;
Ch_num_2=reshape(Ch_num_1,[16,16]);
% Ch_num_3 is the grid layout
Ch_num_3=rot90(rot90(Ch_num_2));

for Fi=1:5
    figure;
    set(gcf, 'Position', [100, 100, 1000, 800]);
    suptitle(['coupled-angle per channel; Finger: ',num2str(Fi)]);
    
    for i=1:256
        subplot(16,16,i)
        r=ceil(i/16);
        c=mod(i,16);
        if c==0
            c=8;
        end
        if FingerPAC(Fi).ValueCh(Ch_num_3(r,c)).PValue<0.05
            polarplot([0 FingerPAC(Fi).ValueCh(Ch_num_3(r,c)).RealAngle],[0 1],'-r','linewidth',1.5)
            rticks('')
            thetaticks('')
        else
            polarplot([0 FingerPAC(Fi).ValueCh(Ch_num_3(r,c)).RealAngle],[0 1],'-b','linewidth',1.5)
            rticks('')
            thetaticks('')
        end
        
        hold on
          
    end
end
end

% showing significant channels on brain map per finger
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')

% for flexion-extension; for direct HG
for Fi=1
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValuesFE.WholeHG.ValueCh(ch)];
    end
    SigChs_Fingers=find(PValuesDirect<0.05); %significant
    % PValuesDirect(find(PValuesDirect>=0.05))=0; % not significant
    figure;
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesDirect),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs_Fingers+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar     
end

% for flexion-extension; for Avg HG
for Fi=1
    PValuesAvg=[];
    for ch=1:256
        PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValuesFE.AvgHG.ValueCh(ch)];
    end
    SigChs_Fingers=find(PValuesAvg<0.05); %significant
    %PValuesAvg(find(PValuesAvg>=0.05))=0; % not significant
    figure;
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(0*PValuesAvg),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    el_add(elecmatrix(SigChs_Fingers+65,:),'msize',5,'color','r');
    title(['Finger',num2str(Fi)]);
    %colorbar  
end


%% 2-Performing more analysis on recorded significant channels 

%% 2-1: loading information for bad channels and also brain-area based channels

% loading and breaking raw ECoG data into trials
load('E:\ECoGLeapMotion\DataPatientTwo\ECoGData\ECoG_data.mat');

% Final list of bad channels
BadChs = [65 66 128 128+1 128+2 128+16 128+18 128+20 128+24 128+30 ...
    128+31 128+32 128+44 128+62 128+64];

% make channels Ready for median or mean
All_Index = ones(size(lfp,1),1);
All_Index(BadChs') = 0;

% only 256 Channels were used for real recording
Nch_Record = 256;
Selected_Chs = logical(All_Index(1:Nch_Record,1));

% time markers in ecog
[Nch,N] = size(lfp);
ecog_time = (0:N-1)/Fs;
anin = anin(1,1:N);

% from Daniel function: start and end of movements
[trial_start_time,trial_end_time] = get_trial_times(anin,Fs);

trial_start = trial_start_time*Fs;
trial_stop = trial_end_time*Fs;

% Finding the precentral (for Primary Motor=PM) and postcentral (somatosensory=SM) electrodes/channels
 load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')
 
 PM_chs=[];
 SM_chs=[];
 PMSM_chs=[];
 
 for i=65:320
     
     % removing/passing bad channels 
     if ismember(i-64,BadChs)
         continue 
         
      % finding PM channels   
     elseif strcmp(anatomy{i,4},'precentral')
         PM_chs=[PM_chs;(i-64)];
         PMSM_chs=[PMSM_chs;(i-64)];
         
      % finding SM channels   
     elseif strcmp(anatomy{i,4},'postcentral')
         SM_chs=[SM_chs;(i-64)];
         PMSM_chs=[PMSM_chs;(i-64)];
         
     end   
     
     
 end
 
% logic PM channels
All_Index_PM=zeros(size(lfp,1),1);
All_Index_PM(PM_chs)=1;
Selected_PMChs=logical(All_Index_PM(1:Nch_Record,1));
 

% logic SM channels
All_Index_SM=zeros(size(lfp,1),1);
All_Index_SM(SM_chs)=1;
Selected_SMChs=logical(All_Index_SM(1:Nch_Record,1));

% logic PMSM channels
All_Index_PMSM=zeros(size(lfp,1),1);
All_Index_PMSM(PMSM_chs)=1;
Selected_PMSMChs=logical(All_Index_PMSM(1:Nch_Record,1));


% Finding the hand area electrodes/channels with observation
% hand area channels
for i=0:4
    HandChs(((i*7)+1):((i*7)+7))=(122:128)+i*16;
end 

 All_Index_Hand=zeros(Nch_Record,1);
 for i=1:Nch_Record
     
     % removing/passing bad channels 
     if ismember(i,BadChs)
         continue 
         
      % finding hand channels   
     elseif ismember(i,HandChs)
         All_Index_Hand(i,1)=1;
         

         
     end   
     
     
 end
 
 % logic Hand channels
Selected_HandChs=logical(All_Index_Hand);


%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2/BrainAreas_Logical.mat','Selected_Chs','Selected_PMChs','Selected_SMChs','Selected_PMSMChs','Selected_HandChs')


%% 2-2: start of analysis for sig channels
% load the correct data
load('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2\FingerPAC_PA_WholeMove.mat')
Fi=5;
ch=106;
% plot of observations; example
figure;
set(gcf, 'Position', [100, 100, 800, 600]);
hist(FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll,100)
hold on 
vline(FingerPAC(Fi).WholeHG.ValueCh(ch).ZScAll(1),'color','r')
xlabel('Zsc values of PAC')
ylabel('Repetation of Observations')
title (' Permutation test of PAC; Fi=5; Ch=106')
set(gca,'fontsize',14)
HighQualityFigs('PAC_PA_Permutation')

% stem plots of p-values for all channels per finger with hline in 0.05
Fi=1;
% for flexion extension
PValuesDirect=[];
PValuesAvg=[];
for ch=1:256
    PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)];
    PValuesAvg=[PValuesAvg; FingerPAC(Fi).PValues.AvgHG.ValueCh(ch)];
end
figure(1);
subplot(2,1,1)
stem(PValuesDirect)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('HG-Direct-LFO')
subplot(2,1,2)
stem(PValuesAvg)
xlim([0,256])
hline(0.05,'r')
set(gca,'fontsize',16)
title('HG-Avg-LFO')

HighQualityFigs('WholeMovement_PAC_PA')

%% 2-3: sig channels and bad channels and brain area channels 
% for direct HG
SigChs_Fingers=logical(zeros(256,5));
for Fi=1:5
    PValuesDirect=[];
    for ch=1:256
        PValuesDirect=[PValuesDirect; FingerPAC(Fi).PValues.WholeHG.ValueCh(ch)];
    end
    SigChs_Fingers(find(PValuesDirect<=0.05),Fi)=1; %significant    
end

% giting rid of bad channels
SigChs_Fingers_Cleaned=logical(zeros(256,5));
for Fi=1:5
    SigChNum1=sum(SigChs_Fingers(:,Fi)); %significant  
    SigChNum2=sum(Selected_Chs.*SigChs_Fingers(:,Fi));
    SigChs_Fingers_Cleaned(:,Fi)=Selected_Chs.*SigChs_Fingers(:,Fi);
end

% common significant channels across fingers
SigChs_Common=logical(zeros(256,1));
SigChs_Common=SigChs_Fingers_Cleaned(:,1).*SigChs_Fingers_Cleaned(:,2).*SigChs_Fingers_Cleaned(:,3).*...
    SigChs_Fingers_Cleaned(:,4).*SigChs_Fingers_Cleaned(:,5);

% all significant channels across fingers
SigChs_All=logical(zeros(256,1));
SigChs_All=logical(sum(SigChs_Fingers_Cleaned,2));

% showing significant channels on brain map per finger and all finger
% toghether
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\EC171_rh_pial.mat')
load('E:\ECoGLeapMotion\DataPatientTwo\ImagingData\TDT_elecs_all.mat')

% for direct HG
for Fi=1:6
    if Fi<6
        subplot(2,3,Fi)
        ctmr_gauss_plot(cortex,elecmatrix(65:320,:),zeros(256,1),'rh'); % rho is a 256ch vector
        el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
        el_add(elecmatrix(find(SigChs_Fingers_Cleaned(:,Fi)==1)+65,:),'msize',5,'color','r');
        title(['Finger',num2str(Fi)]);
    elseif Fi==6
        subplot(2,3,Fi)
        ctmr_gauss_plot(cortex,elecmatrix(65:320,:),zeros(256,1),'rh'); % rho is a 256ch vector
        el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
        el_add(elecmatrix(find(SigChs_All==1)+65,:),'msize',5,'color','r');
        title(['All Fingers']);
        
    end
    
end

HighQualityFigs('All_SigChs_Cleared')


% number of significant ch per brain areas
SigChs_PM = SigChs_All.*Selected_PMChs;
sum(Selected_PMChs)
sum(SigChs_PM)
SigChs_SM = SigChs_All.*Selected_SMChs;
sum(Selected_SMChs)
sum(SigChs_SM)
SigChs_Hand = SigChs_All.*Selected_HandChs;
sum(Selected_HandChs)
sum(SigChs_Hand)

%save('E:\ECoGLeapMotion\DataPatientTwo\github_Branch_V2/SigChs_Logical.mat','SigChs_Fingers_Cleaned','SigChs_All','SigChs_PM','SigChs_SM','SigChs_Hand')
