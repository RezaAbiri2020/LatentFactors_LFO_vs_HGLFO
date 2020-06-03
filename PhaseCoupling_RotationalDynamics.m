


%% phase-amplitude coupling
% how the phase of LFO is coupled to the power(abs-hilbert) of HG-LFO
load('HG_AllAvgSubBands.mat')
load('Delta_AllSignals.mat')

for Fi=1:5
    Amp=abs(hilbert(HG_DeltaofAvgHG(Fi).Finger));
    Phase=angle(hilbert(Delta_Fingers(Fi).Filtered));
    
    for ch=1:256
        %polarplot(Phase(1:end,ch),Amp(1:end,ch),'-')
        FingerPAC(Fi).ValueCh(ch).Real=(abs(sum(Amp(1:end,ch).*exp(1i.*Phase(1:end,ch)))))/length(Phase);
        
        AmpCh=Amp(1:end,ch);
        for j=1:10000 % producing shuffled results
            Shuffled_AmpCh = AmpCh(randperm(numel(AmpCh)));
            FingerPAC(Fi).ValueCh(ch).Shuffled(j)=(abs(sum(Shuffled_AmpCh(1:end).*exp(1i.*Phase(1:end,ch)))))/length(Phase);
        end
        % z-score per channel
        FingerPAC(Fi).ValueCh(ch).ZScAll=zscore([FingerPAC(Fi).ValueCh(ch).Real, FingerPAC(Fi).ValueCh(ch).Shuffled]);
        
    end
    
end

% plot of observations; example
figure
hist(FingerPAC(Fi).ValueCh(1).ZScAll,100)
set(gca,'fontsize',16)

%Calculating p-value for all channels
for Fi=1:5
    for ch=1:256
        P_value=length(find(FingerPAC(Fi).ValueCh(ch).ZScAll(2:end)>FingerPAC(Fi).ValueCh(ch).ZScAll(1)))...
            /length(FingerPAC(Fi).ValueCh(ch).ZScAll);
        
        FingerPAC(Fi).ValueCh(ch).PValue=P_value;
    end
end

% stem plots of p-values for all channels per finger with hline in 0.05
Fi=1;
PValuesFi=[];
for ch=1:256
    PValuesFi=[PValuesFi; FingerPAC(Fi).ValueCh(ch).PValue];
end
stem(PValuesFi)
hline(0.05,'r')
set(gca,'fontsize',16)

% showing significant channels on brain map per finger

load('..\ImagingData\EC171_rh_pial.mat')
load('..\ImagingData\TDT_elecs_all.mat')

for Fi=1:5
    PValuesFi=[];
    for ch=1:256
        PValuesFi=[PValuesFi; FingerPAC(Fi).ValueCh(ch).PValue];
    end
    PValuesFi(find(PValuesFi<0.05))=-1; %significant
    PValuesFi(find(PValuesFi>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(-1*PValuesFi),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    title(['Finger',num2str(Fi)]);
    %colorbar
    
    
end

save('FingerPAC.mat','FingerPAC')

%% phase-phase coupling
% how the phase of LFO is coupled to the phase of HG-LFO

for Fi=1:5
    Phase1=angle(hilbert(HG_DeltaofAvgHG(Fi).Finger));
    Phase2=angle(hilbert(Delta_Fingers(Fi).Filtered));
    
    for ch=1:256
        PPC_Value=(sum(exp(1i.*(Phase1(1:end,ch)-Phase2(1:end,ch)))))/length(Phase1);
        FingerPPC(Fi).ValueCh(ch).RealMag=abs(PPC_Value); 
        FingerPPC(Fi).ValueCh(ch).RealAngle= angle(PPC_Value);
        
        Phase2Ch=Phase2(1:end,ch);
        for j=1:10000 % producing shuffled results
            Shuffled_Phase2Ch = Phase2Ch(randperm(numel(Phase2Ch)));
            FingerPPC(Fi).ValueCh(ch).Shuffled(j)=abs((sum(exp(1i.*(Phase1(1:end,ch)-Shuffled_Phase2Ch(1:end)))))/length(Phase1));
        end
        
    end
    
end

% plot of observations; example
figure
hist([FingerPPC(1).ValueCh(1).RealMag, FingerPPC(1).ValueCh(1).Shuffled],100)
set(gca,'fontsize',16)

%Calculating p-value for all channels
for Fi=1:5
    for ch=1:256
        P_value=length(find(FingerPPC(Fi).ValueCh(ch).Shuffled>FingerPPC(Fi).ValueCh(ch).RealMag))...
            /(length(FingerPPC(Fi).ValueCh(ch).Shuffled)+1);
        
        FingerPPC(Fi).ValueCh(ch).PValue=P_value;
    end
end

% stem plots of p-values for all channels per finger with hline in 0.05
Fi=1;
PValuesFi=[];
for ch=1:256
    PValuesFi=[PValuesFi; FingerPPC(Fi).ValueCh(ch).PValue];
end
stem(PValuesFi)
hline(0.05,'r')
set(gca,'fontsize',16)

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
        if FingerPPC(Fi).ValueCh(Ch_num_3(r,c)).PValue<0.05
            polarplot([0 FingerPPC(Fi).ValueCh(Ch_num_3(r,c)).RealAngle],[0 1],'-r','linewidth',1.5)
            rticks('')
            thetaticks('')
        else
            polarplot([0 FingerPPC(Fi).ValueCh(Ch_num_3(r,c)).RealAngle],[0 1],'-b','linewidth',1.5)
            rticks('')
            thetaticks('')
        end
        
        hold on
          
    end
end

% showing significant channels on brain map per finger

load('..\ImagingData\EC171_rh_pial.mat')
load('..\ImagingData\TDT_elecs_all.mat')

for Fi=1:5
    PValuesFi=[];
    for ch=1:256
        PValuesFi=[PValuesFi; FingerPPC(Fi).ValueCh(ch).PValue];
    end
    PValuesFi(find(PValuesFi<0.05))=-1; %significant
    PValuesFi(find(PValuesFi>=0.05))=0; % not significant
    subplot(2,3,Fi)
    ctmr_gauss_plot(cortex,elecmatrix(65:320,:),(-1*PValuesFi),'rh'); % rho is a 256ch vector
    el_add(elecmatrix(65:320,:),'msize',1.7); % for plotting channels on brain
    title(['Finger',num2str(Fi)]);
    %colorbar
    
    
end

save('FingerPPC.mat','FingerPPC')


%% connection of significant channels to Weights of vel/pos prediction
load('HG_AllWeights.mat');



%% rotation dynamics for high gamma in Ch 186, 187
% State space model estimation:
% Xdot=A*X
Samples=amp(:,[49 50]);
% magnify the values
X=1*Samples;
Xdot=diff(X);
X=X(1:(end-1),:);
X=X';
Xdot=Xdot';
A=Xdot*X'*pinv(X*X'); 
% solving for phase plane
save('AMatrix.mat','A')
tspan=[0,10000];
icond={[4, 0]};
figure;
PhasePlane(@system1,tspan,icond,'Xlim',[-5 5],'Ylim',[-5 5],...
    'ArrowHeads',true,'ArrowSize',1);

% rotation dynamics for delta in Ch 186, 187
% State space model estimation:
% Xdot=A*X
Samples=Delta_Fingers(Fi).Hilbert(:,[49 50]);
% magnify the values
X=1*Samples;
Xdot=diff(X);
X=X(1:(end-1),:);
X=X';
Xdot=Xdot';
A=Xdot*X'*pinv(X*X'); 
% solving for phase plane
save('AMatrix.mat','A')
tspan=[0,10000];
icond={[4, 0]};
figure;
PhasePlane(@system1,tspan,icond,'Xlim',[-5 5],'Ylim',[-5 5],...
    'ArrowHeads',true,'ArrowSize',1);


