load('ECG_database');
%-------------------------------preparation-------------------------------
%convert from raw units to the physical units,Gain=200,base=0,length=5000
%clean ECG 
Data1 = Data1/200;
figure
plot(Data1);ylim([-0.5 1]);
title('Clean ECG Signal');
xlabel('Samples Index');
ylabel('Amplititude (mV)') ;grid
% White Gaussian Noise(WN)
noise_wn = wn/10;
wn_data = noise_wn+Data1;
% Baseline Wander Noise(BWN)
noise_bwn = bwn/200;
bwn_data=BWN_data/200;
% Electrode Movement(EMN)
noise_emn = emn/200;
emn_data=EMN_data/200;
% Muscle Artifacts(MAN)
noise_man = man/200;
man_data=MAN_data/200;
% 50 Hz Power Line Interference(PLI)
noise_pli = 0.1*sin(2*pi*50*(1:5000)/500);
pli_data=noise_pli+Data1;

figure;
subplot(511);plot(wn_data); 
title('ECG corrupted by White Gaussian Noise');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid
subplot(512);plot(bwn_data); 
title('ECG corrupted by Baseline Wander Noise');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid
subplot(513);
plot(emn_data); 
title('ECG corrupted by Electrode Movement Noise');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid
subplot(514);
plot(man_data); 
title('ECG corrupted by Muscle Artifacts');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid
subplot(515);
plot(PLI_data/200); 
title('ECG corrupted by Power Line Interference');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid
%----------------------------------filter----------------------------------
%------White Gaussian Noise------
%[en,yn,wn] = LMSfilter(dn,xn,mu,p)
%mu=0.02,p=2
[en_WN_LMS,yn_WN_LMS,wn_WN_LMS] = LMSfilter(wn_data,noise_wn, 0.02, 2);
%[en,yn,wn]=NLMSfilter(dn,xn,mu,p,a)
%mu=0.02,p=2,a=0.1
[en_WN_NLMS,yn_WN_NLMS,wn_WN_NLMS] = NLMSfilter(wn_data,noise_wn,0.02,2,0.1);
%[en,yn,wn] = RLSfilter(dn,xn,p,lamda)
%p=2,lamda=1
[en_WN_RLS,yn_WN_RLS,wn_WN_RLS] = RLSfilter(wn_data,noise_wn,2,1);

figure;
subplot(411); plot(wn_data); 
title('ECG corrupted by White Gaussian Noise');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(412); plot(en_WN_LMS); ylim([-1 1]);
title('LMS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(413); plot(en_WN_NLMS); ylim([-1 1]);
title('NLMS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(414); plot(en_WN_RLS); ylim([-1 1]);
title('RLS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

%------Baseline Wander Noise------
%[en,yn,wn] = LMSfilter(dn,xn,mu,p)
%mu=0.02,p=2
[en_BWN_LMS,yn_BWN_LMS,wn_BWN_LMS] = LMSfilter(bwn_data,noise_bwn, 0.02, 2);
%[en,yn,wn]=NLMSfilter(dn,xn,mu,p,a)
%mu=0.02,p=2,a=0.1
[en_BWN_NLMS,yn_BWN_NLMS,wn_BWN_NLMS] = NLMSfilter(bwn_data,noise_bwn,0.02,2,0.1);
%[en,yn,wn] = RLSfilter(dn,xn,p,lamda)
%p=2,lamda=1
[en_BWN_RLS,yn_BWN_RLS,wn_BWN_RLS] = RLSfilter(bwn_data,noise_bwn,2,1);
figure;
subplot(411); plot(bwn_data); 
title('ECG corrupted by Baseline Wander Noise');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(412); plot(en_BWN_LMS); ylim([-1 1]);
title('LMS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(413); plot(en_BWN_NLMS); ylim([-1 1]);
title('NLMS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(414); plot(en_BWN_RLS); ylim([-1 1]);
title('RLS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

%------Electrode Movement Noise------
%[en,yn,wn] = LMSfilter(dn,xn,mu,p)
%mu=0.02,p=2
[en_EMN_LMS,yn_EMN_LMS,wn_EMN_LMS] = LMSfilter(emn_data,noise_emn, 0.02, 2);
%[en,yn,wn]=NLMSfilter(dn,xn,mu,p,a)
%mu=0.02,p=2,a=0.1
[en_EMN_NLMS,yn_EMN_NLMS,wn_EMN_NLMS] = NLMSfilter(emn_data,noise_emn,0.02,2,0.1);
%[en,yn,wn] = RLSfilter(dn,xn,p,lamda)
%p=2,lamda=1
[en_EMN_RLS,yn_EMN_RLS,wn_EMN_RLS] = RLSfilter(emn_data,noise_emn,2,1);

figure;
subplot(411); plot(emn_data); 
title('ECG corrupted by Electrode Movement Noise');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(412); plot(en_EMN_LMS); ylim([-1 1]);
title('LMS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(413); plot(en_EMN_NLMS); ylim([-1 1]);
title('NLMS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(414); plot(en_EMN_RLS); ylim([-1 1]);
title('RLS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid



%------Muscle Artifacts Noise------
%[en,yn,wn] = LMSfilter(dn,xn,mu,p)
%mu=0.02,p=2
[en_MAN_LMS,yn_MAN_LMS,wn_MAN_LMS] = LMSfilter(man_data,noise_man, 0.02, 2);
%[en,yn,wn]=NLMSfilter(dn,xn,mu,p,a)
%mu=0.02,p=2,a=0.1
[en_MAN_NLMS,yn_MAN_NLMS,wn_MAN_NLMS] = NLMSfilter(man_data,noise_man,0.02,2,0.1);
%[en,yn,wn] = RLSfilter(dn,xn,p,lamda)
%p=2,lamda=1
[en_MAN_RLS,yn_MAN_RLS,wn_MAN_RLS] = RLSfilter(man_data,noise_man,2,1);


figure;
subplot(411); plot(man_data); 
title('ECG corrupted by Muscle Artifacts Noise');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(412); plot(en_MAN_LMS); ylim([-1 1]);
title('LMS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(413); plot(en_MAN_NLMS); ylim([-1 1]);
title('NLMS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(414); plot(en_MAN_RLS); ylim([-1 1]);
title('RLS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid


%------Power Line Interference------
%[en,yn,wn] = LMSfilter(dn,xn,mu,p)
%mu=0.02,p=2
[en_PLI_LMS,yn_PLI_LMS,wn_PLI_LMS] = LMSfilter(pli_data,noise_pli, 0.02, 2);
%[en,yn,wn]=NLMSfilter(dn,xn,mu,p,a)
%mu=0.02,p=2,a=0.1
[en_PLI_NLMS,yn_PLI_NLMS,wn_PLI_NLMS] = NLMSfilter(pli_data,noise_pli,0.02,2,0.1);
%[en,yn,wn] = RLSfilter(dn,xn,p,lamda)
%p=2,lamda=1
[en_PLI_RLS,yn_PLI_RLS,wn_PLI_RLS] = RLSfilter(pli_data,noise_pli,2,1);

figure;
subplot(411); plot(pli_data); 
title('ECG corrupted by Power Line Interference');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(412); plot(en_PLI_LMS); ylim([-1 1]);
title('LMS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(413); plot(en_PLI_NLMS); ylim([-1 1]);
title('NLMS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid

subplot(414); plot(en_PLI_RLS); ylim([-1 1]);
title('RLS Filter Response');
xlabel('Samples (n)');
ylabel('Amplititude (mV)') ;grid



%-----------------------------SNR Calculation------------------------------
%White Gaussian Noise
imp_snr_wn_lms = abs(snr(Data1, noise_wn) - snr(Data1, Data1-en_WN_LMS'));
imp_snr_wn_nlms = abs(snr(Data1, noise_wn) - snr(Data1, Data1-en_WN_NLMS'));
imp_snr_wn_rls = abs(snr(Data1, noise_wn) - snr(Data1, Data1-en_WN_RLS')) ;
%Baseline Wander Noise
imp_snr_bwn_lms = abs(snr(Data1, noise_bwn) - snr(Data1, Data1-en_BWN_LMS'));
imp_snr_bwn_nlms = abs(snr(Data1, noise_bwn) - snr(Data1, Data1-en_BWN_NLMS'));
imp_snr_bwn_rls = abs(snr(Data1, noise_bwn) - snr(Data1, Data1-en_BWN_RLS')) ;
%Electrode Movement Noise
imp_snr_emn_lms = abs(snr(Data1, noise_emn) - snr(Data1, Data1-en_EMN_LMS'));
imp_snr_emn_nlms = abs(snr(Data1, noise_emn) - snr(Data1, Data1-en_EMN_NLMS'));
imp_snr_emn_rls = abs(snr(Data1, noise_emn) - snr(Data1, Data1-en_EMN_RLS')) ;
%Muscle Artifacts Noise
imp_snr_man_lms = abs(snr(Data1, noise_man) - snr(Data1, Data1-en_MAN_LMS'));
imp_snr_man_nlms = abs(snr(Data1, noise_man) - snr(Data1, Data1-en_MAN_NLMS'));
imp_snr_man_rls = abs(snr(Data1, noise_man) - snr(Data1, Data1-en_MAN_RLS')) ;
% 50 Hz Power Line Interference
imp_snr_pli_lms = abs(snr(Data1, noise_pli) - snr(Data1, Data1-en_PLI_LMS'));
imp_snr_pli_nlms = abs(snr(Data1, noise_pli) - snr(Data1, Data1-en_PLI_NLMS'));
imp_snr_pli_rls = abs(snr(Data1, noise_pli) - snr(Data1, Data1-en_PLI_RLS')) ;

%Mean Squared Error
%--------------------------MSE of reconstructed ECG-------------------------
%White Gaussian Noise
mse_wn_lms = (1/5000)*sum((Data1-en_WN_LMS').^2);
mse_wn_nlms = (1/5000)*sum((Data1-en_WN_NLMS').^2);
mse_wn_rls = (1/5000)*sum((Data1-en_WN_RLS').^2); 
%Baseline Wander Noise
mse_bwn_lms = (1/5000)*sum((Data1-en_BWN_LMS').^2);
mse_bwn_nlms = (1/5000)*sum((Data1-en_BWN_NLMS').^2);
mse_bwn_rls = (1/5000)*sum((Data1-en_BWN_RLS').^2);
%Electrode Movement Noise
mse_emn_lms = (1/5000)*sum((Data1-en_EMN_LMS').^2);
mse_emn_nlms = (1/5000)*sum((Data1-en_EMN_NLMS').^2);
mse_emn_rls = (1/5000)*sum((Data1-en_EMN_RLS').^2);
%Muscle Artifacts Noise 
mse_man_lms = (1/5000)*sum((Data1-en_MAN_LMS').^2);
mse_man_nlms = (1/5000)*sum((Data1-en_MAN_LMS').^2);
mse_man_rls = (1/5000)*sum((Data1-en_MAN_RLS').^2);
% 50 Hz Power Line Interference
mse_pli_lms = (1/5000)*sum((Data1-en_PLI_LMS').^2);
mse_pli_nlms = (1/5000)*sum((Data1-en_PLI_NLMS').^2);
mse_pli_rls = (1/5000)*sum((Data1-en_PLI_RLS').^2);

%percentage root-mean-squared difference
%--------------------------PRD of reconstructed ECG-------------------------
%White Gaussian Noise
prd_wn_lms = sqrt((1/5000)*sum((Data1-en_WN_LMS').^2))*100;
prd_wn_nlms = sqrt((1/5000)*sum((Data1-en_WN_NLMS').^2))*100;
prd_wn_rls = sqrt((1/5000)*sum((Data1-en_WN_RLS').^2))*100;
%Baseline Wander Noise
prd_bwn_lms = sqrt((1/5000)*sum((Data1-en_BWN_LMS').^2))*100;
prd_bwn_nlms = sqrt((1/5000)*sum((Data1-en_BWN_NLMS').^2))*100;
prd_bwn_rls = sqrt((1/5000)*sum((Data1-en_BWN_RLS').^2))*100;
%Electrode Movement Noise
prd_emn_lms = sqrt((1/5000)*sum((Data1-en_EMN_LMS').^2))*100;
prd_emn_nlms = sqrt((1/5000)*sum((Data1-en_EMN_NLMS').^2))*100;
prd_emn_rls = sqrt((1/5000)*sum((Data1-en_EMN_RLS').^2))*100;
%Muscle Artifacts Noise 
prd_man_lms = sqrt((1/5000)*sum((Data1-en_MAN_LMS').^2))*100;
prd_man_nlms = sqrt((1/5000)*sum((Data1-en_MAN_NLMS').^2))*100;
prd_man_rls = sqrt((1/5000)*sum((Data1-en_MAN_RLS').^2))*100;
% 50 Hz Power Line Interference
prd_pli_lms = sqrt((1/5000)*sum((Data1-en_PLI_LMS').^2))*100;
prd_pli_nlms = sqrt((1/5000)*sum((Data1-en_PLI_NLMS').^2))*100;
prd_pli_rls = sqrt((1/5000)*sum((Data1-en_PLI_RLS').^2))*100;


%-----------------------Power Spectrum Density Analysis--------------------
FS=500;
nfft = 2^nextpow2(length(pli_data));
Pxx = abs(fft(pli_data,nfft)).^2/length(pli_data)/FS;
PSD_PLI = dspdata.psd(Pxx(1:length(Pxx)/2),'Fs',FS);  


nfft_lms = 2^nextpow2(length(en_PLI_LMS));
Pxx_lms = abs(fft(en_PLI_LMS,nfft_lms)).^2/length(en_PLI_LMS)/FS;
PSD_PLI_lms = dspdata.psd(Pxx_lms(1:length(Pxx_lms)/2),'Fs',FS);  


nfft_nlms = 2^nextpow2(length(en_PLI_NLMS));
Pxx_nlms = abs(fft(en_PLI_NLMS,nfft_nlms)).^2/length(en_PLI_NLMS)/FS;
PSD_PLI_nlms = dspdata.psd(Pxx_nlms(1:length(Pxx_nlms)/2),'Fs',FS);  

nfft_rls = 2^nextpow2(length(en_PLI_RLS));
Pxx_rls = abs(fft(en_PLI_RLS,nfft_rls)).^2/length(en_PLI_RLS)/FS;
PSD_PLI_rls = dspdata.psd(Pxx_rls(1:length(Pxx_rls)/2),'Fs',FS);  

figure
subplot(2,2,1);
plot(PSD_PLI); 
title('PSD plot of ECG corrupted by 50Hz PLI');
xlabel('Frequency(HZ)');
ylabel('Power/Frequency (dB/Hz)');
subplot(2,2,2);
plot(PSD_PLI_lms); 
title('PSD plot of ECG filtered by LMS algorithm');
xlabel('Frequency(HZ)');
ylabel('Power/Frequency (dB/Hz)');
subplot(2,2,3);
plot(PSD_PLI_nlms); 
title('PSD plot of ECG filtered by NLMS algorithm');
xlabel('Frequency(HZ)');
ylabel('Power/Frequency (dB/Hz)');
subplot(2,2,4);
plot(PSD_PLI_rls); 
title('PSD plot of ECG filtered by RLS algorithm');
xlabel('Frequency(HZ)');
ylabel('Power/Frequency (dB/Hz)');

%--------------------------error Magnitude---------------------------------
figure;
plot(yn_WN_LMS'-noise_wn,'r'); hold on;
plot(yn_WN_NLMS'-noise_wn,'g'); hold on;
plot(yn_WN_RLS'-noise_wn,'b'); hold on;
title('Plot of Error Magnitude vs. Sample Index of three algorithms for White Noise');
legend('LMS','NLMS','RLS');
xlabel('Samples (n)');
ylabel('Error Magnitude');

figure;
plot(yn_PLI_LMS'-noise_pli,'r'); hold on;
plot(yn_PLI_NLMS'-noise_pli,'g'); hold on;
plot(yn_PLI_RLS'-noise_pli,'b'); hold on;
title('Plot of Error Magnitude vs. Sample Index of three algorithms for 50HZ PLI');
legend('LMS','NLMS','RLS');
xlabel('Samples (n)');
ylabel('Error Magnitude');


