clear;
close all;
%Initialization
L=3;
lambda1=2;
lambda2=0.01;
SNR=5;
iteration_times=500;
parameter.num_of_code=4096;
parameter.frequency_offset=0.3;
parameter.phase_offset=0;

signal_type=[2,2,2,3,3,3,6,6,6;8,16,32,8,16,64,4,8,16];%1.ASK;2.PSK;3.QAM;4.FSK;5.MSK;6.PAM.
% signal_type=[3,3,3,;4,16,64];%1.ASK;2.PSK;3.QAM;4.FSK;5.MSK;6.PAM.
M=size(signal_type,2);
%     for SNR=-18:2:0
        parameter.SNR=SNR;
% for num_of_code=1024:1024:4096*2.5
    disp('==================================================');
    parameter.signal_type=repmat(signal_type,1,parameter.num_of_code);
number_of_dictionary=size(parameter.signal_type,2);
%     parameter.num_of_code=num_of_code;
[signal]=generating_signals(parameter);
% load signal_snr_m10;
signal=[real(signal);imag(signal)];
%High cumulant analysis
% high_cumulants=High_Cumulants(signal);
%=======================================================================
% for lambda=0:0.1:1
%     for SNR=[-20:2:-12,-7,0,2]
param.dictionary=signal(:,1:M);
% param.step_1=50000;
% param.step_2=40000;
param.lambda1=lambda1;
param.lambda2=lambda2;
param.M=M;
param.L=L;
param.numIteration=iteration_times;
param.coeff=OMP(param.dictionary,signal,L);
disp(['lambda=',num2str(lambda1)]);
sdlc_ori_start=clock;
output_sdlc_ori=SDLC_proxi_ori(signal,param,signal_type);%*****************************
sdlc_ori_end=clock;
coeff_sdlc_ori=full(output_sdlc_ori.coeff);
[~,label_sdlc_ori]=max(abs(coeff_sdlc_ori));
label=repmat(1:size(signal_type,2),1,size(label_sdlc_ori,2)/size(signal_type,2));
error_ratio_sdlc_ori=error_ratio(label_sdlc_ori,label);
time_sdlc_ori=etime(sdlc_ori_end,sdlc_ori_start);

%---------------------------------------------------------------------------------
sdlc_my_start=clock;
output_sdlc_my=SDLC_proxi(signal,param,signal_type);%*********************************
% output_sdlc=SDLC_advanced(signal,param,signal_type);
sdlc_my_end=clock;
coeff_sdlc_my=full(output_sdlc_my.coeff);
[~,label_sdlc_my]=max(abs(coeff_sdlc_my));
error_ratio_sdlc_my=error_ratio(label_sdlc_my,label);
time_sdlc_my=etime(sdlc_my_end,sdlc_my_start);

%---------------------------------------------------------------------------------
sdlc_advanced_start=clock;
output_sdlc_advanced=SDLC_proxi_advanced(signal,param,signal_type);
sdlc_advanced_end=clock;
coeff_sdlc_advanced=full(output_sdlc_advanced.coeff);
[~,label_sdlc_advanced]=max(abs(coeff_sdlc_advanced));
error_ratio_sdlc_advanced=error_ratio(label_sdlc_advanced,label);
time_sdlc_advanced=etime(sdlc_advanced_end,sdlc_advanced_start);

%----------------------------------------------------------------------------------
% error_ratio_sdlc=sum(sum(label_sdlc~=label))/size(label_sdlc,2);
% error_ratio_ksvd=sum(sum(label_ksvd~=label))/size(label_ksvd,2);
disp(['SNR=',num2str(SNR)]);
filename_record=['counter_october_25_','record_lambda_',num2str(lambda1),'_size_',num2str(numel(signal)),'_snr',num2str(SNR),'.mat'];
save(filename_record);
disp(['Time of original SDLC=',num2str(time_sdlc_ori)]);
disp(['Time of my SDLC=',num2str(time_sdlc_my)]);
disp(['Time of advanced SDLC=',num2str(time_sdlc_advanced)]);
disp(['Error Ratio of origianl SDLC=',num2str(error_ratio_sdlc_ori)]);

disp(['Error Ratio of my SDLC=',num2str(error_ratio_sdlc_my)]);

disp(['Error Ratio of my advanced SDLC=',num2str(error_ratio_sdlc_advanced)]);
%     end
% end
load chirp;
sound(y,Fs)
%=====================================================================================================
%-----------------------------------------------------------------------------------------------------
% paramksvd.L = L;
% paramksvd.K = M;% number of dictionary elements
% paramksvd.numIteration = iteration_times; % number of iteration to execute the K-SVD algorithm.
% paramksvd.errorFlag = 0; % decompose signals until a certain error is reached. do not use fix number of coefficients.
% paramksvd.preserveDCAtom = 0;
% paramksvd.InitializationMethod =  'GivenMatrix';%Initialization
% paramksvd.initialDictionary=param.dictionary;
% % Training Dictionary Set
% ksvd_start=clock;
% error_ratio_ksvd=error_ratio(label_ksvd,label);
% [dictionary,output_ksvd]=KSVD(signal,paramksvd);
% ksvd_end=clock;
% time_ksvd=etime(ksvd_end,ksvd_start);
% %----------------------------------------------
% coeff_ksvd=full(output_ksvd.CoefMatrix);
% [~,label_ksvd]=max(abs(coeff_ksvd));
% disp(['Error Ratio of K-SVD=',num2str(error_ratio_ksvd)]);
% disp(['Time of K-SVD=',num2str(time_ksvd)]);