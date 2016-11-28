function [output]= generating_signals(param)
%SAMPLES_TRANING Generate the test signal
%Initinalization
number_of_code=param.num_of_code;%Number of codes
frequency_offset=param.frequency_offset;
phase_offset=param.phase_offset;
SNR=param.SNR;
signal_type=param.signal_type;
output=Inf*ones(number_of_code,size(signal_type,2));
%The ith element in the first column denotes the modulation order, the jth element in the second column denotes the modulation type.
%Generate the signals
%--------------------------------------------------------------------------
for i=1:size(signal_type,2)
    M=signal_type(2,i);
code=repmat((0:M-1)',number_of_code/M,1);
%Select the proper channel
modulation_type=signal_type(1,i);
switch modulation_type%Select the modulation type
    case 1;%%ASK
        code_ASK=code-(M-1)/2;%
        s_ASK=code_ASK.*exp(-1i*phase_offset);
        signal=s_ASK;
%         title_str=[num2str(M),'ASK '];
%--------------------------------------------------------------------------
    case 2;%%PSK
        PSK_mod=comm.PSKModulator('ModulationOrder',M);
%         PSK_Demod=comm.PSKDemodulator('ModulationOrder',M);
        PSK_mod.PhaseOffset = phase_offset;
        s_PSK=step(PSK_mod,code);
        signal=s_PSK;
%         title_str=[num2str(M),'PSK '];
%--------------------------------------------------------------------------
    case 3;%%QAM
        QAM_mod=comm.RectangularQAMModulator('ModulationOrder',M);
%         QAM_Demod=comm.RectangularQAMDemodulator('ModulationOrder',M);
        QAM_mod.PhaseOffset = phase_offset;
        s_QAM=step(QAM_mod,code);
        signal=s_QAM;
%         title_str=[num2str(M),'QAM '];
%--------------------------------------------------------------------------
    case 4;%%FSK
        FSK_mod=comm.FSKModulator('ModulationOrder',M);
%         FSK_Demod=comm.FSKDemodulator('ModulationOrder',M);
        s_FSK=step(FSK_mod,code);
        signal=s_FSK;
%         title_str=[num2str(M),'FSK '];
%--------------------------------------------------------------------------
    case 5;%%MSK
        MSK_mod=comm.MSKModulator('BitInput', true,'InitialPhaseOffset', phase_offset);
%         MSK_Demod=comm.MSKDemodulator('BitOutput', true, 'InitialPhaseOffset', phase_offset);
        s_MSK=step(MSK_mod,code_MSK);
        signal=s_MSK;
%         title_str=[num2str(M),'MSK '];
%--------------------------------------------------------------------------
    case 6;%%PAM     
        PAM_mod=comm.PAMModulator('ModulationOrder',M);
%         PAM_Demod=comm.PAMDemodulator('ModulationOrder',M);
        s_PAM=step(PAM_mod,code).*exp(-1i*phase_offset);
        signal=s_PAM;
%         title_str=[num2str(M),'PAM '];
end
% signal_send_transmitter=Windows_send(signal,sps);%The signals sent through transmitter, the up sampling rate is sps
% %---------------------------------------------------------------------------
% signal_send_channel=channel_function(signal_send_transmitter,channel_type);%The signals through the channel which defined before.
% %---------------------------------------------------------------------------
% signal_receive=Windows_receive(signal_send_channel,sps);%The signals received by receiver, the down sampling rate is sps.
% signal_receive=signal_receive(sps:end);
% rng('default');rng(1);
signal=awgn(signal,SNR);
signal=frequencyoffset(signal,frequency_offset);
output(:,i)=signal;
% signal_receive=awgn(signal,snr);
% signal_train(:,iteration)=signal_receive;
%--------------------------------------------------------------------------
%Analysis and Illustration
%--------------------------------------------------------------------------
%Illustrat the constellation of transmitted signals
% figure(1);
% scatter(real(signal),imag(signal),'k*');
% grid on;
% axis equal;
% title(['Constellation of transmitted signal ',title_str,num2str(number_of_code),' points']);
%Illustrat the constellation of received signals
end
%END OF PROGRAM