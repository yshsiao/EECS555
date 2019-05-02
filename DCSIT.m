% for 2-user OFDM with delay channel information
clear;
batch_size = 64;
times = 10000;
% let N0 = 1
EbdB = 0:0.5:20;
Eb = 10.^(EbdB/10);
N0=1;
e_rate1 = zeros(1,length(Eb));
e_rate2 = zeros(1,length(Eb));
recode_time = ones(1,length(Eb))*times;
for qq = 1:length(Eb)
    for kk = 1:times
        infobit_i1 = rand(batch_size,2)>0.5;% first column for symbol 1, second for symbol 2
        infobit_q1 = rand(batch_size,2)>0.5;% first user denote 1
        infobit_i2 = rand(batch_size,2)>0.5;% first column for symbol 1, second for symbol 2
        infobit_q2 = rand(batch_size,2)>0.5;% second user denote 2
        % modulation m-file
        infofreq1 = (infobit_i1*2-1 + (infobit_q1*2-1)*1j);
        infofreq2 = (infobit_i2*2-1 + (infobit_q2*2-1)*1j);
        %%%%%%
        %define channel here, first is the channel length, 3 time slot
        chlength = 5;
        constant = 0; % define if the impulse response is poisson and decay or constant(both are multiple by Rayleigh coefficient)
        cha_fadingt11_1 = channel(chlength, constant,5, 0.9,1)';% channel_length, constant or not, arrival rate, threshold of the first path, decay speed
        cha_fadingt11_2 = channel(chlength, constant,5, 0.9,1)';
        cha_fadingt12_1 = channel(chlength, constant,5, 0.9,1)';
        cha_fadingt12_2 = channel(chlength, constant,5, 0.9,1)';
        cha_long1_1 = [cha_fadingt11_1;zeros(batch_size-chlength,1)];
        cha_long1_2 = [cha_fadingt11_2;zeros(batch_size-chlength,1)];
        cha_long2_1 = [cha_fadingt12_1;zeros(batch_size-chlength,1)];
        cha_long2_2 = [cha_fadingt12_2;zeros(batch_size-chlength,1)];
        Ht11_1 = fft(cha_long1_1);
        Ht11_2 = fft(cha_long1_2);
        Ht12_1 = fft(cha_long2_1);
        Ht12_2 = fft(cha_long2_2);
        cha_fadingt21_1 = channel(chlength, constant,5, 0.9,1)';% channel_length, constant or not, arrival rate, threshold of the first path, decay speed
        cha_fadingt21_2 = channel(chlength, constant,5, 0.9,1)';
        cha_fadingt22_1 = channel(chlength, constant,5, 0.9,1)';
        cha_fadingt22_2 = channel(chlength, constant,5, 0.9,1)';
        cha_long1_1 = [cha_fadingt21_1;zeros(batch_size-chlength,1)];
        cha_long1_2 = [cha_fadingt21_2;zeros(batch_size-chlength,1)];
        cha_long2_1 = [cha_fadingt22_1;zeros(batch_size-chlength,1)];
        cha_long2_2 = [cha_fadingt22_2;zeros(batch_size-chlength,1)];
        Ht21_1 = fft(cha_long1_1);
        Ht21_2 = fft(cha_long1_2);
        Ht22_1 = fft(cha_long2_1);
        Ht22_2 = fft(cha_long2_2);
        cha_fadingt31_1 = channel(chlength, constant,5, 0.9,1)';% channel_length, constant or not, arrival rate, threshold of the first path, decay speed
        cha_fadingt31_2 = channel(chlength, constant,5, 0.9,1)';
        cha_fadingt32_1 = channel(chlength, constant,5, 0.9,1)';
        cha_fadingt32_2 = channel(chlength, constant,5, 0.9,1)';
        cha_long1_1 = [cha_fadingt31_1;zeros(batch_size-chlength,1)];
        cha_long1_2 = [cha_fadingt31_2;zeros(batch_size-chlength,1)];
        cha_long2_1 = [cha_fadingt32_1;zeros(batch_size-chlength,1)];
        cha_long2_2 = [cha_fadingt32_2;zeros(batch_size-chlength,1)];
        Ht31_1 = fft(cha_long1_1);
        Ht31_2 = fft(cha_long1_2);
        Ht32_1 = fft(cha_long2_1);
        Ht32_2 = fft(cha_long2_2);
        %%%%% We can make Tx know the channel or not
        % ifft for OFDM
        %trans_sig = ifft(infofreq)*sqrt(Eb(20))*sqrt(batch_size); % this
        %is no precoding for transmitted signals
        % precoding here if needed
        trans_sig1 = ifft(infofreq1)*sqrt(Eb(qq))*sqrt(batch_size);
        trans_sig2 = ifft(infofreq2)*sqrt(Eb(qq))*sqrt(batch_size);
        trans_f3_1 = Ht12_1.*infofreq1(:,1) + Ht12_2.*infofreq1(:,2);
        trans_f3_2 = Ht21_1.*infofreq2(:,1) + Ht21_2.*infofreq2(:,2);
        trans_sig3 = ifft(trans_f3_1 + trans_f3_2)*sqrt(Eb(qq))*sqrt(batch_size);
        CP = 16;
        trans_sig_cp1 = [trans_sig1(end-CP+1:end,:);trans_sig1];
        trans_sig_cp2 = [trans_sig2(end-CP+1:end,:);trans_sig2];
        trans_sig_cp3 = [trans_sig3(end-CP+1:end,:);trans_sig3];
        % define the received signals here
        pass_channel_sigt11_1 =  conv(trans_sig_cp1(:,1),cha_fadingt11_1);
        pass_channel_sigt11_2 =  conv(trans_sig_cp1(:,2),cha_fadingt11_2);
        pass_channel_sigt12_1 =  conv(trans_sig_cp1(:,1),cha_fadingt12_1);
        pass_channel_sigt12_2 =  conv(trans_sig_cp1(:,2),cha_fadingt12_2);
        received_sigt11 = pass_channel_sigt11_1 + pass_channel_sigt11_2 + (randn(size(pass_channel_sigt11_1)) + randn(size(pass_channel_sigt11_1))*1j)*(sqrt(N0/2));
        received_sigt12 = pass_channel_sigt12_1 + pass_channel_sigt12_2 + (randn(size(pass_channel_sigt12_2)) + randn(size(pass_channel_sigt12_2))*1j)*(sqrt(N0/2));
        pass_channel_sigt21_1 =  conv(trans_sig_cp2(:,1),cha_fadingt21_1);
        pass_channel_sigt21_2 =  conv(trans_sig_cp2(:,2),cha_fadingt21_2);
        pass_channel_sigt22_1 =  conv(trans_sig_cp2(:,1),cha_fadingt22_1);
        pass_channel_sigt22_2 =  conv(trans_sig_cp2(:,2),cha_fadingt22_2);
        received_sigt21 = pass_channel_sigt21_1 + pass_channel_sigt21_2 + (randn(size(pass_channel_sigt21_1)) + randn(size(pass_channel_sigt21_1))*1j)*(sqrt(N0/2));
        received_sigt22 = pass_channel_sigt22_1 + pass_channel_sigt22_2 + (randn(size(pass_channel_sigt22_2)) + randn(size(pass_channel_sigt22_2))*1j)*(sqrt(N0/2));
        pass_channel_sigt31_1 =  conv(trans_sig_cp3(:,1),cha_fadingt31_1);
        %pass_channel_sigt31_2 =  conv(trans_sig_cp3(:,2),cha_fadingt31_2);
        pass_channel_sigt32_1 =  conv(trans_sig_cp3(:,1),cha_fadingt32_1);
        %pass_channel_sigt32_2 =  conv(trans_sig_cp3(:,2),cha_fadingt32_2);
        received_sigt31 = pass_channel_sigt31_1 + (randn(size(pass_channel_sigt31_1)) + randn(size(pass_channel_sigt31_1))*1j)*(sqrt(N0/2));
        received_sigt32 = pass_channel_sigt32_1 + (randn(size(pass_channel_sigt32_1)) + randn(size(pass_channel_sigt32_1))*1j)*(sqrt(N0/2));
        % remove CP and tail of the signals and do fft
        rmcp_ch_sigt11 = received_sigt11(CP+1:end-chlength+1,1);
        received_symbolt11 = fft(rmcp_ch_sigt11);
        rmcp_ch_sigt12 = received_sigt12(CP+1:end-chlength+1,1);
        received_symbolt12 = fft(rmcp_ch_sigt12);
        rmcp_ch_sigt21 = received_sigt21(CP+1:end-chlength+1,1);
        received_symbolt21 = fft(rmcp_ch_sigt21);
        rmcp_ch_sigt22 = received_sigt22(CP+1:end-chlength+1,1);
        received_symbolt22 = fft(rmcp_ch_sigt22);
        rmcp_ch_sigt31 = received_sigt31(CP+1:end-chlength+1,1);
        received_symbolt31 = fft(rmcp_ch_sigt31);
        rmcp_ch_sigt32 = received_sigt32(CP+1:end-chlength+1,1);
        received_symbolt32 = fft(rmcp_ch_sigt32);
        % uncoded hard decision
        %decode1_i = real(received_symbol1./H1_1)>0; % need to change if do
        %some precoding
        %decode1_q = imag(received_symbol1./H1_1)>0;
        L2_u1 = received_symbolt31./Ht31_1 - received_symbolt21;
        L1_u2 = received_symbolt32./Ht32_1 - received_symbolt12;
        L1_u1 = received_symbolt11;
        L2_u2 = received_symbolt22;
        decode1_1 = (L1_u1.*Ht12_2 - L2_u1.*Ht11_2)./(Ht11_1.*Ht12_2 - Ht11_2.*Ht12_1);
        decode1_2 = -1*(L1_u1.*Ht12_1 - L2_u1.*Ht11_1)./(Ht11_1.*Ht12_2 - Ht11_2.*Ht12_1);
        decode2_1 = (L1_u2.*Ht22_2 - L2_u2.*Ht21_2)./(Ht21_1.*Ht22_2 - Ht21_2.*Ht22_1);
        decode2_2 = -1*(L1_u2.*Ht22_1 - L2_u2.*Ht21_1)./(Ht21_1.*Ht22_2 - Ht21_2.*Ht22_1);
        decode1_i = real([decode1_1 decode1_2])>0;
        decode1_q = imag([decode1_1 decode1_2])>0;
        decode2_i = real([decode2_1 decode2_2])>0;
        decode2_q = imag([decode2_1 decode2_2])>0;
        error1 = [(decode1_i~=infobit_i1);(decode1_q~=infobit_q1)];
        error2 = [(decode2_i~=infobit_i2);(decode2_q~=infobit_q2)];
        e_rate1(qq) = e_rate1(qq) + sum(sum(error1));
        e_rate2(qq) = e_rate2(qq) + sum(sum(error2));
        if e_rate1(qq) + e_rate2(qq) > 50000
            recode_time(qq) = kk;
            break;
        end
    end
end
e_rate1 = e_rate1./recode_time/batch_size/4;
e_rate2 = e_rate2./recode_time/batch_size/4;
RayBPSK_e_rate = 1/2-1/2*sqrt(Eb./(1+Eb));
AWGNBPSK_e_rate = qfunc(sqrt(2*Eb));
figure(1);
semilogy(EbdB,e_rate1,'LineWidth',2);
hold on;
semilogy(EbdB,e_rate2,'LineWidth',2);
semilogy(EbdB,RayBPSK_e_rate,'LineWidth',2);
semilogy(EbdB,AWGNBPSK_e_rate,'LineWidth',2);
axis([-inf inf 1e-4 1]);
hline(0.01);
legend('OFDM QPSK error rate(user1)','OFDM QPSK error rate(user2)','Theoretical Rayleigh fading QPSK error rate','AWGN BER');
title('Simulation and theoretical error rate');
hold off;
