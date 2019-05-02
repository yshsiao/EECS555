% simulate rayleigh fading and OFDM signals
%parameter definition
batch_size = 64;% sub-carrier number
times = 1000;% simulation times
% let N0 = 1
EbdB = 0:0.5:10;
Eb = 10.^(EbdB/10);
N0=1;
e_rate = zeros(1,length(Eb));
recode_time = ones(1,length(Eb))*times;
for qq = 1:length(Eb)
    for kk = 1:times
        % define the transmit signals
        infobit_i = rand(batch_size,1)>0.5;
        infobit_q = rand(batch_size,1)>0.5;
        infofreq= (infobit_i*2-1 + (infobit_q*2-1)*1j);
        trans_sig = ifft(infofreq)*sqrt(Eb(qq))*sqrt(batch_size);
        % add CP before sending
        %define CP = 600, channel delay spread = 500, both are adjustable
        CP = 16;
        trans_sig_cp = [trans_sig(end-CP+1:end,1);trans_sig];
        %define channel here
        chlength = 5;
        constant = 1; % define if the impulse response is poisson and decay or constant(both are multiple by Rayleigh coefficient)
        channel_fading = channel(chlength, constant,5, 0.9,1)';% channel_length, constant or not, arrival rate, threshold of the first path, decay speed
        pass_channel_sig =  conv(trans_sig_cp,channel_fading);
        received_sig = pass_channel_sig + (randn(size(pass_channel_sig)) + randn(size(pass_channel_sig))*1j)*(sqrt(N0/2));
        % remove CP and tail of the signals and do fft
        rmcp_ch_sig = received_sig(CP+1:end-chlength+1,1);
        decode = fft(rmcp_ch_sig);
        % uncoded hard decision
        % we know the channel fading:
        cha_long = [channel_fading;zeros(batch_size-chlength,1)];
        H = fft(cha_long);
        decode_i = real(decode.*conj(H))>0;
        decode_q = imag(decode.*conj(H))>0;
        error = [(decode_i~=infobit_i);(decode_q~=infobit_q)];
        e_rate(qq) = e_rate(qq) + sum(error);
        if e_rate(qq) > 5000
            recode_time(qq) = kk;
            break;
        end
    end
end
e_rate = e_rate./recode_time/batch_size/2;
RayBPSK_e_rate = 1/2-1/2*sqrt(Eb./(1+Eb));
figure(1);
semilogy(EbdB,e_rate,'b');
hold on;
semilogy(EbdB,RayBPSK_e_rate,'r');
axis([-inf inf 1e-2 1]);
legend('Simulation OFDM QPSK error rate','Theoretical Rayleigh fading QPSK error rate');
title('Simulation and theoretical error rate');
hold off;
