function [channel_sig] = channel(chlength, constant, mu, threshold,decay)  
    if constant == 0
        h = exp(-[0:chlength-1]*decay);    % Power Delay Profile
        gaps = exprnd(mu,1,length(h));
        not_zero_loc   = [1, round(cumsum(gaps))+1];  
        nzero_idx      = find(not_zero_loc <= length(h));
        zero_loc_flag                          = ones(1, length(h));
        zero_loc_flag(not_zero_loc(nzero_idx)) = 0;
        h(zero_loc_flag == 1)                  = 0;
        channel_sig                            = 1/sqrt(2)*(randn(size(h)) + 1i*randn(size(h))).*h;
        channel_sig                            = channel_sig/(sqrt(sum(abs(channel_sig).^2)));% make the variance = 1
        % Force the first path to have moderate power
        while abs(channel_sig(1))<threshold*max(abs(channel_sig))
             channel_sig = 1/sqrt(2)*(randn(size(h)) + 1i*randn(size(h))).*h;
             channel_sig = channel_sig/(sqrt(sum(abs(channel_sig).^2)));% make the variance = 1
        end
    else
        channel_sig = 1/sqrt(2)*(randn(1,chlength) + 1i*randn(1,chlength));
        channel_sig = channel_sig/(sqrt(sum(abs(channel_sig).^2)));
    end
end
