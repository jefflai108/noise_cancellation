function y=SpectralDenoising(filepath, Tms, a, b)

    %% locate the file
%     filepath = '/Users/jefflai108/signals_systems/project2/section 3/pinkNoise_sample.wav';
%     Tms = 20; % ms per frame 
    [x,Fs] = audioread(filepath); 
    frame = Tms*Fs/1000; % samples per frame 
    FFT_coeff = 160;
    
    %% segment the signal
    % The input signal x(t) is analyzed over short segments xw(t) of size Tw
    n_iteration = ceil(length(x)/frame); % total number of frames 
    x_w = zeros(n_iteration,frame);
    
    assert(n_iteration*frame >= length(x) & length(x) >= (n_iteration-1)*frame);
    
    for i=1:n_iteration-1
        x_w(i,:) = x((i-1)*frame+1:i*frame);
    end
    for j=1:length(x(i*frame+1:end))
        x_w(i+1,j) = x(i*frame+j); % pad with zeros at the end 
    end
    
    %% speech/noise classification
    % decide if each frame is speech or noise. The first second of the
    % signal is noise. The rest is speech. 
    n_frames = 1*1000/Tms; % number of noise frames 
    s_frames = n_iteration-n_frames; % number of speech frames 
    n_w = zeros(n_frames,frame);
    s_w = zeros(s_frames,frame);
    
    assert(n_iteration == n_frames+s_frames);
    
    for i=1:n_frames 
        n_w(i,:) = x_w(i,:);
    end 
    
    for i=1:s_frames
        s_w(i,:) = x_w(n_frames+i,:); 
    end 
    
    %% FFT for n_w 
    % Magnitude spectrum N_w is the average fft for each noise frames
    N_w = zeros(1,FFT_coeff);
    for i=1:n_frames
        N_w = N_w + abs(fft(n_w(i,:),FFT_coeff)); % update N_w as more frames come in 
    end 
    N_w = N_w/n_frames; % average of all noise frames 
    
    %% Construct the denoised signal - PART 1
    % First, compute FFT for s_w.  
    % Second, subtract the noise magnitude spectrum from the speech magnitude spectrum 
    % to estimate the magnitude spectrum of the denoised signal.
    % Third, multiply step two with the original phase of speech signal 
    % Inverse Fourier transform to get the time signal yw(t).
    y_1 = [];
    for i=1:s_frames
        freq_dom = max((abs(fft(s_w(i,:),FFT_coeff))-N_w),0).*exp(1j*angle(fft(s_w(i,:),FFT_coeff)));
        y_1 = [y_1, real(ifft(freq_dom',FFT_coeff,'symmetric'))'];
    end 
    
    %% Construct the denoised signal - PART 2
    y = [];
    for i=1:s_frames
        speech_freq = fft(s_w(i,:),FFT_coeff);
        recons_freq = [];
        for i=1:length(speech_freq)
            if abs(speech_freq(i)) > a*N_w(i)
                recons_freq = [recons_freq, abs(speech_freq(i)) - a*N_w(i)];
            else 
                recons_freq = [recons_freq, b*N_w(i)];
            end 
        end 
        freq_dom = recons_freq.*exp(1j*angle(speech_freq));
        y = [y, real(ifft(freq_dom',FFT_coeff,'symmetric'))'];
    end  
    
    %% SNR - PART 1 & 2
    noise_est = 0;
    for i=1:n_frames
        noise_est = noise_est + 10*log10(1/frame*sum(power(n_w(i,:),2)));
    end 
    noise_est = noise_est/n_frames; %noise_power average over the first second
    
    SNR1 = [];
    for i=1:s_frames
        SNR1 = [SNR1, 10*log10(1/frame*sum(power(y_1((i-1)*FFT_coeff+1:i*FFT_coeff),2)))-noise_est]; %SNR for one frame of signal 
    end 
    SNR1 = mean(SNR1);
    
    SNR2 = [];
    for i=1:s_frames
        SNR2 = [SNR2, 10*log10(1/frame*sum(power(y_2((i-1)*FFT_coeff+1:i*FFT_coeff),2)))-noise_est]; %SNR for one frame of signal 
    end 
    SNR2 = mean(SNR2);
    
    %% Plot
%     [s1,s2] = size(s_w);
%     [n1,n2] = size(n_w);
%     s_w = reshape(s_w,1,s1*s2);
%     n_w = reshape(n_w,1,n1*n2);
%     subplot(3,1,1)
%     plot(s_w)
%     title('original f16 Noise signal');
%     xlabel('time ms');
%     ylabel('magnitude');
%     subplot(3,1,2)
%     plot(n_w)
%     title('noise signal');
%     xlabel('time ms');
%     ylabel('magnitude');
%     subplot(3,1,3)
%     plot(y_2)
%     title('Denoised f16 Noise signal');
%     xlabel('time ms');
%     ylabel('magnitude');
end 