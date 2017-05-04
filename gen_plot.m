function gen_plot()
    
%% Part 2: Weighted Noise cancellation system
% Modify the original noise cancellation system to improve reconstruction
% quality
   a = 2.8
    b = 0.8;
    
    %% locate the file
    filepath1 = '/Users/jefflai108/signals_systems/project2/section 3/whiteNoise_sample.wav';
    filepath2 = '/Users/jefflai108/signals_systems/project2/section 3/violetNoise_sample.wav';
    filepath3 = '/Users/jefflai108/signals_systems/project2/section 3/pinkNoise_sample.wav';
    filepath4 = '/Users/jefflai108/signals_systems/project2/section 3/greyNoise_sample.wav';
    filepath5 = '/Users/jefflai108/signals_systems/project2/section 3/f16Noise_sample.wav';
    Tms = 20; % ms per frame 
    [x1,Fs] = audioread(filepath1); % (40713,8000)
    [x2,Fs] = audioread(filepath2);
    [x3,Fs] = audioread(filepath3);
    [x4,Fs] = audioread(filepath4);
    [x,Fs] = audioread(filepath5);
    frame = Tms*Fs/1000;  % samples per frame 
    FFT_coeff = 160;
    
    %% segment the signal
    % The input signal x(t) is analyzed over short segments xw(t) of size Tw
    n_iteration = ceil(length(x)/frame); % total number of frames 
    x_w = zeros(n_iteration,frame);
    
    for i=1:n_iteration-1
        x_w(i,:) = x((i-1)*frame+1:i*frame);
    end
    
    %% speech/noise classification
    % decide if each frame is speech or noise. The first second of the
    % signal is noise. The rest is speech. 
    n_frames = 1*1000/Tms; % number of noise frames  
    n_w = zeros(n_frames,frame);
    
    for i=1:n_frames 
        n_w(i,:) = x_w(i,:);
    end 
    

    %% FFT for n_w 
    % Magnitude spectrum N_w is the average fft for each noise frames
    N_w = zeros(1,FFT_coeff);
    for i=1:n_frames
        N_w = N_w + abs(fft(n_w(i,:),FFT_coeff)); % update N_w as more frames come in 
    end 
    N_w = N_w/n_frames; % average of all noise frames 
    
    f = 1:50:Fs 
    plot(f,N_w);
    title('f16 noise frequency spectrum');
    xlabel('frequency Hz');
    ylabel('magnitude');
    
%         subplot(3,1,1)
%     plot(y_2)
%     title('denoised signal with alpha = 2.8 beta = 0.1');
%     xlabel('time ms');
%     ylabel('sample data');
%     subplot(3,1,2)
%     plot(y_3)
%     title('denoised signal with alpha = 2.8 beta = 0.5');
%     xlabel('time ms');
%     ylabel('sample data');
%     subplot(3,1,3)
%     plot(y_4)
%     title('denoised signal with alpha = 2.8 beta = 0.8');
%     xlabel('time ms');
%     ylabel('sample data');
    %   subplot(4,1,4)
%     plot(y_2)
%     title('PART 2 denoised signal');
end 