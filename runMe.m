% Copyright: Written by Sijia Zhao, University of Oxford, 2020

addpath(genpath('Functions'));
% Load sound files from the folder "input"
path_in = fullfile('input'); % where is the data?
fileList = dir(fullfile(path_in,'*.wav'));
numSound = length(fileList);

T = table;

for s = 1:numSound
    
    filename = fullfile(fileList(s).folder,fileList(s).name);
    [signal, fs] = audioread(filename);
    
    signal = mean(signal,2);
    
    %% Compute ERB based loudness
    % Use 6 functions written by Alain de Cheveigné
    %compute the ERB-based loudness over the entire waveform (because loudness is one of the most obvious salience factors)
    erb = mean2(ERBpower(signal,fs).^0.3);
    
    %% Compute roughness
    % compute roughness of the entire waveform in the exact same way I did in the 2019 jneurosci paper. I strongly suggest this, as it has been found to be strongly correlated with subjective salience as well as the physiological index of auditory salience.
    filter.method = 3;           % 1: Notch Filter/ 2: High pass/ 3: Low pass/ 4: No filter
    filter.fband = 32;           % fband is the width of the frequency band - use 32 Hz for speech - 125 Hz for zebra finch song
    filter.wf_high = Inf;      % Frequencial cutoff frequency
    filter.wt_high = 30;       % Temporal cutoff frequency
    filter.wf_it = Inf;          % Width of the spectral mask (Notch filter)
    filter.wt_it = Inf ;          % Width of the temporal mask (Notch filter)
    
    filter.song_path = filename;
    filter.mod_song_name = '18Inf30Hz';
    filter.mod_song_path= ['./FilteredFiles/18Inf30Hz'] ;
    
    %     signal = signal(:,1); % signal should be an Nx1 array
    [roughness,~,~] = mean_mps(signal, fs, filter.fband, filter.method, filter.wf_high, filter.wt_high, filter.wf_it, filter.wt_it, filter);
    
    %% Use Kayser's model to compute salience map
    % Require Kayser's salience model
    % extract the key features from Kayser et al.'s salience model output. Briefly speaking, the Kayser model turns a sound into a salience map (looks like the attached image). You can decrease the dimension of this salience map by computing its key features, such as the max amplitude (how warm is the warmest colour in that image), mean amplitude, max change in time, mean change in time, max change in amplitude etc etc.
    
    %     nfft = min(length(signal),256);
    %     window = hanning(nfft);
    %     noverlap = ceil(length(window)/2);
    
    nfft = 1024;
    window = 800;
    noverlap = 778;
    
    dur = length(signal)/fs;
    
    [signalK,f,t1] = specgram(signal,nfft,fs,window,noverlap); % compute spectrogram of this sound
    signalK =log(abs(signalK)); % make intensity map
    SALIENCY = Saliency_map(signalK,4); % compute saliency map
    SAL = SALIENCY.eo + SALIENCY.esi + SALIENCY.epi; % combine saliency maps from the three different filters
    
    % % %     % Plot Kayser salience map for this sound
    % % %     figure(1);clf;
    % % %     subplot(2,1,1);
    % % %     spectrogram(signal,'yaxis');
    % % %     %     axis square;
    % % %     title(fileList(s).name);
    % % %     subplot(2,1,2);
    % % %
    % % %     imagesc(downsample(t1,2),downsample(f,2),SAL);
    % % %     xlabel('Time(s)');
    % % %     set(gca,'YDir','normal');
    % % %     %     axis square;
    
    
    % Extract key stats of the Kayser salience map
    maxSal = max(SAL,[],'all');
    meanSal = mean2(SAL);
    x = mean(SAL);
    maxGradient = max(diff(x)); %max change over time
    meanGradient = mean(diff(x));
    
    T = vertcat(T,...
        cell2table({fileList(s).name, erb, roughness, maxSal,meanSal,maxGradient,meanGradient},...
        'VariableNames',{'filename','ERB','roughness','maxSalience','meanSalience','maxGradient','meanGradient'}));
    
end
writetable(T,'output.csv');

figure(1);clf;
histogram(T.ERB);
