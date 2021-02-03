function [roughness,o1,o2] = mean_mps(sound_in, samprate, fband, method, wf_high, wt_high, wf_it, wt_it, filter)
% % Filters the sound by filtering the modulation spectrum of the sound.
% % Sound_in is the sound pressure waveform.
% % samprate is the sampling rate
% % fband is the width of the frequency band - use 32 Hz for speech - 125
% % Hz for zebra finch song
% % sound_out is the filtered sound
% % method 1 is notch filter where wf_high is the lower bound for spectral
% modulations in cycles/Hz and wf_high+wf_it is the upper bound in
% cycles/Hz (similarly for wt in Hz).
% For a spectral notch, wt filtering limits extend from 0 to at least 100.
% For a temporal notch, wf filtering limits extend from 0 to at least .016.
% Taffeta used 500 cyc/Hz for wt_it on spectral notches; 0.032 Hz wf_it on
% temporal notches.
% % method 2 is high-pass filter where wf_high and wt_high are the cutoff
% frequencies (higher bound of the filter).  Set wf_it and wt_it to 0.
% % method 3 is low-pass filter where wf_high and wt_high are the cutoff
% frequencies (lower bound of the filter).
% % method 4 is unfiltered for control purposes
% % filter can be left blank: it is used by makemodfilterfiles


%% Some parameters: these could become flags
nstd = 6;
twindow = 1000*nstd/(fband*2.0*pi);        % Window length in ms - 6 times the standard dev of the gaussian window
winLength = fix(twindow*samprate/1000.0);  % Window length in number of points
winLength = fix(winLength/2)*2;            % Enforce even window length
increment = fix(0.001*samprate);           % Sampling rate of spectrogram in number of points - set at 1 kHz
ampsrate = samprate/increment;
f_low=250;                                 % Lower frequency bounds to get average amplitude in spectrogram
f_high=10000;                              % Upper frequency bound to get average amplitude in spectrogram
sil_len=500;                               % Amount of silence added at each end of the sound and then subracted
tone_ramp=25;                              % tone ramp in ms

debug_fig = 0;                             % Set to 1 to see spectrograms, 2 to see just resulting spectra, 0 to run many files without creating figures.
save_fig = 0;                              % Set to 1 to save data to recreate images
data_path = './DataFiles/';                % Specify folder if you wish to save data
no_it = 20;                                % The number of iterations for the spectrum inversion
logflg = 1;                                % Perform the mod spectrum in log or linear amplitudes
DBNOISE = 80;                              % dB in Noise for the log compression - values below will be set to zero.
% This defines the ramp of the gain from 0 to 1
dfi=0.0001;                                % 0.1 cycle per kHz ramp in frequency
dti=1;                                     % One Hz ramp in time


soundlen = length(sound_in);
% Force sound_in to have zero mean.
sound_in = sound_in - mean(sound_in);

% fprintf(1,'Length of sound is %f (ms)\n', soundlen*1000.0/samprate);

% find the length of the spectrogram and get a time label in ms
maxlenused = soundlen+fix(sil_len*2.0*(samprate/1000.0));
maxlenint = ceil(maxlenused/increment)*increment;
w = hamming(maxlenint);
frameCount = floor((maxlenint-winLength)/increment)+1;
%t = 0:frameCount-1;
%t = t + (winLength*1000.0)/(2*samprate);

% Pad the sound with silence (silence added at beginning of sound)
input = zeros(1,maxlenint);
nzeros = fix((maxlenint - soundlen)/2);
input(1+nzeros:nzeros+soundlen) = sound_in;

% Get the spectrogram
% input = input .* w';
[s, to, fo, pg] = GaussianSpectrum(input, increment, winLength, samprate);  % Gaussian Spectrum called here to get size of s and fo
fstep = fo(2)-fo(1);
fl = floor(f_low/fstep)+1;        % low frequency index to get average spectrogram amp
fh = ceil(f_high/fstep)+1;        % upper frequency index to get average spectrogram amp
sabs = abs(s);
%% TEST
% sabs(:,1:27)=0;
% sabs(:,end-27:end)=0;
% sabs(1:27, :)=0;
% sabs(end-27:end, :)=0;
%%
sphase = s;
%clear s;
maxsabs = max(max(sabs));

sabsDisp = 20*log10(sabs./maxsabs)+DBNOISE;
sabsDisp(find(sabsDisp<0.0)) = 0.0;

if logflg
    sabs = sabsDisp;
end

nb = size(sabs,1);
nt = length(to);

% variables to save figure data
if save_fig
    toInitial = to;
    foInitial = fo;
end



% Display the spectrogram
if (debug_fig==1)
    figure;
    imagesc(to, fo, sabsDisp, [0 DBNOISE]);
    axis xy;
    cmap = spec_cmap();
    colormap(cmap);
end


% Find the mean level of the amplitude envelope as a function of frequency
meanf = mean(sabs');
meant = mean(sabs);

% calculate the 2D fft
fabs = fft2(sabs);      % Modulation domain

% calculate amplitude and phase
amp_fabs = abs(fabs);
phase_fabs = angle(fabs);

% Calculate the labels for temporal and spectral frequencies in physical
% units
% f_step is the separation between frequency bands
fstep = fo(2);

for ib=1:ceil((nb+1)/2)
    dwf(ib)= (ib-1)*(1/(fstep*nb));
    if (ib > 1)
        dwf(nb-ib+2)=-dwf(ib);
    end
end

% ampsrate is the sampling rate of the amplitude enveloppe
for it=1:ceil((nt+1)/2)
    dwt(it) = (it-1)*(ampsrate/nt);
    if (it > 1 )
        dwt(nt-it+2) = -dwt(it);
    end
end

if (debug_fig==1)  % Display the amplitude and phase spectrum
    figure;
    imagesc(fftshift(dwt), fftshift(dwf).*1000, log(fftshift(amp_fabs)));
    title('Amplitude Spectrum');
    axis xy;
    axis([-50 50 -16 16]);
    colorbar;
    figure;
    %imagesc(dwt, dwf.*1000, unwrap(unwrap(fftshift(phase_fabs)),1));
    imagesc(fftshift(dwt), fftshift(dwf).*1000, fftshift(phase_fabs));
    
    title('Phase Spectrum');
    axis xy;
    axis([-50 50 -16 16]);
    cmap = colormap('HSV');
    cmap(1,:) = [1 1 1];
    colormap(cmap);
    colorbar()
end
if save_fig
    dwtInitial = dwt;
    dwfInitial = dwf;
end

% Filter the sound
tic;
if (method==1)   % Notch pass filtering
    
    newamp_fabs = amp_fabs;
    newphase_fabs = phase_fabs;
    gainmap=ones(size(phase_fabs)); % Define a gain by which to multiply the mod spectrum
    clear ib it;
    
    for ib=1:nb
        for it=1:nt
            % Define the regions to set to zero
            if (abs(dwf(ib))>=wf_high) && (abs(dwf(ib))<(wf_high+wf_it)) && (abs(dwt(it))>=wt_high) && (abs(dwt(it))<(wt_high+wt_it))
                newphase_fabs(ib,it) = (rand(size(phase_fabs(ib,it)))-0.5)*2*pi;     % Randomize the phase
                newphase_fabs(1,1) = 0.0;  % The phase of the DC value has to be zero
                gainmap(ib, it) = 0.0;  % Set gain map to zero
            end
            
            % Perform notch filter with smoothed lead into an empty region within the mod spectrum
            if ((abs(dwf(ib))>=(wf_high-dfi)) && (abs(dwf(ib))<=wf_high) && (abs(dwt(it))>=(wt_high-dti)) && (abs(dwt(it))<=(wt_high+wt_it+dti)))
                gainmap(ib,it)=gainmap(ib,it)*cos((((abs(dwf(ib))-(wf_high-dfi))./dfi))*(pi/2))^2;
            end
            
            if ((abs(dwf(ib))>=(wf_high+wf_it)) && (abs(dwf(ib))<=(wf_high+wf_it+dfi)) && (abs(dwt(it))>=(wt_high-dti)) && (abs(dwt(it))<=(wt_high+wt_it+dti)))
                gainmap(ib,it)=gainmap(ib,it)*cos((((abs(dwf(ib))-(wf_high+wf_it+dfi))./dfi))*(pi/2))^2;
            end
            
            if ((abs(dwt(it))>=(wt_high-dti)) && (abs(dwt(it))<=wt_high) && (abs(dwf(ib))>=(wf_high-dfi)) && (abs(dwf(ib))<=(wf_high+wf_it+dfi)))
                if ((abs(dwf(ib))>=(wf_high-dfi)) && (abs(dwf(ib))<=wf_high) && (abs(dwt(it))>=(wt_high-dti)) && (abs(dwt(it))<=(wt_high+wt_it+dti)) || ...
                        (abs(dwf(ib))>=(wf_high+wf_it)) && (abs(dwf(ib))<=(wf_high+wf_it+dfi)) && (abs(dwt(it))>=(wt_high-dti)) && (abs(dwt(it))<=(wt_high+wt_it+dti)) )
                    gainmap(ib,it)=max(gainmap(ib,it),cos((((abs(dwt(it))-(wt_high-dti))./dti))*(pi/2))^2);
                    
                else
                    gainmap(ib,it)=gainmap(ib,it)*cos((((abs(dwt(it))-(wt_high-dti))./dti))*(pi/2))^2;
                end
            end
            
            if ((abs(dwt(it))>=(wt_high+wt_it)) && (abs(dwt(it))<=(wt_high+wt_it+dti)) && (abs(dwf(ib))>=(wf_high-dfi)) && (abs(dwf(ib))<=(wf_high+wf_it+dfi)))
                if ((abs(dwf(ib))>=(wf_high-dfi)) && (abs(dwf(ib))<=wf_high) && (abs(dwt(it))>=(wt_high-dti)) && (abs(dwt(it))<=(wt_high+wt_it+dti)) || ...
                        (abs(dwf(ib))>=(wf_high+wf_it)) && (abs(dwf(ib))<=(wf_high+wf_it+dfi)) && (abs(dwt(it))>=(wt_high-dti)) && (abs(dwt(it))<=(wt_high+wt_it+dti)) )
                    gainmap(ib,it)=max(gainmap(ib,it),cos((((abs(dwt(it))-(wt_high+wt_it+dti))./dti))*(pi/2))^2);
                    
                else
                    gainmap(ib,it)=gainmap(ib,it)*cos((((abs(dwt(it))-(wt_high+wt_it+dti))./dti))*(pi/2))^2;
                end
            end
        end
    end
elseif (method==2)   % High pass filtering
    newamp_fabs = amp_fabs;
    newphase_fabs = phase_fabs;
    gainmap=ones(size(phase_fabs)); % Define a gain by which to multiply the mod spectrum
    for ib=1:nb
        for it=1:nt
            % Define the regions to set to zero gain - first along the
            % wf asix
            if ((abs(dwf(ib)))<wf_high)
                gainmap(ib, it) = 0.0;
                newphase_fabs(ib,it) = (rand(size(phase_fabs(ib,it)))-0.5)*2*pi;     % Randomize the phase
                
            end
            
            if (wf_high~=0)
                if ((abs(dwf(ib))>=wf_high) && (abs(dwf(ib))<=(wf_high+dfi)))
                    gainmap(ib,it)=gainmap(ib,it)*cos((((abs(dwf(ib))-(wf_high-dfi))./dfi))*(pi/2))^2;
                end
            end
            
            % Define the regions to set to zero - along the wt axis
            if ((abs(dwt(it)))<wt_high)
                gainmap(ib, it) = 0.0;
                newphase_fabs(ib,it) = (rand(size(phase_fabs(ib,it)))-0.5)*2*pi;     % Randomize the phase
            end
            
            if (wt_high~=0)
                if ((abs(dwt(it))>=(wt_high) && (abs(dwt(it))<=(wt_high+dti))))
                    gainmap(ib,it)=gainmap(ib,it)*cos((((abs(dwt(it))-(wt_high-dti))./dti))*(pi/2))^2;
                end
            end
        end
    end
    newphase_fabs(1,1) = 0.0;  % The phase of the DC value has to be zero
elseif (method==3)   % Low pass filtering
    newamp_fabs = amp_fabs;
    newphase_fabs = phase_fabs;
    gainmap=ones(size(phase_fabs)); % Define a gain by which to multiply the mod spectrum
    for ib=1:nb
        for it=1:nt
            % Define the regions to set to zero gain - first along the
            % wf asix
            if ((abs(dwf(ib)))>wf_high+dfi)
                gainmap(ib, it) = 0.0;
                newphase_fabs(ib,it) = (rand(size(phase_fabs(ib,it)))-0.5)*2*pi;     % Randomize the phase
                newphase_fabs(1,1) = 0.0;  % The phase of the DC value has to be zero
            end
            
            if (wf_high~=0)
                if ((abs(dwf(ib))>=wf_high) && (abs(dwf(ib))<=(wf_high+dfi)))
                    gainmap(ib,it)=gainmap(ib,it)*cos((((abs(dwf(ib))-wf_high)./dfi))*(pi/2))^2;
                end
            end
            
            % Define the regions to set to zero - along the wt axis
            if ((abs(dwt(it)))>wt_high+dti)
                gainmap(ib, it) = 0.0;
                newphase_fabs(ib,it) = (rand(size(phase_fabs(ib,it)))-0.5)*2*pi;     % Randomize the phase
            end
            
            if (wt_high~=0)
                if ((abs(dwt(it))>=(wt_high) && (abs(dwt(it))<=(wt_high+dti))))
                    gainmap(ib,it)=gainmap(ib,it)*cos((((abs(dwt(it))-wt_high)./dti))*(pi/2))^2;
                end
            end
        end
    end
    newphase_fabs(1,1) = 0.0;  % The phase of the DC value has to be zero
elseif (method == 4) %No filtering
    newamp_fabs = amp_fabs;
    newphase_fabs = phase_fabs;
    gainmap=ones(size(phase_fabs));
    wf_high = 0;
    wt_high = 0;
    wf_it = 0;
    wt_it = 0;
end
toc;



% gainmap = abs(gainmap-1);
% newamp_fabs=gainmap.*newamp_fabs;
% A = log(fftshift(newamp_fabs));
% res = mean(mean(A(find(A>-inf))));


gainmap = abs(gainmap-1);
newamp_fabs1 = gainmap.*newamp_fabs;
A = log(fftshift(newamp_fabs1));
o1 = mean(mean(A(find(A>-inf))));

newamp_fabs2 = abs(gainmap-1).*newamp_fabs;
A = log(fftshift(newamp_fabs2));
o2 =  mean(mean(A(find(A>-inf))));
roughness = o1/o2;

% newamp_fabs1=gainmap.*newamp_fabs;
% newamp_fabs2 = (1-gainmap).*newamp_fabs1;
% 
% o1 = mean(mean(newamp_fabs1));
% o2 = mean(mean(newamp_fabs2));
% roughness = o1/o2;

% mean_mps = roughness;

if (debug_fig==1 || debug_fig==2)  % Display the amplitude and phase spectrum
    figure;
    imagesc(fftshift(dwt), fftshift(dwf).*1000, log(fftshift(newamp_fabs)));
    title('Resultant Amplitude Spectrum');
    axis xy;
    axis([-500 500 -15 15]);
    colorbar;
    figure;
    %imagesc(dwt, dwf.*1000, unwrap(unwrap(fftshift(phase_fabs)),1));
    imagesc(fftshift(dwt), fftshift(dwf).*1000, fftshift(newphase_fabs));
    title('Resultant Phase Spectrum');
    axis xy;
    axis([-50 50 -16 16]);
    cmap = colormap('HSV');
    cmap(1,:) = [1 1 1];
    colormap(cmap);
    colorbar();
end
if save_fig
    savednewamp_fabs = newamp_fabs;
    savednewphase_fabs = newphase_fabs;
end

