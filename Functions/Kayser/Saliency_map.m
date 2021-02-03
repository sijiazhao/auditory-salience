function SALIENCY = Saliency_map(img,N)

% function SALIENCY = Saliency_map(img,N)
%
% compute a saliency map of the input spectrogram using a  N-level pyramid.
% as described in Kayser et al CurrBiol 2005
%
% The following types of filters are used:
%   - EO. excitatory-only filters
%   - ESI. excitatory plus inhibitory side-flanks 
%   - EPI. excitatory plus  post-inhibition
%
% output: SALIENCY.eo, SALIENCY.esi, SALIENCY.epi
% the final saliency map is obtained as the sum: 
% SAL = SALIENCY.eo + SALIENCY.esi + SALIENCY.epi

% important local variables
% IMG contains the spectrograms sampled at different scales
% PYR is the feature pyramid. 
% LEVEL 1 contains the filtered images (feature maps)
%       2 the first order center surround interaction
%
% Copyright: Written by C. Kayser, Max Planck Institute for Biological
% Cybernetics, Tuebingen, 2005

% parameters  ------------------------------------------------------------------------
% window for temporal normalization
WINNORM = 150;
% the range of scales at which center surround interactions are computed
scale_interact =[1,2];
% Use on_off and off_on maps?
ONOFF = [1,0];

% prepare images ------------------------------------------------------------------

% insert an image border to avoid edge artifacts
img2 = replicate_img(img);
% resample spectrogram on different scales
for n=1:N
  IMG{n} = imresize(img2,1/(2^(n-1)),'nearest');
  % IMG{n} = imresize(img2,1/(2^(n-1)));
end
% the size at which the maps are stored
MapSize = size(img)/2;

% extract the features on first order --------------------------------------
LEV = 1;  % the level of the pyramid. 1 - filter responses
for n=1:N
  % EO) The intensity feature
  Filter = auditory_RF(0,0);
  PYR.data{1}{LEV}{n} = abs(conv2(IMG{n},Filter,'same'));
  PYR.data{1}{LEV}{n} = imresize(cropimg(PYR.data{1}{LEV}{n}),MapSize,'nearest');
  % ESI) the contrast features
  Filter = auditory_RF(1,0);
  PYR.data{2}{LEV}{n} = abs(conv2(IMG{n},Filter,'same'));
  PYR.data{2}{LEV}{n} = imresize(cropimg(PYR.data{2}{LEV}{n}),MapSize,'nearest');
  % EPI) post inhibition
  Filter = auditory_RF(0,1);
  PYR.data{3}{LEV}{n} = abs(conv2(IMG{n},Filter,'same'));
  PYR.data{3}{LEV}{n} = imresize(cropimg(PYR.data{3}{LEV}{n}),MapSize,'nearest');
  % report if pyramid is too high (i.e. if filters get to big)
  if sum(size(Filter)>size(IMG{n}/2))
    fprintf('Filter size has reached half the data size\n');
    fprintf('%d and %d\n',size(Filter,1),size(IMG{n},1)/2);
  end
end

% center surround interactions on different spatial scales ---------------
PYR = CenterSurroundPyramid(LEV,PYR,scale_interact,N,ONOFF);

% now fuse the different spatial scales within each pyramid ------------
% first use the normalization of each map with respect to local and
% global maxima. Then  sum the different maps across spatial scales.

LEV = 2;
SALIENCY.eo = zeros(size(PYR.data{1}{LEV}{1}));
SALIENCY.esi = zeros(size(PYR.data{2}{LEV}{1}));
SALIENCY.epi = zeros(size(PYR.data{3}{LEV}{1}));
for O=find(ONOFF)
  % EO)
  for n=1:size(PYR.data{1}{LEV},2)
    SALIENCY.eo = SALIENCY.eo + normalizemap(PYR.data{1}{LEV}{O,n},WINNORM);
  end
  % ESI
  for n=1:size(PYR.data{2}{LEV},2)
    SALIENCY.esi = SALIENCY.esi + normalizemap(PYR.data{2}{LEV}{O,n},WINNORM);
  end
  % EPI)
  for n=1:size(PYR.data{2}{LEV},2)
    SALIENCY.epi = SALIENCY.epi + normalizemap(PYR.data{3}{LEV}{O,n},WINNORM);
  end
end

SALIENCY.eo = SALIENCY.eo / (n*length(find(ONOFF)));
SALIENCY.esi= SALIENCY.esi / (n*length(find(ONOFF)));
SALIENCY.epi = SALIENCY.epi / (n*length(find(ONOFF)));

return;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function PYR = CenterSurroundPyramid(LEVin,PYR,offsets,N,ONOFF)
% center surround interactions on different spatial scales
% offsets: indicates the differences in spatial scale used for the
%          interaction
LEV = LEVin+1;
% for each feature map
cnt=1;
for n1=1:N
  for n2=n1+offsets
    if n2<=N
      % loop the features
      for f=1:length(PYR.data)
        if ONOFF(1)
          % on-off
          dummy = PYR.data{f}{LEVin}{n1}-PYR.data{f}{LEVin}{n2};
          PYR.data{f}{LEV}{1,cnt} =  dummy.*(dummy>0);
          PYR.helper{f}{LEV}(1,cnt) = abs(n2-n1);
        end
        if ONOFF(2)
          % off-on
          dummy = PYR.data{f}{LEVin}{n2}-PYR.data{f}{LEVin}{n1};
          PYR.data{f}{LEV}{2,cnt} =  dummy.*(dummy>0);
          PYR.helper{f}{LEV}(cnt) = abs(n2-n1);
        end
      end
      cnt=cnt+1;
    end
  end
end

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = cropimg(img)
% crop the image border. Undo replicate_img
% the input consists of 0.5,1,0.5 times the real image.
s = ceil(size(img)/4);
out = img(s(1)+1:end-s(1)+1,s(2)+1:end-s(2)+1);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = replicate_img(img)
% replicate the image to its borders to later avoid
% edge artifacts. This is undone by cropimg
out = zeros(size(img,1)*3,size(img,2)*3,size(img,3));
S = size(img);
if length(S)==2
  S(3) = 1;
end

for s=1:S(3)
  % top
  k=1; out([1:S(1)],[1:S(2)]+(S(2)*(k-1)),s) = fliplr(flipud(img(:,:,s)));
  k=3; out([1:S(1)],[1:S(2)]+(S(2)*(k-1)),s) = fliplr(flipud(img(:,:,s)));
  k=2; out([1:S(1)],[1:S(2)]+(S(2)*(k-1)),s) = (flipud(img(:,:,s)));
  % bottom
  k=2;  out([1:S(1)]+(S(1)*2),[1:S(2)]+(S(2)*(k-1)),s) = flipud(img(:,:,s));
  k=1;  out([1:S(1)]+(S(1)*2),[1:S(2)]+(S(2)*(k-1)),s) = fliplr(flipud(img(:,:,s)));
  k=3;  out([1:S(1)]+(S(1)*2),[1:S(2)]+(S(2)*(k-1)),s) = fliplr(flipud(img(:,:,s)));
  % left right
  k=1;  out([1:S(1)]+(S(1)),[1:S(2)],s) = fliplr(img(:,:,s));
  k=3;  out([1:S(1)]+(S(1)),[1:S(2)]+(S(2)*(k-1)),s) = fliplr(img(:,:,s));
  out([1:S(1)]+S(1),[1:S(2)]+S(2),s) = img(:,:,s);
end
out = out(ceil(S(1)/2):end-ceil(S(1)/2)+1,ceil(S(2)/2):end-ceil(S(2)/2)+1,:);
return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = normalizemap(in,win)
% function that normalizes a feature map with respect to 
% local maxima. In the visual case, this is a glocal spatial operation. 
% Here we need an interaction global in frequency but more localized in temporal
% dimension.
% use a sliding window analysis that normalizes windows of length win
% independently, but taking the data from a 3*win window into account.

warning off
ntp = size(in,2);
winborder = [1:win:ntp];
out = zeros(size(in));
if mod(ntp,win)~=0
  winborder = [winborder,ntp];
end

% before normalizing, apply a mask to avoid edge-effects
scale = 11;
h = hanning(scale*2);
mask1 = ones(size(in));
mask2 = mask1;
mask1(1:scale,:) = h(1:scale)*ones(1,size(mask1,2));
mask1(end-scale+1:end,:) = h(scale+1:end)*ones(1,size(mask1,2));
mask2(:,1:scale) = ones(size(mask1,1),1)*(h(1:scale)');
mask2(:,end-scale+1:end) = ones(size(mask1,1),1)*(h(scale+1:end)');

in = in.*mask1.*mask2;
% get the local extrema
in = in-min(in(:));
in = in./max(in(:));
[LocMa,LocMi] = localextrema(in);

for W=1:length(winborder)-1
  I_win = [winborder(W):winborder(W+1)];
  I_all = round([winborder(W)-win:winborder(W+1)+win/3]);
  % clip this interval to the data range
  I_all = I_all(find(I_all>0));
  I_all = I_all(find(I_all<=(ntp-2)));
  data = in(:,I_all);
  max_data = max(data(:));
  data = data./max_data;
  globmax = max(data(:));
  
  LocMaX = find(LocMa(:,I_all));
  LocMaX = data(LocMaX(:));
  % cancel the global maxima from this list.
  LocMaX = LocMaX(find((LocMaX~=1)));
  LocMaX = LocMaX./max_data;
  % now normalize the data in the smaller interval
  if isempty(LocMaX)
    % no local maximum left. Don't normalize but report
    fprintf('problem with normalization\n');
    keyboard;
    out_n = (in(:,I_win))./max_data;
  else
    out_n = (in(:,I_win))./max_data;
    out_n = out_n*((1-mean(LocMaX))^2);
  end
  out(:,I_win) = out_n;
end
warning on;
return;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function out = auditory_RF(sideband,postinh)

% function out = auditory_RF(LAT,BW,sideband,postinh)
%
% Computes and auditory RF model consisting of a Gabor, with or
% without sidebands in frequency direction and a possible post-inhibition.
% the parameters are:
%
% LAT = [Latency Gabor, Latency post-inh]
% BW  = [Bandwidth Gabor, Bandwidth post-inh]
% sideband if 1 sidebands appear
% postinh if 1 a post-inhibition is added
%
% time in 1 ms steps and frequency in quarter tones!

% the parameters
LAT = [24,50];
BW = [0.035,0.04];
Freq = 1;
DUR = 5;

if sideband
  BW(1) = 0.08;
else
  BW(2) = 0.035;
end

% this is the frequency vector
OCTAVES = round(1.3458*2.^[5:1/8:13]);

% time in 1 ms steps and frequency in quarter tones!
Fax = [0.3:1/32:1.7-0.000001];
Tax = [1:80];
[T,F] = meshgrid(Tax,Fax);
LAT = LAT;
BW = BW*4;

% the gabor
env = exp( -( ((T-LAT(1)).^2)/(2*DUR^2) + ((F-1).^2)/(2*BW(1)^2)));
osc = cos(2*pi*F*Freq);
out = osc.*env;

if postinh
  DUR2 = DUR*1.3;
  env = exp( -( ((T-LAT(2)).^2)/(2*DUR2^2) + ((F-1).^2)/(2*BW(2)^2)));

out = out - osc.*env/2;  
end
out = out/abs(sum(out(:)));

return;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Maxima,Minima] = localextrema(in)

% function [Maxima,Minima] = localextrema(in)
% 
% find the local maxima and minima of a function

d1 = diff(in,1,1);
d2 = diff(in,1,2);
% find zeros by detecting zero crossings
dum1 = d1(1:end-1,:);
dum2 = d1(2:end,:);
% downwards
cross11 = (dum1>0).*(dum2<=0);
% upward crossing
cross21 =(dum1<0).*(dum2>=0);
dum1 = d2(:,1:end-1);
dum2 = d2(:,2:end);
% downwards
cross12 =(dum1>0).*(dum2<=0);
% upward crossing
cross22 = (dum1<0).*(dum2>=0);

s = size(in)-2;
cross11 = cross11([1:s(1)],[1:s(2)]);
cross12 = cross12([1:s(1)],[1:s(2)]);
cross21 = cross21([1:s(1)],[1:s(2)]);
cross22 = cross22([1:s(1)],[1:s(2)]);

cross11 = [zeros(s(1),1),cross11];
cross12 = [zeros(s(1),1),cross12];
cross11 = [zeros(1,s(2)+1);cross11];
cross12 = [zeros(1,s(2)+1);cross12];

cross21 = [zeros(s(1),1),cross21];
cross22 = [zeros(s(1),1),cross22];
cross21 = [zeros(1,s(2)+1);cross21];
cross22 = [zeros(1,s(2)+1);cross22];

% a local maximum occurs if both derivatives cross zero downwards
Maxima = cross11.*cross12;
Minima = cross21.*cross22;



