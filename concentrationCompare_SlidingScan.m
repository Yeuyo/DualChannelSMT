% DO NOT COMMENT THESE OUT
clear; clc; close all;

% Folder containing all the tracked files from the dashboard
c1Path = 'G:\Bitong\OCT4 SMT_Mode 1\640';
% Folder containing mask files
c2Path = 'G:\Bitong\OCT4 SMT_Mode 1\Masks';
resultPath = 'G:\Bitong\OCT4 SMT_Mode 1\Result';
% Plot figures?
toPlot = 0;
% Seconds to be defined as long binding
tLong = 3;
distTol = 4;
lI = 8; mI = 13; hI = 16;
fL = 1; % frame to stay long to be defined as long
pxSize = 0.11;
totalFrames = 300;

%% Parameters - BECAREFUL WITH CHANGING VALUES BELOW THIS LINE
analysis_type = "percentage"; % percentage or number for traj used
clip_factor = 100; % 0.8; % percentage or number of tracks in a trajectory trajectory to use for fitting of MSD
traj_length = 1; % Length of traj to keep (traj appear with less than this number of frame will be discarded)
min_traj = 1; % Minimum trajectorys in a file to be accepted into the analysis
tol = 12; % Numbers of decimals to keep for rounding
LocalizationError = -6.5; % Localization Error: -6 = 10^-6
EmissionWavelength = 580; % wavelength in nm; consider emission max and filter cutoff
ExposureTime = 500; % in milliseconds
NumDeflationLoops = 0; % Generaly keep this to 0; if you need deflation loops, you are imaging at too high a density;
MaxExpectedD = 0.3; % The maximal expected diffusion constant for tracking in units of um^2/s;
NumGapsAllowed = 1; % the number of gaps allowed in trajectories

%%% DEFINE STRUCTURED ARRAY WITH ALL THE SPECIFIC SETTINGS FOR LOC AND TRACK
% imaging parameters
impars.PixelSize=0.11; % um per pixel
impars.psf_scale=1.35; % PSF scaling
impars.wvlnth= EmissionWavelength/1000; %emission wavelength in um
impars.NA=1.49; % NA of detection objective
impars.psfStd= impars.psf_scale*0.55*(impars.wvlnth)/impars.NA/1.17/impars.PixelSize/2; % PSF standard deviation in pixels
impars.FrameRate= ExposureTime/1000; %secs
impars.FrameSize= ExposureTime/1000; %secs

% localization parameters
locpars.wn=9; %detection box in pixels
locpars.errorRate= LocalizationError; % error rate (10^-)
locpars.dfltnLoops= NumDeflationLoops; % number of deflation loops
locpars.minInt=0; %minimum intensity in counts
locpars.maxOptimIter= 50; % max number of iterations
locpars.termTol= -2; % termination tolerance
locpars.isRadiusTol=false; % use radius tolerance
locpars.radiusTol=50; % radius tolerance in percent
locpars.posTol= 1.5;%max position refinement
locpars.optim = [locpars.maxOptimIter,locpars.termTol,locpars.isRadiusTol,locpars.radiusTol,locpars.posTol];
locpars.isThreshLocPrec = false;
locpars.minLoc = 0;
locpars.maxLoc = inf;
locpars.isThreshSNR = false;
locpars.minSNR = 0;
locpars.maxSNR = inf;
locpars.isThreshDensity = false;

% tracking parameters
trackpars.trackStart=1;
trackpars.trackEnd=inf;
trackpars.Dmax= MaxExpectedD;
trackpars.searchExpFac=1.2;
trackpars.statWin=10;
trackpars.maxComp=3;
trackpars.maxOffTime=NumGapsAllowed;
trackpars.intLawWeight=0.9;
trackpars.diffLawWeight=0.5;

%% DO NOT CHANGE ANYTHING BELOW THIS LINE
c1List = dir([c1Path, filesep, '*.tif']);
c2List = dir([c2Path, filesep, '*.tif']);

addpath(genpath(['.' filesep 'Batch_MTT_code' filesep])); % MTT & BioFormats
fileTypeList = strings(0);
fileTypes = [];
for n = 1 : length(c1List)
  c1Name = split(c1List(n).name, '_');
  c1Name = [c1Name{1}, '_', c1Name{2}];
  c2Name = split(c2List(n).name, '_');
  c2Name = [c2Name{1}, '_', c2Name{2}];
  if c1Name == c2Name
    % work SMT on c1
    [c1Stack, nbImages] = tiffread([c1Path, filesep, c1List(n).name]);
    
    imgs_3d_double = double(reshape([c1Stack.data], size(c1Stack(1).data, 1), size(c1Stack(1).data, 2), nbImages));
    data = localizeParticles_ASH(0, impars, locpars, imgs_3d_double);
    data=buildTracks2_ASH(0, 0, data, impars, locpars, trackpars, data.ctrsN, imgs_3d_double);
    data_cell_array = data.tr;
    trackedPar = struct;
    for i = 1:length(data_cell_array)
%       trackedPar(1,i).xy =  impars.PixelSize .* data_cell_array{i}(:,1:2);
      trackedPar(1,i).xy =  data_cell_array{i}(:,1:2);
      trackedPar(i).Frame = data_cell_array{i}(:,3);
      trackedPar(i).TimeStamp = impars.FrameRate.* data_cell_array{i}(:,3);
    end
    tracksE = struct('data', []);
    for m = 1 : numel(trackedPar)
      tracksE(m) = struct('data', [trackedPar(m).TimeStamp, trackedPar(m).xy, trackedPar(m).Frame, repelem(m, size(trackedPar(m).xy, 1))']);
    end
    tracksE = cell2mat(reshape(struct2cell(tracksE), [], 1));
    tracksE = [tracksE, zeros(size(tracksE, 1), 2)]; % time, x, y, frame, traj, int, long
    % Debug - check traj localised
    % imtool(c1Stack(1).data)
    % % pause
    % figure;
    % imshow(c1); hold on;
    % ind = find(tracksE(:, 4) == 1);
    % plot(tracksE(ind, 2), tracksE(ind, 3), 'x');
    % End Debug


    % New Pdist approach
    pts = tracksE;
    pts(:, 7) = [1 : size(pts, 1)];
    clusterID = 0;
    idx = find(pts(:, 6) == 0);
    while length(idx) > 1
      clusterID = clusterID + 1;
      D = pdist2(pts(idx, 2:3), pts(idx, 2:3));
      ptsInRange = sum(D <= 4);
      ptsMax = find(ismember(ptsInRange, max(ptsInRange(:))));
      if length(ptsMax) > 1
        row = find(D(:, ptsMax(1)) <= 4);
      else
        row = find(D(:, ptsMax) <= 4);
      end
      tempPts = pts(idx, 7);
      pts(tempPts(row), 6) = clusterID;
      idx = find(pts(:, 6) == 0);
    end

    clusterOccupancy = zeros(clusterID, 1);
    for i = 1 : clusterID
      clusterFrame = pts(find(pts(:, 6) == i), 4);
      clusterOccupancy(i) = length(unique(clusterFrame)) / totalFrames;
    end
    writecell(clusterOccupancy, [resultPath, filesep, c1Name, '_Occupancy.xls'])
    % End Pdist approach

    % break c2 into 3 layers of intensity
    c2Stack = tiffread([c2Path, filesep, c2List(n).name]);
    c2Data = zeros(size(c2Stack(1).data, 1), size(c2Stack(1).data, 2), length(c2Stack)); % intensity in each frame
    c2RH = zeros(size(c2Stack(1).data, 1), size(c2Stack(1).data, 2)); % number of frames in high intensity
    for m = 1 : length(c2Stack)
      intRaw = c2Stack(m).data;
      % intFil = filterStack(intRaw, [150, inf]);
      intFil = intRaw;
      intL = imquantize(intFil, multithresh(intFil, 20));
      % Debug - check traj on high concentration
      % imtool(intL);
      % pause
      % figure;
      % fig = differenceOfGaussian(intRaw, 0.5, 2);
      % imtool(fig);
      % fig = detectFaintBlobs(intRaw, 0.5, 0.2, false);
      % imtool(fig);
      % imshow(c2); hold on;
      % plot(tracksE(ind, 2), tracksE(ind, 3), 'x');
      % End Debug

      ind = sub2ind(size(c2Stack(1).data), round(tracksE(:, 3)), round(tracksE(:, 2)));
      int = intL(ind);
      fInd = find(tracksE(:, 4) == m);
      tracksE(fInd, 6) = int(fInd);

      c2Region = zeros(size(intL, 1), size(intL, 2));
      c2Region(find(intL >= hI)) = 3;
      c2Region(find(intL >= mI & intL <hI)) = 2;
      c2Region(find(intL >= lI & intL <mI)) = 1;
      c2Data(:, :, m) = c2Region;
      c2RH(find(intL >= hI)) = c2RH(find(intL >= hI)) + 1;
    end
    % Quantify how many trajectories are high intensity
    % 16+ is H, 13-15 M, 8+ in cell
    perctH = numel(int(int>=hI)) / length(int);
    highInt = multithresh(intFil, 20);
    highInt(hI)
    highInt(mI)
    perctM = (numel(int(int>=mI)) - numel(int(int>=hI))) / length(int);
    perctL = (numel(int(int>=lI)) - numel(int(int>=mI))) / length(int);
    perctT = perctL + perctM + perctH;
    realH = perctH / perctT;
    realM = perctM / perctT;
    realL = perctL / perctT;

    % Only accounting for the long binding trajectories
    fLong = tLong / (ExposureTime / 1000);
    for m = 1 : length(trackedPar)
      if length(trackedPar(m).Frame) >= fLong
        ind = find(tracksE(:, 5) == m);
        tracksE(ind, 7) = 1;
      end
    end
    tracksL = tracksE(find(tracksE(:, 7) == 1), :);
    int = tracksL(:, 6);
    perctH = numel(int(int>=16)) / length(int);
    highInt = multithresh(intFil, 20);
    highInt(16)
    highInt(13)
    perctM = (numel(int(int>=mI)) - numel(int(int>=hI))) / length(int);
    perctL = (numel(int(int>=lI)) - numel(int(int>=mI))) / length(int);
    perctT = perctL + perctM + perctH;
    realH = perctH / perctT
    realM = perctM / perctT
    realL = perctL / perctT
  else
    % SHOULD NOT HAPPEN
    continue
  end
  % plot for collaborator
%   frameToDraw = 30;
%   imStack = c2Stack(frameToDraw).data;
%   ind = find(c2Data(:, :, frameToDraw) >= 2);
%   imDraw = zeros(size(imStack, 1), size(imStack, 2));
%   imDraw(ind) = imStack(ind);
% %   imtool(imDraw);
%   imtool(c2Data(:, :, frameToDraw));
%   % pause here to export imBack
%   1
%   figure;
%   imshow(imBack); hold on;
%   [col, row] = ind2sub(size(imDraw), ind);
%   frame = find(tracksE(:,1) == frameToDraw);
%   tracksInFrame = tracksE(frame, 2:3);
%   plot(tracksInFrame(:, 1), tracksInFrame(:, 2), 'rx', 'MarkerSize', 9, 'LineWidth', 2);
%   1

  % find the sox2 in high intensity region
  regionToDot = zeros(size(c2Stack(1).data, 1), size(c2Stack(1).data, 2));
  regionToDot(find(c2RH >= fL)) = 1;
%   imshow(regionToDot)
  regionToDot = ~regionToDot;
  D = bwdist(regionToDot);
  figure; imshow(regionToDot); hold on;
  [col, row] = ind2sub(size(regionToDot), find(D>=4));
  plot(row, col, 'rx')
  sox2Appear = zeros(nbImages, 1);
  sox2Coord = [];
  sox2AppearLong = zeros(nbImages, 1);
  sox2Number = zeros(nbImages, 1);
  minDist = zeros(nbImages, 1);
  soxInFrame = [];
  coord = [row + 30, col - 30];
  saveData = zeros(length(col), 2);
  % PAUSE HERE, GET COORDINATES AND PUT IN NEXT LINE
  % col = []; row = []; % col is Y and row is X
  for m = 1 : length(col)
    sox2Appear = zeros(nbImages, 1);
    % figure;
    % imshow(regionToDot); hold on;
    % plot(row(m), col(m), 'rx');
    % q = 0:0.01:2*pi;
    % r = 4;
    % x = r*cos(q);
    % y = r*sin(q);
    % plot(x+row(m), y+col(m), 'r');
    % coord = [row, col];
    distSox = squareform(pdist([coord(m, :); tracksE(:, 2:3)]));
    distSox = distSox(1, 2:end);
    soxNum = length(unique(tracksE(find(distSox < distTol), 5))); % stop display
    soxFrame = length(unique(tracksE(find(distSox < distTol), 4))); % stop display
    for o = 1 : nbImages
      fTime = o * ExposureTime / 1000;
      frame = find(tracksE(:,1) == fTime);
      tracksInFrame = tracksE(frame, 2:3);
      if ~isempty(frame)
        distSox = squareform(pdist([coord(m, :); tracksInFrame]));
        distSox = distSox(1, 2:end);
        minDist(o) = min(distSox);
        if sum(distSox < distTol) > 0
          sox2Appear(o) = 1;
          sox2Number(o) = sum(distSox < distTol);
          sox2Coord = [sox2Coord; repmat(o, size(tracksInFrame(distSox < distTol, :), 1), 1), tracksInFrame(distSox < distTol, :)];
          if size(soxInFrame, 1) > 0
            for p = 1 : size(soxInFrame, 1)
              distSoxH = squareform(pdist([soxInFrame(p, :); tracksInFrame(distSox < distTol, :)]));
              distSoxH = distSoxH(1, 2:end);
              if sum(distSoxH < 1) > 0
                % a sox last longer than 2 frames!
                sox2AppearLong(o) = sox2AppearLong(o) + 1;
              end
            end
          end
          soxInFrame = tracksInFrame(distSox < distTol, :);
        else
          soxInFrame = [];
        end
      end
    end
    % save data
    saveData(m, :) = [soxNum, soxFrame];
    figure;
    plot(1:nbImages, sox2Appear, 'Color', '#BB8FCE', 'linewidth', 1.25); % Sox 2 on last seen GFP spot
    title('Sox2 Occur', 'interpreter', 'latex');
    xlabel('Frame', 'interpreter', 'latex');
    ylabel('Presence of Proteins', 'interpreter', 'latex');
    legend('Sox 2', 'interpreter', 'latex') 
    ylim([0, 1.2]);
    saveas(gcf, [resultPath, filesep, c1Name, '_', num2str(m), '.png']);
    close all;
    1
  end
  cellSummaryData = cell(1, 2);
  cellSummaryData{1} = 'Total number of Sox detected:';
  cellSummaryData{2} = 'Number of frame with Sox:';
  cellSummaryData = [cellSummaryData; arrayfun(@num2str,saveData,'un',0)];
  writecell(cellSummaryData, [resultPath, filesep, c1Name, '_Summary.xls'])
  1
end

%% functions
function stack = filterStack(stackIn, filterRange)
  stackF = stackIn > filterRange(1);
  stackIn = double(stackIn) .* stackF;
  stackF = stackIn < filterRange(2);
  stack = double(stackIn) .* stackF;
end

function dogImage = differenceOfGaussian(image, sigma1, sigma2)
% DIFFERENCEOFGAUSSIAN - Applies Difference of Gaussians to an image
%
% Inputs:
%   image    - Input grayscale image (can be uint8, uint16, or double)
%   sigma1   - Standard deviation of the first Gaussian kernel (smaller)
%   sigma2   - Standard deviation of the second Gaussian kernel (larger)
%
% Output:
%   dogImage - Difference of Gaussians result (double precision)

    % Convert image to double for processing
    % img = im2double(image);
    img = image;
    
    % Ensure sigma1 < sigma2
    if sigma1 >= sigma2
        error('sigma1 must be smaller than sigma2');
    end
    
    % Create Gaussian kernels
    % We use 2*ceil(3*sigma) + 1 as the kernel size for sufficient coverage
    size1 = 2 * ceil(3 * sigma2) + 1;  % Use larger sigma for size
    kernel1 = fspecial('gaussian', size1, sigma1);
    kernel2 = fspecial('gaussian', size1, sigma2);
    
    % Apply Gaussian blurs
    blurred1 = imfilter(img, kernel1, 'replicate');
    blurred2 = imfilter(img, kernel2, 'replicate');
    
    % Compute Difference of Gaussians
    dogImage = blurred1 - blurred2;
    
    % Optional: Normalize to [0,1] for display
    % dogImage = mat2gray(dogImage);
    
end

function blobImage = detectFaintBlobs(img, sigma, thresh, enhance)
% DETECTFAINTBLOBS - Detect faint blobs in low-contrast images
%
% Inputs:
%   img      - Input grayscale image (uint8/double)
%   sigma    - Scale of blobs to detect (e.g., 2–5)
%   thresh   - Detection threshold (small, e.g., 0.001–0.01)
%   enhance  - (Optional) Whether to enhance contrast [true/false]

if nargin < 4, enhance = true; end

% Step 1: Convert to double and enhance contrast if needed
% img = im2double(img);

if enhance
    % Option 1: CLAHE - Contrast Limited Adaptive Histogram Equalization
    img = adapthisteq(img, 'ClipLimit', 0.02);
    
    % Optional: Normalize intensity range
    img = mat2gray(img);
end

% Step 2: Apply Laplacian of Gaussian (LoG) for blob detection
% LoG response is strong at blob centers (positive or negative)
logFilter = fspecial('log', round(6*sigma)+1, sigma);
logImage = imfilter(img, logFilter, 'replicate');

% Step 3: Find local maxima in LoG response (blob centers)
% We look for strong negative or positive peaks
% Since blobs can be bright/dark, use absolute response
absLog = abs(logImage);

% Find local maxima above threshold
blobImage = imregionalmax(absLog) & (absLog > thresh);

% Optional: Clean up with small morphological opening
blobImage = bwareaopen(blobImage, 3);  % Remove tiny detections

% Display results
figure;
subplot(1,3,1); imshow(img); title('Enhanced Image');
subplot(1,3,2); imshow(logImage, []); title('LoG Response');
subplot(1,3,3); imshow(label2rgb(bwlabel(blobImage))); title('Detected Blobs');

% Optional: Mark blob centers on original
hold on;
[y, x] = find(blobImage);
plot(x, y, 'r+', 'MarkerSize', 10, 'LineWidth', 2);
hold off;

end