%% Basic analysis script that plots cosine tuning curves for each electrode channel
%% Set up directory system
% set direct path to data folder
dataPath = './Data';
addpath(genpath(dataPath));

% designate the data to be analyzed
dataName = 't5_06_02_2021'; 

%% Load data
dat = load([dataPath filesep dataName '.mat']);

% remove block-wise means
blkNums = [dat.blockNum];
uniqueBlkNums = unique(blkNums);
tx_meanRemoved = zeros(size(dat.tx)); % initialize a new matrix to store block-wise mean removed tx features
for b = 1 : length(uniqueBlkNums)
    blkInds = find(blkNums == uniqueBlkNums(b));
    curr_tx = dat.tx(:,blkInds);
    tx_meanRemoved(:,blkInds) = curr_tx - nanmean(curr_tx,2);
end

features = tx_meanRemoved*50; % multiply by 50 to get into seconds (20ms * 50 = 1s)
N = size(features,1); % number of channels

%% Cosine tuning curves
window = [15 35]; % window (in # of 20ms bins) relative to go cue to avg. spikes for tuning curve, e.g., [15 35]->[300 700]ms window after go cue  

% discretize the 2D target space into 8 equally sized wedges for trial-averaging
targetMapping = [157:201; 202:246; 247:291; 292:336; [337:359 0:21]; 22:66; 67:111; 112:156]'; % maps a target's location (position vector's angle wrt x-axis) to a wedge ID'd by column index 

% define the trial condition in the following form:
% [move type, right hand target ID, left hand target ID]
% move type : 1 = unimanual right trial, 2 = unimanual left trial, 3 = bimanual trial
% target ID is the corresponding column index of targetMapping
conditions = zeros(numel(dat.goCue),3); % initialize n-dimensional array
for trial = 1 : numel(dat.goCue)
    curr_goCue = dat.goCue(trial); % get current trial's go cue
    currTarg = dat.Tp(:,curr_goCue); % get position of both targets 
    currCursor = dat.Cp(:,curr_goCue); % get position of both cursors

    positionError = currTarg - currCursor; % compute the vector pointing from each cursor to its target
    TargID_rightHand = whichTarget_bimanual(targetMapping,positionError(1:2,:)'); % get the target ID for the right hand
    TargID_leftHand = whichTarget_bimanual(targetMapping,positionError(3:4,:)'); % get the target ID for the left hand

    mvType = dat.moveType(trial); % move type can be pulled from data structure
    if mvType == 1 % if a unimanual right hand trial, force left hand target ID to 9
        TargID_leftHand = 9;
    elseif mvType == 2 % do the same for the right hand if uni left trial
        TargID_rightHand = 9;
    end
    conditions(trial,:) = [mvType TargID_rightHand TargID_leftHand];
end

% tile features matrix so that we can easily access windows of activity for
% each desired factorization
featureVals = computeFeatureMatrices(features', conditions(:,2:3), dat.goCue, window);
% featureVals : N x RH x LH x T x Trials

% collect all time bins in each window for all 8 movement directions
% unimanual case only
uniR_tc = [];
uniL_tc = [];
for d = 1 : 8 % iterate over movement directions
    curr_uniR_featvals = squeeze(featureVals(:,d,9,:,:)); % extract data associated with unimanual right trials only (F2 of 9 means 'no movement' for left hand)
    curr_uniL_featvals = squeeze(featureVals(:,9,d,:,:)); % do the same for unimanual left trials (F1 of 9 means 'no movement' for right hand)

    % tile all bins
    tiled_uniR = [];
    tiled_uniL = [];
    for t = 1 : size(curr_uniR_featvals,3)
        tiled_uniR = [tiled_uniR squeeze(curr_uniR_featvals(:,:,t))];
        tiled_uniL = [tiled_uniL squeeze(curr_uniL_featvals(:,:,t))];
    end

    % compute the mean across all collected time bins 
    uniR_tc = [uniR_tc nanmean(tiled_uniR,2)];
    uniL_tc = [uniL_tc nanmean(tiled_uniL,2)];
end
% mean-subtract to center about zero
uniR_tc = uniR_tc - mean(uniR_tc,2); % uniR_tc = unimanual right hand tuning curve 
uniL_tc = uniL_tc - mean(uniL_tc,2); % uniL_tc = unimanual left hand tuning curve

% now for the bimanual context
inds_bi = find(conditions(:,1) == 3);
goCues_bi = dat.goCue(inds_bi); % get go cues for the bimanual trials only
conds_bi_R = conditions(inds_bi,2); % extract movement conditions for the right hand during bimanual context
featureVals_biR = computeFeatureMatrices(features', conds_bi_R, goCues_bi, window);

% extract movement conditions for the left hand during bimanual context
conds_bi_L = conditions(inds_bi,3);
featureVals_biL = computeFeatureMatrices(features', conds_bi_L, goCues_bi, window);

% collect all time bins in each window for all 8 movement directions
biR_tc = [];
biL_tc = [];
for d = 1 : 8
    curr_biR_featvals = squeeze(featureVals_biR(:,d,:,:));
    curr_biL_featvals = squeeze(featureVals_biL(:,d,:,:));

    % tile all 20ms bins
    tiled_biR = [];
    for t = 1 : size(curr_biR_featvals,3)
        tiled_biR = [tiled_biR squeeze(curr_biR_featvals(:,:,t))];
    end

    tiled_biL = [];
    for t = 1 : size(curr_biL_featvals,3)
        tiled_biL = [tiled_biL squeeze(curr_biL_featvals(:,:,t))];
    end

    % average across all collected time bins
    biR_tc = [biR_tc nanmean(tiled_biR,2)];
    biL_tc = [biL_tc nanmean(tiled_biL,2)];
end
% mean-subtract to center about zero
biR_tc = biR_tc - mean(biR_tc,2); % BiR_tc = bimanual right hand tuning curve
biL_tc = biL_tc - mean(biL_tc,2); % BiL_tc = bimanual left hand tuning curve

%% Plot tuning curve for each channel
numPerPage = 10;
rowInds = repmat([1:2:numPerPage*2],1,ceil(N/numPerPage));
figIncs = [1 : numPerPage: N]; % array to help define a new figure when page fills up
degree_order = {'180'   '225'   '270'   '315'     '0'    '45'    '90'   '135'}; % x-axis as seen in Fig.1b of journal paper
for c = 1 : N % iterate over each channel
    if ismember(c,figIncs) % decide whether a new page is needed
        figure('Position',[1346 45 335 902]);
    end
    condCol = [1 0 0; 0 0 1; 1 0 1; 0 1 1];

    % plot right hand tuning curves for both unimanual and bimanual contexts
    subplot(numPerPage,2,rowInds(c));
    plot(uniR_tc(c,:), 'r*-');
    hold on;
    plot(biR_tc(c,:), 'mo-');
    title('RH');
    xticks([1:8]);
    xticklabels(degree_order);
    ylabel(sprintf('CH. %s',num2str(c)));

    % plot left hand tuning curves for both unimanual and bimanual contexts
    subplot(numPerPage,2,rowInds(c)+1);
    plot(uniL_tc(c,:), 'b*-');
    hold on;
    plot(biL_tc(c,:), 'co-');
    title('LH');
    xticks([1:8]);
    xticklabels(degree_order);
end


%% -----------------------------------------------------   Helper functions

function TargID = whichTarget_bimanual(Targets, X)
    % X : Tx2 array of x and y-coordinates (position vectors) for each trial, T is the number of trials 
    % Targets is the target map as described on line 31
    % This function computes the angle (relative to origin) of each 2D position vector with the x -axis
    % Angles will range from 0 - 359 deg, and will return the associated target ID as governed by the target mapping 
    % TargID is returned, which is an array of the Target IDs for each position vector
    
    TargID = [];
    for i = 1 : size(X,1)
        % compute the angle for the current position vector
        theta = floor(angleFromXaxis(X(i,1), X(i,2)));
    
        if (X(i,1) == 0 ) || (X(i,2) == 0) % catch any position vectors without a x- or y-dimension
            TargID(i) = NaN;
        else
            for j = 1 : size(Targets, 2)
                if ismember(theta, Targets(:,j))
                    TargID(i) = j; % check the target map for which wedge the current theta resides
                end
            end
        end
    end
end


%% 
function angle = angleFromXaxis(x,y)
% Computes the angle (in degrees) from the x-axis of the position vector [x y]
    u = [1 0 0];
    v = [x y 0];
    
    norm_u = norm(u);
    norm_v = norm(v);
    
    dotprod = dot(u,v);
    
    theta = acosd(dotprod/(norm_u*norm_v));
    
    angle = 0;
    if (x >= 0 && y >=0) || (x <=0 && y>=0) 
        angle = theta;
    else
        angle = 360-theta;
    end
    
    if (x==0) && (y == 0)
        angle = 0;
    end
end



%%
function featureVals = computeFeatureMatrices(features, trialCodes, eventIdx, timeWindow)
    % This function tiles the features matrix into a N x F1 x F2 x W x R tensor
    % N : # of channels
    % F1 : # of conditions for Factor 1 (in our case 8 movement directions with an additional 'No movement' condition for the right hand)
    % F2 : # of conditions for Factor 2 (same as factor 1 but for the left hand)
    % W : # of timesteps within the timeWindow array
    % R : # of repetitions for each condition set
    
    % Inputs: 
    % features : TxN feature matrix, T is number of time bins, N is # of channels
    % trialCodes: Kx2 array of factor conditions for each trial (K),  column 1 is Factor 1 conditions, column 2 is Factor 2 conditions 
    % eventIdx : Kx1 array of event indices (go cues) where a window (timeWindow) relative to these event indices will be extracted 
    % timeWindow : 1x2 array of timesteps relative to the eventIdx that will be extracted, first element is the start index, second element is the end index

    N = size(features,2);                           % number of features
    T = length(timeWindow(1):timeWindow(2));        % number of time steps in a trial
    nFactors = size(trialCodes,2);                  % number of factors

    % find unique trial codes
    [codeList,~,muxCodes] = unique(trialCodes, 'rows');
    maxRep = max(hist(muxCodes, max(muxCodes)));
    nCodes = size(codeList,1);
    nCons = zeros(nFactors, 1);
    for f=1:nFactors
        nCons(f) = length(unique(trialCodes(:,f)));
    end
    
    % initialize tensors 
    matrixSize_singleTrial = [N, nCons', T, maxRep];
    matrixSize_numTrials = [N, nCons'];
    featureVals = nan(matrixSize_singleTrial); 

    % extract appropriate windows for each trial and tile into the respective factors
    for codeIdx = 1:nCodes
        trlIdx = find(muxCodes==codeIdx);
        indOp = '(:';
        for f=1:nFactors
            indOp = [indOp, ',' num2str(codeList(codeIdx,f))];
        end
        indOp = [indOp, ')'];
        eval(['trialNum' indOp ' = length(trlIdx);']);
        
        indOp_neural = '(:';
        for f=1:nFactors
            indOp_neural = [indOp_neural, ',' num2str(codeList(codeIdx,f))];
        end
        indOp_neural = [indOp_neural, ',:,e)'];
        
        for e = 1:length(trlIdx)
            loopIdx = (eventIdx(trlIdx(e))+timeWindow(1)):(eventIdx(trlIdx(e))+timeWindow(2));
            if loopIdx(end)>size(features,1) || loopIdx(1)<1
                eval(['trialNum' indOp ' = trialNum' indOp ' - 1;']);
                continue;
            end
            
            eval(['featureVals' indOp_neural ' = features(loopIdx,:)'';']);
        end
    end

end

