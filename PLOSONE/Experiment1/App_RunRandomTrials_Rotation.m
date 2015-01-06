clear, clc

% Run test sequence, binning test into different ranges of xfm offsets

isoInit = 0;
anisBothSets = 1;

maxIter = 20;
threshAng = 0.0001;
threshPos = 0.0001;
nTrials = 1000;  % nTrials per bin
numPts = 50;

writeToFile = 1;
outputFile = ['TestResults_TLS-Rotation_isoInit',num2str(isoInit),...
  '_anisBothSets', num2str(anisBothSets),'.txt'];
fid = [];
if (writeToFile)
  fid = fopen(outputFile,'w');
  if fid == -1
    error('Open Output File Failed!')
  end
end

addpath('./Balachandran');
addpath('./mtimesx');

rand('seed',10987)
randn('seed',10988)

% magnitude of rotation offsets
RotBins = [
   0,15
   15,45
   45,90
   90,150
   150,180
  ];
nRotBins = size(RotBins,1);

% % magnitude of translation offsets
% TransBins = [
%   10,20
%   90,100
% ];
% nTransBins = size(TransBins,1);
nTransBins = 1;

nAlgs = 3;

% standard output (all methods)
std_field_names = {...
  'AlgName',...
  'FRE',...
  'numIter',...
  'time',...  
  'ErrAng',...
  'ErrTrans'
  % 'ErrAxis',...
  % 'ErrTransVec',...
  % 'dTheta',...
  % 'dAxis',...
  % 'dTrans'
  };
% preallocate struct for storing iteration data
%  (this data is overwritten by every bin)
empty_cells = repmat(cell({0}),1,numel(std_field_names));
entries = {std_field_names{:}; empty_cells{:}};
trialData = struct(entries{:});
trialData(nAlgs,nTrials).ErrAng = 0;   % initialize struct size
trialData(:) = trialData(1);       % initialize struct data to zero

% custom output Estepar et al.
estepar_field_names = {...
  'numLoopIter',...
  'LoopData',...
  };
% preallocate struct for storing data
empty_cells = repmat(cell({0}),1,numel(estepar_field_names));
entries = {estepar_field_names{:}; empty_cells{:}};
EData = struct(entries{:});
EData(nTrials).ErrAng = 0;   % initialize struct size
EData(:) = EData(1);         % initialize struct data to zero

dispCapture(sprintf('\nLegend:  (numIter runTime FRE ErrAng)'), fid, writeToFile );
%dispCapture(sprintf('\nLegend:  (numIter runTime FRE ErrAng ErrTrans)'), fid, writeToFile );
% dispCapture(sprintf('\nTLS Legend:           (numIter runTime FRE ErrAng ErrTrans)'), fid, writeToFile );
% dispCapture(sprintf('Estepar Legend:       (numIter runTime FRE ErrAng ErrTrans numLoops)\n'), fid, writeToFile );
% dispCapture(sprintf('Balachandran Legend:  (numIter runTime FRE ErrAng ErrTrans)\n'), fid, writeToFile );

% for storing transform offsets
Bin_R = cell(nTransBins,nRotBins);
Bin_Axis = cell(nTransBins,nRotBins);
Bin_Angle = cell(nTransBins,nRotBins);
Bin_T = cell(nTransBins,nRotBins);
% for storing iteration data
Bin_Data = cell(nTransBins,nRotBins,nAlgs);

bin_stats_fields = {...
  'numIterAvg',...
  'timeAvg',...
  'numUnstable',...
  'FREAvg',...
  'ErrAngAvg',...
  'ErrTransAvg'};
% preallocate struct for storing data
empty_cells = repmat(cell({0}),1,numel(bin_stats_fields));
entries = {bin_stats_fields{:}; empty_cells{:}};
Bin_Stats = struct(entries{:});
Bin_Stats(nTransBins,nRotBins,nAlgs).numIterAvg = 0; % initialize struct size
Bin_Stats(:) = Bin_Stats(1);         % initialize struct data to zero


dispCapture(sprintf('Running %u sets of binned trials', nRotBins), fid, writeToFile );
%dispCapture(sprintf('Running %u sets of binned trials', nTransBins*nRotBins), fid, writeToFile );

t = zeros(3,1);
dt = zeros(3,1);

for transBin=1:1 %transBin=1:nTransBins

%   minTrans = TransBins(transBin,1);
%   maxTrans = TransBins(transBin,2);
  
  % Generate random translations for this bin
  %Tact = minTrans*ones(nTrials,3) + (maxTrans-minTrans)*rand(nTrials,3);
  PosNeg = rand(nTrials,3);
  PosInd = find(PosNeg > 0.5);
  NegInd = find(PosNeg <= 0.5);
  PosNeg(PosInd) = 1;
  PosNeg(NegInd) = -1;
  %Tact = Tact.*PosNeg;
  
  Tact = zeros(nTrials,3);
  
  for rotBin=1:nRotBins

    minAng = RotBins(rotBin,1);
    maxAng = RotBins(rotBin,2);

    dispCapture(sprintf('\n\nRunning %u trials for rot range [%u %u]...',...
      nTrials,minAng,maxAng), fid, writeToFile );
%     dispCapture(sprintf('\n\nRunning %u trials for trans/rot range [%u %u]/[%u %u]...',...
%       nTrials,minTrans,maxTrans,minAng,maxAng), fid, writeToFile );
    
    % Generate random rotations for each trial in this bin
    Ract = GenerateRandomRotations( minAng, maxAng, nTrials );
    Bin_R{transBin,rotBin} = Ract;
    actAxis = zeros(nTrials,3);
    actAngle = zeros(nTrials,1);
    for i=1:nTrials
      [actAxis(i,:), actAngle(i)] = rot2AxisAngle(Ract{i});
    end
    Bin_Axis{transBin,rotBin} = actAxis;
    Bin_Angle{transBin,rotBin} = actAngle;
    Bin_T{transBin,rotBin} = Tact;

    % Generate random noise covariance matrices for each point
    Cx = cell(nTrials,1);
    Cy = cell(nTrials,1);
    Rx = GenerateRandomRotations(0,180,nTrials);
    Ry = GenerateRandomRotations(0,180,nTrials);
    eigVal = [0.5,0.5,2.0];
    for i=1:nTrials
      Cy{i} = Ry{i}*diag(eigVal)*Ry{i}';            
      if (anisBothSets)
        Cx{i} = Rx{i}*diag(eigVal)*Rx{i}';
      else
        Cx{i} = diag([0.25,0.25,0.25]);
      end            
    end            
    %Nx = GenerateRandomCovarianceMatrices( 0.1, 1.0, nTrials );
    %Ny = GenerateRandomCovarianceMatrices( 0.1, 1.0, nTrials );  

    for trial=1:nTrials

%       % Generate random noise covariance matrix for each point
%       Mx = zeros(3,3,numPts);
%       My = zeros(3,3,numPts);
%       Rx = GenerateRandomRotations(0,180,numPts);
%       Ry = GenerateRandomRotations(0,180,numPts);
%       eigVal = [0.5,1.0,2.0];
%       for i=1:numPts
%         Mx(:,:,i) = Rx{i}*diag(eigVal)*Rx{i}';
%         My(:,:,i) = Ry{i}*diag(eigVal)*Ry{i}';
%         %Weight{i} = Rot{i}*sqrt(1./eigValue)*Rot{i}';
%       end            
%       %Nx = GenerateRandomCovarianceMatrices( 0.1, 1.0, nTrials );
%       %Ny = GenerateRandomCovarianceMatrices( 0.1, 1.0, nTrials );
%       
      % generate random points      
      spread = diag([1 1 1]).*100;
      xact = (-0.5*ones(numPts,3) + rand(numPts,3))*spread;
      yact = xact*Ract{trial}';% + repmat(Tact(trial,:),[numPts,1]);
      % add noise of random covariance
      dx = mvnrnd([0,0,0], Cx{trial}, numPts);
      xobs = xact + dx;
      dy = mvnrnd([0,0,0], Cy{trial}, numPts);
      yobs = yact + dy;            
      
      % define assumed covariance matrices
      % Note: Mx & My are block-diagonal, with all sub-blocks in each
      %       matrix being equal to each other
      Mxi = repmat(Cx{trial},1,1,numPts);
      Myi = repmat(Cy{trial},1,1,numPts);

      % Registration ground truth is inverse transform of applied
      %  misalignment
      Rgt = Ract{trial};
      tgt = Tact(trial,:)';
      % Rgt = Ract{trial}';
      % tgt = -Ract{trial}'*Tact(trial,:)';
      
      % Iso Method
      name = 'Iso';
      tic
      [R] = Compute3dRotation_Iso(xobs,yobs);
      time = toc;
      numIter = 1;
      
      dR = R*Rgt';
      %dt = tgt - t;
      [~, ErrAng] = rot2AxisAngle(dR);
      ErrAng = ErrAng*180/pi;
      ErrTrans = norm(dt);
      xactReg = bsxfun(@plus, xact*R', t');
      dist = sqrt(sum((yact - xactReg).^2,2));
      FRE = sum(dist,1)/numPts;      
      trialData(1,trial).AlgName = name;
      trialData(1,trial).FRE = FRE;      
      trialData(1,trial).numIter = numIter;
      trialData(1,trial).time = time;      
      trialData(1,trial).ErrAng = ErrAng;      
      trialData(1,trial).ErrTrans = ErrTrans;
      
      % Kanatani Method 
      name = 'Kanatani';
      tic      
      %data = Compute3dRotation_Kanatani_ForLoops( xobs,yobs,Mxi,Myi,maxIter,threshAng );
      data = Compute3dRotation_Kanatani_Vectorized( xobs,yobs,Mxi,Myi,maxIter,threshAng );
      time = toc;
      R = data{1};
      %R = data{1}(:,1:3);
      %t = data{1}(:,4);
      numIter = data{2};
      %LoopData = data{3};
      %numLoopIter = data{4};
      
      %EData(trial).LoopData = LoopData;
      %EData(trial).numLoopIter = numLoopIter;      
      
      dR = R*Rgt';
      dt = tgt - t;      
      [~, ErrAng] = rot2AxisAngle(dR);
      ErrAng = ErrAng*180/pi;
      ErrTrans = norm(dt);
      xactReg = bsxfun(@plus, xact*R', t');
      dist = sqrt(sum((yact - xactReg).^2,2));
      FRE = sum(dist,1)/numPts;      
      trialData(2,trial).AlgName = name;
      trialData(2,trial).FRE = FRE;      
      trialData(2,trial).numIter = numIter;
      trialData(2,trial).time = time;      
      trialData(2,trial).ErrAng = ErrAng;      
      trialData(2,trial).ErrTrans = ErrTrans;
      
      % Proposed Method
      name = 'Proposed';
      tic
      data = Compute3dRotation_TotalLeastSquares(xobs,yobs,Mxi(:,:,1),Myi(:,:,1),maxIter,threshAng);
      %data = ComputeRigidBodyXfm_TotalLeastSquares( xobs,yobs,Mxi,Myi,maxIter,threshAng,threshPos,isoInit );
      time = toc;
      R = data{1};
%       R = data{1}(:,1:3);
%       t = data{1}(:,4);
      numIter = data{2};
      % dTheta = data{3};
      % dAxis = data{4};
      % dTrans = data{5};

      dR = R*Rgt';
      dt = tgt - t;      
      [~, ErrAng] = rot2AxisAngle(dR);
      ErrAng = ErrAng*180/pi;
      ErrTrans = norm(dt);
      xactReg = bsxfun(@plus, xact*R', t');
      dist = sqrt(sum((yact - xactReg).^2,2));
      FRE = sum(dist,1)/numPts;      
      trialData(3,trial).AlgName = name;
      trialData(3,trial).FRE = FRE;      
      trialData(3,trial).numIter = numIter;
      trialData(3,trial).time = time;      
      trialData(3,trial).ErrAng = ErrAng;      
      trialData(3,trial).ErrTrans = ErrTrans;
      
    end

    header = 'trial   ';
    algDisplay = cell(4,1);
    space = repmat('  ',[nTrials,1]);
    for alg=1:nAlgs
      data = trialData(alg,:);
      header = [header,sprintf('\t\t\t%-20s',data(1).AlgName)];      
      %  round displayed results to desired precision
      numIter = int2str([data.numIter]');
      time = num2str(round(10000*[data.time]')/10000);
      FRE = num2str(round(100*[data.FRE]')/100);
      ErrAng = num2str(round(100*[data.ErrAng]')/100);
      ErrTrans = num2str(round(100*[data.ErrTrans]')/100);
      %algDisplay{alg} = [ numIter space time space FRE space ErrAng space ErrTrans ];    
      algDisplay{alg} = [ numIter space time space FRE space ErrAng ];    
    end    
    dispCapture(header, fid, writeToFile);
    dispCapture([int2str([1:nTrials]') space space space ...
      algDisplay{1} space space space ...
      algDisplay{2} space space space ...
      algDisplay{3} ], fid, writeToFile );
%       algDisplay{3} space space space ...
%       algDisplay{4} ], fid, writeToFile );
     
    % Display and store summary results
    algUnstableTrials = cell(nAlgs);
    for alg=1:nAlgs
      algData = trialData(alg,:);
      % compute trial stats for this bin
      numIter = [algData.numIter];
      time = [algData.time];
      ErrAng = [algData.ErrAng];
      ErrTrans = [algData.ErrTrans];
      FRE = [algData.FRE];
      % find unstable trials
      if strcmp(algData(1).AlgName, 'Estepar')
        % special case for Estepar et al (due to dual iterative method)
        unstableTrials = [];
        for i=1:nTrials
          unstable = find( [EData(i).LoopData.numIterR] >= maxIterR );
          if ( (size(unstable,2) > 0) || (EData(i).numLoopIter >= maxIter) )
            unstableTrials(end+1) = i;
          end
        end
      else
        unstableTrials = find( [algData.numIter] >= maxIter );
      end
      if (size(unstableTrials,2) > 0)
        % dispCapture(sprintf('\n%u unstable TLS trials detected!', size(unstableTLS,2)), fid, writeToFile)
        % dispCapture(unstableTLS, fid, writeToFile)   
        numIter(unstableTrials) = [];   % remove unstable trials before computing averages
        time(unstableTrials) = [];
        ErrAng(unstableTrials) = [];
        ErrTrans(unstableTrials) = [];
        FRE(unstableTrials) = [];
      end
      %  compute averages for this bin
      numStable = nTrials - size(unstableTrials,2);
      % ---
      numIterAvg = sum(numIter)/numStable;
      timeAvg = sum(time)/numStable;
      ErrAngAvg = sum(ErrAng)/numStable;
      ErrTransAvg = sum(ErrTrans)/numStable;
      FREAvg = sum(FRE)/numStable;
      % ---
      dispCapture(sprintf('\n%15s Averages:        %.1f  %.4f  %.2f  %.2f',...
        algData(1).AlgName, numIterAvg, timeAvg, FREAvg, ErrAngAvg), fid, writeToFile)
%       dispCapture(sprintf('\n%15s Averages:        %.1f  %.4f  %.2f  %.2f %.2f',...
%         algData(1).AlgName, numIterAvg, timeAvg, FREAvg, ErrAngAvg, ErrTransAvg ), fid, writeToFile)
      % ---
      Bin_Stats(transBin,rotBin,alg).numIterAvg = numIterAvg;
      Bin_Stats(transBin,rotBin,alg).timeAvg = timeAvg;
      Bin_Stats(transBin,rotBin,alg).ErrAngAvg = ErrAngAvg;
      Bin_Stats(transBin,rotBin,alg).ErrTransAvg = ErrTransAvg;
      Bin_Stats(transBin,rotBin,alg).numUnstable = size(unstableTrials,2);
      Bin_Stats(transBin,rotBin,alg).FREAvg = FREAvg;
      algUnstableTrials{alg} = unstableTrials;
      
      Bin_Data{transBin,rotBin,alg} = algData;
    end
    % display unstable trial numbers
    for alg = 1:nAlgs
      unstableTrials = algUnstableTrials{alg};
      if (size(unstableTrials,2) > 0)
        dispCapture(sprintf('\n%u unstable %s trials detected!', size(unstableTrials,2), trialData(alg,1).AlgName), fid, writeToFile)
        dispCapture(unstableTrials, fid, writeToFile)
      end
    end

  end
end

% displays also ErrAng & ErrPos
% dispCapture( sprintf('\nResults Summary  (avgNumIter, avgTime, avgFRE, avgAngErr, numUnstable):'), fid, writeToFile)
% for i=1:1 %nTransBins
%   dispCapture(' ', fid, writeToFile)
%   for j=1:nRotBins
%     dispCapture( sprintf('Bin [%u %u]',...
%       RotBins(j,1),RotBins(j,2)), fid, writeToFile)
% %     dispCapture( sprintf('Bin [%u,%u]/[%u %u]',...
% %       TransBins(i,1),TransBins(i,2),RotBins(j,1),RotBins(j,2)), fid, writeToFile)
%     for alg=1:nAlgs
%       dispCapture(sprintf('     %15s:       %.1f  %.4f  %.2f  %.2f  \t\t%u',...
%         trialData(alg,1).AlgName,...
%         Bin_Stats(i,j,alg).numIterAvg, Bin_Stats(i,j,alg).timeAvg,...
%         Bin_Stats(i,j,alg).FREAvg, Bin_Stats(i,j,alg).ErrAngAvg,... %Bin_Stats(i,j,alg).ErrTransAvg,...
%         Bin_Stats(i,j,alg).numUnstable), fid, writeToFile)
%     end
%   end
% end


dispCapture( sprintf('\nResults Summary  (avgNumIter, avgTime, avgFRE, percentUnstable):'), fid, writeToFile)
  dispCapture(' ', fid, writeToFile)
i = 1;
for j=1:nRotBins
  dispCapture( sprintf('Bin [%u %u]',...
    RotBins(j,1),RotBins(j,2)), fid, writeToFile)
  for alg=1:nAlgs
    dispCapture(sprintf('     %15s:       %4.1f %.4f %.3f %2.0f',...
      trialData(alg,1).AlgName,...
      Bin_Stats(i,j,alg).numIterAvg, Bin_Stats(i,j,alg).timeAvg,...
      Bin_Stats(i,j,alg).FREAvg, 100 * Bin_Stats(i,j,alg).numUnstable / nTrials), fid, writeToFile)
  end
end


% display transposed data
dispCapture(' ', fid, writeToFile);
dispCapture(' ', fid, writeToFile);
dispCapture(' ', fid, writeToFile);
dispCapture( sprintf('Results Summary Transposed:  (avgNumIter, avgTime, avgFRE, percentUnstable)'), fid, writeToFile);
dispCapture(' ', fid, writeToFile);
header = '';
for alg=1:nAlgs
  header = [header,sprintf('\t\t%-18s',trialData(alg,1).AlgName)];
end
dispCapture(header, fid, writeToFile);
i = 1;
dispCapture(' ', fid, writeToFile)  % extra white space between translation bins
for j=1:nRotBins
  lineData = '';
  for alg=1:nAlgs
    lineData = [lineData, sprintf('%.1f %.4f %.3f %2.0f    ',...
      Bin_Stats(i,j,alg).numIterAvg, Bin_Stats(i,j,alg).timeAvg,...
      Bin_Stats(i,j,alg).FREAvg, 100 * Bin_Stats(i,j,alg).numUnstable / nTrials)];
  end
  dispCapture(lineData, fid, writeToFile);
end

if (writeToFile)
  fclose(fid);
end

% save workspace
[fileDir,fileName,fileExt] = fileparts(outputFile);
save(['./',fileName]);