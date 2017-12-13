function [F, retStruct] = IterateICP( alg, opt, F0, extras, algPayload )
%
% inputs:
%  alg ~  algorithm object
%         (an object derived from class algICP.m)
%  opt ~  registration options object
%         (an object of class objOptICP.m)
%  F0 = [R0,t0]  ~  initial guess for the transformation  [optional]
%                   (a 4x4 homogenous transform)
%  extras ~ everything else                               [optional]
%

%profile on
bEnableRuntimeProfiler = false;
runtime = [];
plotName = [];


% --- Input Handling --- %

if ~exist('algPayload','var')
  algPayload = [];
end

if ~exist('extras','var') || isempty(extras)
  % make extras a struct array
  extras.dummyStruct = 0;
end

if ~isfield(extras,'videoPath') || isempty(extras.videoPath)
  bCreateVideo = false;
else  
  bCreateVideo = true;
  videoFrames = cell(1000,1);
  numVideoFrames = 0;  
  videoPath = extras.videoPath;
end

if ~isfield(extras,'fid') || isempty(extras.fid)
  fid = [];
else
  fid = extras.fid;
end

if ~isfield(extras,'bEnableDebug') || isempty(extras.bEnableDebug)
  bEnableDebug = false;
else
  bEnableDebug = extras.bEnableDebug;
end

if ~isfield(extras,'bDisableCallbacks') || isempty(extras.bDisableCallbacks)
  bDisableCallbacks = false;
else
  bDisableCallbacks = extras.bDisableCallbacks;
end

if ~isfield(extras,'bPlotIter') || isempty(extras.bPlotIter)
  bPlotIter = false;  
else
  bPlotIter = extras.bPlotIter;
end

if isfield(extras,'savePlotFile') && ~isempty(extras.savePlotFile)
  bSavePlot = true;
  [~,plotName,~] = fileparts(extras.savePlotFile);
else
  bSavePlot = false;
end

if ~isfield(extras,'bSavePlotVisible') || isempty(extras.bSavePlotVisible)
  bSavePlotVisible = true;
else
  bSavePlotVisible = extras.bSavePlotVisible;
end

if ~isfield(extras,'AZ') || isempty(extras.AZ)
  AZ = [];
else
  AZ = extras.AZ;
end
if ~isfield(extras,'EL') || isempty(extras.EL)
  EL = [];
else
  EL = extras.EL;    
end

if ~isfield(extras,'bReturnFinalMatches') || isempty(extras.bReturnFinalMatches)
  extras.bReturnFinalMatches = false;
end

Vplot = [];
Tplot = [];
sampsPlot = [];
sampsPlotGT = [];
if bSavePlot || bPlotIter
  if isfield(extras,'Vplot')
    Vplot = extras.Vplot;
  end
  if isfield(extras,'Tplot')
    Tplot = extras.Tplot;
  end
  if isfield(extras,'Xplot')
    sampsPlot = extras.Xplot;
  end
  if isfield(extras,'XplotGT')
    sampsPlotGT = extras.XplotGT;
  end
  subPos1 = [0.05 0.05 0.4 0.9];
  subPos2 = [0.55 0.05 0.4 0.9];    
else

end

% initial guess
if ~exist('F0','var') || isempty(F0)
  F0 = getFrm3(eye(3),zeros(3,1));
end
alg.F0 = F0;
F = F0;

if bPlotIter
  hFig = figure(1);
  % hope this stops the crashes
  %set(gcf, 'Renderer', 'zbuffer');  % only use this for saving figures
  
  %set(hFig,'Position',[0 0 1024 768])
  hPlot1 = subplot(1,2,1,'Position',subPos1);  % for initial offset  
  hPlot2 = subplot(1,2,2,'Position',subPos2);  % for moving plot    
  
  bPlotMatches = false;
  bSaveIter = false;
  bSetView = true;  
  if exist('sampsPlotGT','var') && ~isempty(sampsPlotGT)
    bPlotGroundTruth = true;
  else
    bPlotGroundTruth = false;
  end  
  % initial offset
  plotIteration(F, bPlotMatches,bPlotGroundTruth,bSaveIter, bSetView,hPlot1);
  if ~isempty(plotName)
    title(fixFigureString(plotName));    
  end
  % moving plot
  bPlotGroundTruth = false;
  plotIteration(F, bPlotMatches,bPlotGroundTruth,bSaveIter, bSetView,hPlot2);
  
  VideoFrameCapture();
end


% --- Registration Loop --- %

% Initialize Algorithm
if (bEnableRuntimeProfiler)
  tic;  % start timer
end
utlDebugMsg( bEnableDebug,'initializing ICP parameters\n' );
alg.ICP_InitializeParameters( F0, algPayload );
utlRuntimeProfilerMsg( bEnableRuntimeProfiler, 'Registration Loop Initialize' );

% Run Registration Loop
termCount = 0;
iter = 0;
if (~bEnableRuntimeProfiler)
  tic   % start timer for over-all run-time clock
end
while( iter < opt.maxIter )  
  iter = iter + 1;
  
  utlDebugMsg( bEnableDebug,'Beginning iteration %d\n',iter );
    
  % compute matches
  utlDebugMsg( bEnableDebug,'computing matches\n' );
  matches = alg.ICP_ComputeMatches();  
  utlRuntimeProfilerMsg(bEnableRuntimeProfiler, ' iter %d: compute matches', iter);  
  
  if iter == 1
    % evaluate error function
    utlDebugMsg( bEnableDebug,'evaluate error function\n' );
    fval = alg.ICP_EvaluateErrorFunction();  
    utlRuntimeProfilerMsg(bEnableRuntimeProfiler, ' iter %d: evaluate error function', iter);

    % initial iteration callback
    if ~bDisableCallbacks
      utlDebugMsg( bEnableDebug,'iteration callbacks\n' );
      alg.ICP_IterationCallback( iter-1 );  
      utlRuntimeProfilerMsg(bEnableRuntimeProfiler, ' iter %d: iteration callback', iter-1);
    end
    
    if bPlotIter
      % moving plot
      bPlotMatches = true;
      bPlotGroundTruth = false;
      plotIteration(F, bPlotMatches,bPlotGroundTruth,bSavePlot,bSetView,hPlot2);      
      utlRuntimeProfilerMsg(bEnableRuntimeProfiler, ' iter %d: update plot', iter);
      
      VideoFrameCapture();
    end
  end
  
  % register matches
  utlDebugMsg( bEnableDebug,'registering matches\n' );
  [F,dRrod,dt] = alg.ICP_RegisterMatches();  
  dAng = norm(dRrod);
  dPos = norm(dt);  
  utlRuntimeProfilerMsg(bEnableRuntimeProfiler, ' iter %d: register matches', iter);
  
  if bPlotIter
    bPlotMatches = true;
    bPlotGroundTruth = false;
    bSaveIter = false;    
    bSetView = false;
    plotIteration(F, bPlotMatches,bPlotGroundTruth,bSaveIter,bSetView, hPlot2);
    utlRuntimeProfilerMsg(bEnableRuntimeProfiler, ' iter %d: update plot', iter);
    
    VideoFrameCapture();
  end
  
  % evaluate error function
  utlDebugMsg( bEnableDebug,'evaluate error function\n' );
  fval = alg.ICP_EvaluateErrorFunction();  
  utlRuntimeProfilerMsg(bEnableRuntimeProfiler, ' iter %d: evaluate error function', iter);
  
  % iteration callback
  if ~bDisableCallbacks
    utlDebugMsg( bEnableDebug,'iteration callbacks\n' );
    alg.ICP_IterationCallback( iter, fid );  
    utlRuntimeProfilerMsg(bEnableRuntimeProfiler, ' iter %d: iteration callback', iter);
  end
  
  % termination test
  algTerminate = alg.ICP_Terminate();  
  switch algTerminate
    case 1    % by algorithm request, terminate
      break;
    case -1   % by algorithm request, do not allow termination
      termCount = 0;
    case 0    % general termination test
      if dAng < opt.term_dAng && dPos < opt.term_dPos
        termCount = termCount + 1;
      else
        termCount = 0;
      end
      if termCount >= opt.term_holdIter
        break;  % terminate
      end
    otherwise
      error('unrecognized termination return value from algorithm');
  end
end

if (~bEnableRuntimeProfiler)
  runtime = toc;   % report overall run-time
end

if bSavePlot  
  pause on
  
  if bSavePlotVisible
    hFig = figure('Visible','on');
  else
    hFig = figure('Visible','off');
  end
  % hope this stops the crashes
  set(gcf, 'Renderer', 'zbuffer');
  
  hPlot1 = subplot(1,2,1,'Position',subPos1);  % for initial offset  
  hPlot2 = subplot(1,2,2,'Position',subPos2);  % for moving plot    
  
  bSaveIter = false;
  bSetView = true;
  bPlotMatches = false;
  if exist('sampsPlotGT','var') && ~isempty(sampsPlotGT)
    bPlotGroundTruth = true;
  else
    bPlotGroundTruth = false;
  end  
  % initial offset
  plotIteration(F0, bPlotMatches,bPlotGroundTruth,bSaveIter, bSetView,hPlot1);
  %title(fixFigureString(plotName));
  % registered position
  bPlotGroundTruth = false;
  bPlotMatches = true;
  plotIteration(F, bPlotMatches,bPlotGroundTruth,bSaveIter, bSetView,hPlot2);
    
  suptitle(fixFigureString(plotName));  
  % prevent suptitle from automatically displaying figure
  if ~bSavePlotVisible
    set(hFig,'Visible','off');
  end  
  
  saveas(hFig,extras.savePlotFile,'fig');
  saveas(hFig,extras.savePlotFile,'png');
  close(hFig);
end

% % compute final match residuals
% alg.ComputeMatchDistances();

% write video
if bCreateVideo
  %     'Motion JPEG AVI'  - Compressed AVI file using Motion JPEG codec.
  %                          (default)      
  %     'MPEG-4'           - Compressed MPEG-4 file with H.264 encoding 
  %                          (Windows 7 and Mac OS X 10.7 only)
  %     'Uncompressed AVI' - Uncompressed AVI file with RGB24 video.      
  videoWriter = VideoWriter( videoPath, 'Motion JPEG AVI' );
  %videoWriter = VideoWriter( videoPath, 'Uncompressed AVI' );
  videoWriter.FrameRate = 3;
  open(videoWriter);    
  for i = 1:numVideoFrames
    writeVideo( videoWriter, videoFrames{i} );
  end
  close(videoWriter);
end

% extra returned data
retStruct.runtime = runtime;
retStruct.numIter = iter;

if extras.bReturnFinalMatches
  % compute final matches
  utlDebugMsg( bEnableDebug,'computing final matches\n' );
  retStruct.matches = alg.ICP_ComputeMatches();  
end


  % --- Nested Functions --- %
  
  function plotIteration( F, bPlotMatches, bPlotGroundTruth, bSavePlotIteration, bSetView, hPlot) %, plotPos )
        
    if ~exist('bSetView','var') || isempty(bSetView)
      bSetView = false;
    end
    if bSetView
      if ~isempty(AZ) && ~isempty(EL)
        % initialize point-of-view
        view([AZ,EL]);
      end
    end
    if ~isempty(AZ)
      [AZ,EL] = view(); % save current viewpoint (in case user changed it)
    end
        
    % plot mesh
    subplot(hPlot)
    PlotMesh(Vplot,Tplot);
    if isempty(AZ)
      [AZ,EL] = view(); % initialize viewpoint
    end
    hold on
    % plot ground truth
    if bPlotGroundTruth
      hGT = plot3(hPlot,sampsPlotGT(:,1),sampsPlotGT(:,2),sampsPlotGT(:,3),'ko');
      %hGT = plot3(hPlot,sampsPlotGT(:,1),sampsPlotGT(:,2),sampsPlotGT(:,3),'bo');
    end
    % plot matches
    if bPlotMatches
      hMatches = plot3(hPlot,matches(:,1),matches(:,2),matches(:,3),'bo');
    end    
    % plot samples    
    plotXfmSamps = bsxfun(@plus, sampsPlot*getRot(F)', getPos(F)');
    hSamples = plot3(hPlot,plotXfmSamps(:,1),plotXfmSamps(:,2),plotXfmSamps(:,3),'r.');
    hold off
    axis equal;
    %axis vis3d;    % this prevents plots from filling the entire figure
    view([AZ,EL]);  % restore viewpoint
    %set(hPlot,'Position',plotPos);
    drawnow;
    
    if exist('bSavePlotIteration','var') && ~isempty(bSavePlotIteration) && bSavePlotIteration
      if exist('savePlotFile','var') && ~isempty(savePlotFile)
        saveas(hFig,[extras.savePlotFile,'_iter',num2str(iter)],'fig');
        saveas(hFig,[extras.savePlotFile,'_iter',num2str(iter)],'png');
      end      
    end
  end

  function VideoFrameCapture()    
    if bCreateVideo
      % Customize figure for movie recording
      if numVideoFrames == 0;
        % Setting the Renderer property to zbuffer or Painters works around 
        % limitations of getframe with the OpenGL® renderer on some Windows systems.
        set(gca,'nextplot','replacechildren');
        set(gcf,'Renderer','zbuffer');
        set(gcf, 'Position', [100, 100, 1024, 768]);
      end      

      % capture movie frame
      numVideoFrames = numVideoFrames + 1;
      videoFrames{numVideoFrames} = getframe(hFig);
      %videoFrames{numVideoFrames} = getframe(movFig,movieWinSize);
    end
  end

end
