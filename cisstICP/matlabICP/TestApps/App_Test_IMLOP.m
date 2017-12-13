clear, clc

inputDir = '../../test_data';
meshPath = [inputDir,'/','ProximalFemur.ply'];

% dependencies
dep = {
  '../../matlabICP_build/Release'  % mex interface to cisstICP
  '../m'            % matlab component of cisstICP
  '../m/utilities'  % utilities for transformations, etc.
  };
for i = 1:length(dep)
  addpath( dep{i} );
end

% load mesh
mesh = ply_read2( meshPath );

% create initial offset
rng(5);  % seed
thetas = rand(3,1)*5*pi/180;
R = rotz(thetas(3))*roty(thetas(2))*rotx(thetas(1));
t = rand(3,1)*5;
Fi = getFrm3(R,t);

disp('Initial Offset')
disp(Fi)

% use mesh triangle centers as sample points
samplePts = MeshGetTriangleCenters( mesh.vertices, mesh.faces );
sampleNorms = mesh.face_normals;

% offset the sample points from the mesh
samplePts_xfm = applyFrm3( Fi, samplePts );
sampleNorms_xfm = sampleNorms * getRot(Fi)';


%--- Run ICP ---%

% choose the algorithm
alg = algDirICP_IMLOP();

% set algorithm options
sigma2_init = 1.0;        % initial variance for positional noise model
kinit = 0;                % initial concentration for orientation noise model 
                          %  (0 = orientation doesn't matter for first iteration)
wRpos = 0.5;              % relative influence of sample positions on controlling
                          %  the orientation noise model
bDynamicParamEst = true;  % automatically update noise model parameters as algorithm iterates
bTargetAsMesh = true;     % set to true if registering to a mesh rather than point cloud
bEnableDebug = false;     % enable for verbose output

alg.Initialize( ...
  mesh.vertices, mesh.faces, mesh.face_normals, ...
  samplePts_xfm, sampleNorms_xfm, ...
  sigma2_init, kinit, wRpos, bDynamicParamEst, ...
  bTargetAsMesh, bEnableDebug );

% optimization options
opt = objOptICP();
opt.term_dAng = 0.01*pi/180;    % termination condition: min change in rotation
opt.term_dPos = 0.01;           % termination condition: min change in translation
opt.term_holdIter = 2;          % termination condition: # sequential iterations all termination conditions must be met
opt.maxIter = 200;              % termination condition: max iterations

extras.bEnableDebug = false;
extras.bPlotIter = true;
extras.Vplot = mesh.vertices;
extras.Tplot = mesh.faces;
extras.Xplot = samplePts_xfm;   % unregistered position of sample points
extras.XplotGT = samplePts;     % ground truth position of sample points

[ Freg ] = IterateICP( alg, opt, Fi, extras );
disp('Freg')
disp( Freg )

% compute registration error
dF = invFrm3(Fi)*invFrm3(Freg);
[~,AngErr] = rot2AxisAngle(getRot(dF));
PosErr = norm(getPos(dF));

% compare to initial offset
[~,AngErrInit] = rot2AxisAngle(getRot(Fi));
PosErrInit = norm(getPos(Fi));

disp('AngErrInit  PosErrInit')
disp([AngErrInit*180/pi PosErrInit])

disp('AngErrReg      PosErrReg')
disp([AngErr*180/pi PosErr])

% cleanup
alg.delete();   % calling "clear" also does the trick
