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

% use mesh vertices as samples
samps = applyFrm3( Fi, mesh.vertices );


%--- Run ICP ---%

% choose the algorithm
alg = algICP_StdICP();
algOpts.bEnableDebug = false;
algOpts.bTargetAsMesh = true;
alg.Initialize( mesh.vertices, mesh.faces, samps, algOpts );

% optimization options
opt = objOptICP();
opt.term_dAng = 0.01*pi/180;
opt.term_dPos = 0.01;
opt.term_holdIter = 2;
opt.maxIter = 200;

extras.bEnableDebug = false;
extras.bPlotIter = true;
extras.Vplot = mesh.vertices;
extras.Tplot = mesh.faces;
extras.Xplot = samps;
extras.XplotGT = mesh.vertices;

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
