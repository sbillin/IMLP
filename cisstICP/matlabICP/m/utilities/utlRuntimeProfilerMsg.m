function utlRuntimeProfilerMsg( bRunProfiler, fstring, varargin )

if bRunProfiler
  fs = ['RuntimeProfiler: ', fstring, '\n   '];
  fprintf( fs,varargin{:} );
  toc;  % display runtime
  tic;  % restart runtime timer  
end

end