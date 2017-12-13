function utlDebugMsg( bPrint, fstring, varargin )

if bPrint
  fprintf( fstring,varargin{:} );
end

end