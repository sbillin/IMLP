function printfCapture( str, fid, bWriteToFile, bPrint )
% displays string to terminal and saves to a file if indicated

if ~exist('bPrint','var')
  bPrint = true;
end

if bPrint
  fprintf('%s',str);
end

if exist('fid','var') && ~isempty(fid)
  if ~exist('bWriteToFile','var') || isempty(bWriteToFile)
    bWriteToFile = true;
  end
  if (bWriteToFile)
    fprintf(fid,'%s',str);
  end
end

end