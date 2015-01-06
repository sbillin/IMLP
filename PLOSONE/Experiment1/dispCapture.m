function dispCapture( dispArg, fid, writeToFile )
% displays string to terminal and saves to a file if indicated

disp(dispArg)
if (writeToFile)
  str = evalc('disp(dispArg)');
  fprintf(fid,'%s',str);
end

end