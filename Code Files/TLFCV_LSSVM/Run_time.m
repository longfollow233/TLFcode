function Run_time()
% Run TLFCV_LSSVM.
files = dir(fullfile('.\Processed data\','*.mat'));
dataNames = {files.name};
fid = fopen('RunTimes.txt', 'a');
if ftell(fid) == 0
    fprintf(fid, '%-50s %-25s %-25s %-15s\n', 'Dataset', 'Start Time', 'End Time', 'Run Time');
    fprintf(fid, '%-50s %-25s %-25s %-15s\n', '-------', '----------', '--------', '--------');
end
for d = 1:length(dataNames)
    dataName = dataNames{d};
    disp(dataName);
    startTime = datetime('now');
    tic;
    fprintf('Start Time: %s\n', datestr(startTime, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '%-50s %-25s', dataName, datestr(startTime, 'yyyy-mm-dd HH:MM:SS'));
    x = load(['.\Processed data\', dataName]);
    res = best_PandY(x.Xl, x.Yl, x.Xu, x.Yu); %>>>>>
    save(['./Performance/', dataName], 'res');
    endTime = datetime('now');
    elapsedTime = toc;
    fprintf('End Time: %s\n', datestr(endTime, 'yyyy-mm-dd HH:MM:SS'));
    fprintf(fid, '%-25s %-15s\n', datestr(endTime, 'yyyy-mm-dd HH:MM:SS'), sprintf('%.4f s', elapsedTime));
    runTime = endTime - startTime;
    fprintf('Run Time: %s (%.4f seconds)\n\n', char(runTime), elapsedTime);
end
fclose(fid);

end