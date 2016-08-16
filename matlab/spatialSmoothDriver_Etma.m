clc
clear all
close all

% Add my code folders.
addpath('./functionMfiles/dataManagement',...
    './functionMfiles/analysisFunctions')

[isOctave,matType,matExt] = checkIsOctave();
                    
% Preproccessing parameters
curDataParams = struct();
curDataParams.binSize = 0.25;

disp('---------------------')
disp('      Smoothing      ')
disp('---------------------')

% Dataset
vFolNams = {'EA1','EA2',...
            'EB1','EB2'};
for folnam_idx = 1:length(vFolNams)
curDataParams.folNam = vFolNams{folnam_idx};
disp(' ')
disp(curDataParams.folNam)

mFileNam = matFileNamSelect('Raw',curDataParams);
load([mFileNam matExt],'fExists','emptySpec')
fExists_raw = fExists;
emptySpec_raw = emptySpec;

w_d = 3;
for w_n = 0:(w_d-1)
if w_n == 0
    curDataParams.wiggle = false;
else
    curDataParams.wiggle = true;
end 
disp(' ')
disp(['  wiggle ' num2str(w_n) ' / ' num2str(w_d)])
curDataParams.wiggle_den = w_d;
curDataParams.wiggle_num = w_n;

curDataParams.smoothParam = 0;

mFileNam = matFileNamSelect('Binned',curDataParams);
load([mFileNam matExt],'mcdata','vbincentrs')
mbdata_raw = mcdata > 0;
clear mcdata
vbincentrs_raw = vbincentrs;

for s = [0.15 0.25]
curDataParams.smoothParam = s;

[mbdata,vbincentrs,fExists,emptySpec,nIter,converged] = ...
    spatiallySmooth(mbdata_raw,vbincentrs_raw,fExists_raw,emptySpec_raw,curDataParams.smoothParam);
[mFileNam,regexpVars] = matFileNamSelect('Binned',curDataParams);
save(matType,[mFileNam matExt],'-regexp',regexpVars)
end

end
end












            






