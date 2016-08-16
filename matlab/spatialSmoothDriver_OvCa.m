clc
clear all
close all

% Add my code folders.
addpath('./functionMfiles/dataManagement',...
    './functionMfiles/analysisFunctions',...
    './functionMfiles/plottingFunctions')
% Some of Steve Marrons code, most particularly pcaSM for PCA.
addpath('./functionMfiles/codeSM/General',...
    './functionMfiles/codeSM/BatchAdjust',...
    './functionMfiles/codeSM/Smoothing')
% Some other peoples code I have used.
addpath('./functionMfiles/otherPeoplesCode')        

[isOctave,matType,matExt] = checkIsOctave();
                    
% Preproccessing parameters
curDataParams = struct();
curDataParams.binSize = 0.25;
curDataParams.wiggle = true;

disp('---------------------')
disp('      Smoothing      ')
disp('---------------------')

% Dataset
vFolNams = {'A1', 'A2', 'A3', 'A4', ...
            'B1', 'B2', 'B3', 'B4', ...
            'C1', 'C2', 'C3', 'C4'};
for folnam_idx = 1:length(vFolNams)
curDataParams.folNam = vFolNams{folnam_idx};
disp(' ')
disp(curDataParams.folNam)

    mFileNam = matFileNamSelect('Raw',curDataParams);
    load([mFileNam matExt],'fExists','emptySpec')
    fExists_raw = fExists;
    emptySpec_raw = emptySpec;

    curDataParams.smoothParam = 0;
    mFileNam = matFileNamSelect('Binned',curDataParams);
    load([mFileNam matExt],'mcdata','vbincentrs')
    mbdata_raw = mcdata > 0;
    clear mcdata
    vbincentrs_raw = vbincentrs;

    curDataParams.smoothParam = 0.25;
    [mbdata,vbincentrs,fExists,emptySpec,nIter,converged] = ...
        spatiallySmooth(mbdata_raw, ...
                        vbincentrs_raw, ...
                        fExists_raw, ...
                        emptySpec_raw, ...
                        curDataParams.smoothParam);
    [mFileNam,regexpVars] = matFileNamSelect('Binned',curDataParams);
    save(matType,[mFileNam matExt],'-regexp',regexpVars)
    
end



% Write output for R for the motivating dataset
curDataParams.folNam = 'A3';
curDataParams.smoothParam = 0.25;

mFileNam = matFileNamSelect('Binned',curDataParams);
load([mFileNam matExt],'mbdata','vbincentrs','fExists','emptySpec')

acq = fExists + emptySpec;
acq(acq > 0) = 1:sum(sum(acq));
acq = acq(fExists);

dataTypeNam = dataTypeNamSelect('Binned',curDataParams);
dlmwrite(['./output/' curDataParams.folNam dataTypeNam '_data.csv'],[0 acq' ; vbincentrs' mbdata])










            






