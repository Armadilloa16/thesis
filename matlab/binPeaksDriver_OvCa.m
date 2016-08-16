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
curDataParams.wiggle = false;

disp('---------------------')
disp('       Binning       ')
disp('---------------------')

vFolNams = {'A1', 'A2', 'A3', 'A4', ...
            'B1', 'B2', 'B3', 'B4', ...
            'C1', 'C2', 'C3', 'C4'};
for folnam_idx = 1:length(vFolNams)
curDataParams.folNam = vFolNams{folnam_idx};
disp(' ')
disp(curDataParams.folNam)

    mFileNam = matFileNamSelect('Raw',curDataParams);
    load([mFileNam matExt],'L','LXY')
    sortedPeaks = sortPeaks(L,LXY);
    clear L LXY

    [mcdata,vbincentrs,mdataL] = bin1dFast(sortedPeaks,curDataParams.binSize);
    [mFileNam,regexpVars] = matFileNamSelect('Binned',curDataParams);
    save(matType,[mFileNam matExt],'-regexp',regexpVars)

end


