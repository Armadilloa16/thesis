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

disp('---------------------')
disp('  Reading Peaklists  ')
disp('---------------------')

vFolNams = {'A1', 'A2', 'A3', 'A4', ...
            'B1', 'B2', 'B3', 'B4', ...
            'C1', 'C2', 'C3', 'C4'};
for folnam_idx = 1:length(vFolNams)
    curDataParams.folNam = vFolNams{folnam_idx};
    disp(' ')
    disp(curDataParams.folNam)

    [L,LXY,XYL,X,Y,fExists,Vars,emptySpec,R] = readData(curDataParams.folNam);
    [mFileNam,regexpVars] = matFileNamSelect('Raw',curDataParams);
    save(matType,[mFileNam matExt],'-regexp',regexpVars)
end

