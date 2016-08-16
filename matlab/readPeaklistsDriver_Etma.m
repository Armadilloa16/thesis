clc
clear all
close all

% Add my code folders.
addpath('./functionMfiles/dataManagement',...
    './functionMfiles/analysisFunctions')

% Some other peoples code I have used.
addpath('./functionMfiles/otherPeoplesCode')        

[isOctave,matType,matExt] = checkIsOctave();
                    
% Preproccessing parameters
curDataParams = struct();

disp('---------------------')
disp('  Reading Peaklists  ')
disp('---------------------')

vFolNams = {'EA1','EA2',...
            'EB1','EB2'};
for folnam_idx = 1:length(vFolNams)
    curDataParams.folNam = vFolNams{folnam_idx};
    disp(' ')
    disp(curDataParams.folNam)

    [L,LXY,XYL,X,Y,fExists,Vars,emptySpec,R] = readData(curDataParams.folNam);
    [mFileNam,regexpVars] = matFileNamSelect('Raw',curDataParams);
    save(matType,[mFileNam matExt],'-regexp',regexpVars)
end

