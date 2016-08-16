clc
clear all
close all

% Add my code folders.
addpath('./functionMfiles/dataManagement',...
    './functionMfiles/analysisFunctions')

[isOctave,matType,matExt] = checkIsOctave();
                    
% Preproccessing parameters
curDataParams = struct();
curDataParams.ppmNormTol = true;
curDataParams.normTol = 10;

% Calibrant Details.
calibrantMZ = [1296.685,1570.677,2147.199,2932.588];
calibrantName = {'Angiotensin I','[Glu]-Fibrinopeptide B', ... 
    'Dynorphin A','ACTH fragment (1-24)'};
calibrantRatio = [0.4 0.4 2 2];

disp('---------------------')
disp('    Normalisation    ')
disp('---------------------')

% Dataset
vFolNams = {'EA1','EA2',...
            'EB1','EB2'};
for folnam_idx = 1:length(vFolNams)
curDataParams.folNam = vFolNams{folnam_idx};
disp(' ')
disp(curDataParams.folNam)

mFileNam = matFileNamSelect('Raw',curDataParams);
load([mFileNam matExt],'L','LXY','Vars')

for var_idx = [2 5 6]
var_to_extract = Vars{var_idx};
switch var_to_extract
    case 'intensity'
        curDataParams.dataType = 'Intensity';
    case 'area'
        curDataParams.dataType = 'Area';
    case 'SN'
        curDataParams.dataType = 'SN';
    otherwise
        error('Problemo')
end
disp(['    ' curDataParams.dataType])



[Shat,nCal,meanCalIntensity,CV] = normaliseDataset(L,LXY,Vars,calibrantMZ,curDataParams.normTol,var_to_extract,true,curDataParams.ppmNormTol);

[mFileNam,regexpVars] = matFileNamSelect('Normalisation',curDataParams);
save(matType,[mFileNam matExt],'-regexp',regexpVars)

end
end






