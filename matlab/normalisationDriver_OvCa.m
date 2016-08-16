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
                    
% Calibrant Details.
calibrantMZ = [1296.685,1570.677,2147.199,2932.588];
calibrantName = {'Angiotensin I','[Glu]-Fibrinopeptide B', ... 
    'Dynorphin A','ACTH fragment (1-24)'};
calibrantRatio = [0.4 0.4 2 2];

% Preproccessing parameters
curDataParams = struct();
curDataParams.ppmNormTol = true;
curDataParams.normTol = 10;

% Dataset
curDataParams.folNam = 'A3';

mFileNam = matFileNamSelect('Raw',curDataParams);
load([mFileNam matExt],'L','LXY','Vars')

var_to_extract = Vars{5};
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

[Shat,nCal,meanCalIntensity,CV] = normaliseDataset(L,LXY,Vars,calibrantMZ,curDataParams.normTol,var_to_extract,true,curDataParams.ppmNormTol);

[mFileNam,regexpVars] = matFileNamSelect('Normalisation',curDataParams);
save(matType,[mFileNam matExt],'-regexp',regexpVars)







