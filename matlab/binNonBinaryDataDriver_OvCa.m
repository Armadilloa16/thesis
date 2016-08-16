clc
close all
clear all

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

curDataParams = struct();

curDataParams.binSize = 0.25;
curDataParams.smoothParam = 0;
curDataParams.wiggle = false;

disp('---------------------')
disp('   Non-Binary Data   ')
disp('---------------------')

curDataParams.folNam = 'A3';
disp(' ')
disp(curDataParams.folNam)
    
mFileNam = matFileNamSelect('Raw',curDataParams);
load([mFileNam matExt],'L','Vars')

curDataParams.dataType = 'Binary';
mFileNam = matFileNamSelect('Binned',curDataParams);
load([mFileNam matExt],'mcdata','mdataL')

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

    mdata = zeros(size(mdataL));
    peaks = find(mdataL>0);
    for peak_idx = 1:length(peaks)
        [var,spec] = ind2sub(size(mdataL),peaks(peak_idx));
        mdata(var,spec) = L{spec}(mdataL(var,spec),strcmp(Vars,var_to_extract));
        if mdata(var,spec) == 0
            disp('zero value? replaced with 10^-10')
            mdata(var,spec) = 10^(-10);
        end
    end

    groups = find(mcdata>1);
    groupSize = mcdata(groups);
    for group_idx = 1:length(groups)
        [var,spec] = ind2sub(size(mcdata),groups(group_idx));
        mdata(var,spec) = sum(L{spec}(mdataL(var,spec):(mdataL(var,spec)+groupSize(group_idx)-1),strcmp(Vars,var_to_extract)));
    end

    [mFileNam,regexpVars] = matFileNamSelect('Binned',curDataParams);
    save(matType,[mFileNam matExt],'-regexp',regexpVars)

end






   