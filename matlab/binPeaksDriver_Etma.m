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
disp('       Binning       ')
disp('---------------------')

vFolNams = {'EA1','EA2',...
            'EB1','EB2'};
for folnam_idx = 1:length(vFolNams)
curDataParams.folNam = vFolNams{folnam_idx};
disp(' ')
disp(curDataParams.folNam)

mFileNam = matFileNamSelect('Raw',curDataParams);
load([mFileNam matExt],'L','LXY')
sortedPeaks = sortPeaks(L,LXY);
clear L LXY

curDataParams.wiggle = false;
[mcdata,vbincentrs,mdataL] = bin1dFast(sortedPeaks,curDataParams.binSize);
[mFileNam,regexpVars] = matFileNamSelect('Binned',curDataParams);
save(matType,[mFileNam matExt],'-regexp',regexpVars)

curDataParams.wiggle = true;
w_d = 3;
curDataParams.wiggle_den = w_d;
for w_n = 1:(w_d-1)
        
disp(['  wiggle ' num2str(w_n) ' / ' num2str(w_d)])
curDataParams.wiggle_num = w_n;

if curDataParams.wiggle
    if isfield(curDataParams,'wiggle_den')
        if curDataParams.wiggle_den == 2
            [mcdata,vbincentrs,mdataL] = bin1dFast(sortedPeaks,curDataParams.binSize,(curDataParams.binSize/2));
        elseif isfield(curDataParams,'wiggle_num')
            if curDataParams.wiggle_num <= (curDataParams.wiggle_den/2)
                [mcdata,vbincentrs,mdataL] = bin1dFast(sortedPeaks,curDataParams.binSize,curDataParams.wiggle_num*(curDataParams.binSize/curDataParams.wiggle_den));
            else                
                [mcdata,vbincentrs,mdataL] = bin1dFast(sortedPeaks,curDataParams.binSize,(curDataParams.wiggle_num*(curDataParams.binSize/curDataParams.wiggle_den))-curDataParams.binSize);
            end
        end
    else
        [mcdata,vbincentrs,mdataL] = bin1dFast(sortedPeaks,curDataParams.binSize,(curDataParams.binSize/2));
    end
else
    [mcdata,vbincentrs,mdataL] = bin1dFast(sortedPeaks,curDataParams.binSize);
end

[mFileNam,regexpVars] = matFileNamSelect('Binned',curDataParams);
save(matType,[mFileNam matExt],'-regexp',regexpVars)

end

end


