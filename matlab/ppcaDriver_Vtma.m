clc
clear all
close all

%%
% Add my code folders.
addpath('./functionMfiles/dataManagement')
% Some of Steve Marrons code, most particularly pcaSM for PCA.
addpath('./functionMfiles/codeSM/General',...
    './functionMfiles/codeSM/BatchAdjust',...
    './functionMfiles/codeSM/Smoothing')

%%
[isOctave,matType,matExt] = checkIsOctave();
  
vNcal = [0];
% vNcal = [0 3 4];
v_minNspec_r = [1];
% v_minNspec_r = [1 100];
v_minNspec_bin = [1];
% v_minNspec_bin = [1 100];

%%
% Calibrant Details.
% calibrantMZ = [1296.685,1570.677,2147.199,2932.588];
% calibrantName = {'Angiotensin I','[Glu]-Fibrinopeptide B', ... 
%     'Dynorphin A','ACTH fragment (1-24)'};


%%
% Proccessing parameters
curDataParams = struct();
curDataParams.binSize = 0.25;
curDataParams.wiggle = true;
curDataParams.wiggle_den = 3;
curDataParams.ppmNormTol = true;
curDataParams.normTol = 10;

disp('---------------------')
disp(' PCA on patient data ')
disp('---------------------')

vFolNams = {'Vtma'};
vDatTypes = {'SN','LogIntensity','Binary','Intensity','Area'};

%% Dataset
for folnam_idx = 1:length(vFolNams)
curDataParams.folNam = vFolNams{folnam_idx};
disp(curDataParams.folNam)

%% Data Type
for datType_idx = 1:length(vDatTypes) 
curDataParams.dataType = vDatTypes{datType_idx};
disp(['  ' curDataParams.dataType ' Data'])

%% Wiggled Binning
for w_n = 0:(curDataParams.wiggle_den-1)
curDataParams.wiggle_num = w_n;
disp(['    Wiggle ' num2str(w_n) ' / ' num2str(curDataParams.wiggle_den)])    

%% Annotated Cancer Spectra
for use_cancer_annot = [false true]
curDataParams.use_cancer_annot = use_cancer_annot;
if use_cancer_annot
    disp(['      ' 'Using only annotated cancer spectra'])
else
    disp(['      ' 'Using all spectra'])
end

%% Calibrant QA filter
for nCal_tol = vNcal
curDataParams.nCal_tol = nCal_tol;
disp(['        ' num2str(curDataParams.nCal_tol) ' calibrant QA filter'])

%% Number of spectra per region filter
for nSpec_r = v_minNspec_r
curDataParams.nSpec_r_tol = nSpec_r;
disp(['          ' 'only using regions with at least ' num2str(nSpec_r) ' spectra'])

%% Number of spectra per bin filter
for nSpec_bin = v_minNspec_bin
curDataParams.nSpec_bin_tol = nSpec_bin;
disp(['            ' 'only using bins with at least ' num2str(nSpec_bin) ' spectra'])

switch curDataParams.dataType
    case 'Binary'
        vNorm = [false];
        vIncludeEmptyVals = [true];
        vSmoothParams = [0 0.15 0.25];

    case {'LogIntensity','Intensity','SN','Area'}
        vNorm = [false true];
        vIncludeEmptyVals = [false true];
        vSmoothParams = [0];
    
end

%% Binary Smoothing
for s = vSmoothParams
curDataParams.smoothParam = s;
disp(['              ' num2str(curDataParams.smoothParam) '-smooth'])
    
%% Normalisation
for n = vNorm
curDataParams.normalisation = n;
if ~n
    disp('                No normalisation')
else
    disp('                With normalisation')
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load Patient Summary
if ~strcmp(curDataParams.dataType,'Binary')
    curDataParams.dataType = 'Binary';
    curDataParams.normalisation = false;
    [mFileNam,regexpVars] = matFileNamSelect('pSummary',curDataParams);
    load([mFileNam matExt],'-regexp',regexpVars)
    pbdata = pdata;
    curDataParams.dataType = vDatTypes{datType_idx};
    curDataParams.normalisation = n;
    
    % Number of spectra per bin filter
    pbdata = pbdata(nSpec >= curDataParams.nSpec_bin_tol,:);
    
end
[mFileNam,regexpVars] = matFileNamSelect('pSummary',curDataParams);
load([mFileNam matExt],'-regexp',regexpVars)

% Number of spectra per bin filter
pdata = pdata(nSpec >= curDataParams.nSpec_bin_tol,:);

% Include empty values when taking averages?
for includeEmptyVals = vIncludeEmptyVals
curDataParams.includeEmptyVals = includeEmptyVals;
if ~includeEmptyVals
    disp('                Not including empty values when averaging')
else
    disp('                Including empty values when averaging')
end
    pdata_temp = pdata;
    if ~curDataParams.includeEmptyVals
        pdata_temp(pbdata > 0) = pdata(pbdata > 0)./pbdata(pbdata > 0);
    end

    for suit_p_only = [false true]
    curDataParams.restrict_p_suit = suit_p_only; 
        if curDataParams.restrict_p_suit
            pdata_temp = pdata_temp(:,p_suit);
        end

        % PCA
        paramStruct_pca = struct(   'npc',      rank(cov(pdata_temp')),    ...
                                    'viout',    [1 1 1 0 1]);
        outStruct_pca = pcaSM(pdata_temp,paramStruct_pca); 

        mFileNam = matFileNamSelect('pPCA',curDataParams);
        save(matType,[mFileNam matExt],'-struct','outStruct_pca')
    
    end
    
    
end
end
end     
end
end
end
end
end
end
end














































