clc
clear all
close all

%%
% Add my code folders.
addpath('./functionMfiles/dataManagement')

%%
[isOctave,matType,matExt] = checkIsOctave();
  
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

vFolNams = {'Etma'};        

vDatTypes = {'Binary','Intensity','LogIntensity','SN','Area'};
vVRmethods = {'None' 'PCA' 'CCA'};
% vVRmethods = {'None' 'PCA' 'CCA' 'altCCA'};

vUseAnnot = [false true];

vNcal = 0;
% vNcal = [0 3 4];
v_minNspec_r = 1;
% v_minNspec_r = [1 100];
v_minNspec_bin = 1;
% v_minNspec_bin = [1 100];

fid = fopen('.\output\Etma_detailed_LOOclassification_results.txt','wt');

p_cols = [1   2   4   5   6   7   8   9  10  11  12  13  14  17  18  20  21  22  29  30  32  33  35  36  37  38  40  44  46  47  48  49  51  52  53  54  55  56  57  58  60  61  65];

% Print Header
% disp(['VRmethod, DAmethod, DataType, Normalisation, includeEmptyValues, ' ...
%     'Smooth, CancerAnnotation, minNcal, minNspecPerCore, minNspecPerBin, ' ...
%     'nComponents, Wiggle'])
fprintf(fid,['VRmethod,DAmethod,DataType,Normalisation,includeEmptyValues,' ...
    'Smooth,CancerAnnotation,minNcal,minNspecPerCore,minNspecPerBin,' ...
    'nComponents,Wiggle' num2str(p_cols(1:8),',p%i') ...
    num2str(p_cols(9:end),',p%i') '\n']);

%% Dataset
for folnam_idx = 1:length(vFolNams)
curDataParams.folNam = vFolNams{folnam_idx};

%% Data Type
for datType_idx = 1:length(vDatTypes) 
curDataParams.dataType = vDatTypes{datType_idx};
disp([curDataParams.dataType ' Data'])

%% Annotated Cancer Spectra
for use_cancer_annot = [false true]
curDataParams.use_cancer_annot = use_cancer_annot;
if curDataParams.use_cancer_annot
    disp(['  ' 'Annotated Spectra Only'])
else
    disp(['  ' 'All Spectra'])
end    

%% Calibrant QA filter
for nCal_tol = vNcal
curDataParams.nCal_tol = nCal_tol;

%% Number of spectra per region filter
for nSpec_r = v_minNspec_r
curDataParams.nSpec_r_tol = nSpec_r;

%% Number of spectra per bin filter
for nSpec_bin = v_minNspec_bin
curDataParams.nSpec_bin_tol = nSpec_bin;

switch curDataParams.dataType
    case 'Binary'
        vSmoothParams = [0 0.15 0.25];
        vNorm = [false];
        vUseEmptyVals = [true];
        
    case {'LogIntensity','Intensity','SN','Area'}
        vSmoothParams = [0];
        vNorm = [false true];
        vUseEmptyVals = [true false];
        
end


%% Binary Smoothing
for s = vSmoothParams
curDataParams.smoothParam = s;

%% Normalisation
for n = vNorm
curDataParams.normalisation = n;

%% Include empty values when taking averages?
for includeEmptyVals = vUseEmptyVals
curDataParams.includeEmptyVals = includeEmptyVals;

%% Variable Reduction Method
for vr_method_idx = 1:length(vVRmethods)
vr_method_cur = vVRmethods{vr_method_idx};
disp(['  ' '  ' vr_method_cur])
    
switch vr_method_cur
    case 'None'
        vDAmethods = {'NB','DWD'};
    case {'PCA','CCA','altCCA'}
        vDAmethods = {'NB','LDA','DWD'};
end

%% Classification Method
for da_method_idx = 1:length(vDAmethods)
da_method_cur = vDAmethods{da_method_idx};
disp(['  ' '  ' '  ' da_method_cur])
    
switch vr_method_cur
    case 'None'
        curDataParams.classificationMethod = da_method_cur;
        vNcomponents = -1;
        
    case 'PCA'
        curDataParams.classificationMethod = ['pca' da_method_cur];
        curDataParams.restrict_p_suit = true;
        [mFileNam,~] = matFileNamSelect('pPCA',curDataParams);
        load([mFileNam matExt],'veigval')
        vNcomponents = 1:length(veigval);
        clear veigval
        
    case 'CCA'
        curDataParams.classificationMethod = ['cca2' da_method_cur];
        vNcomponents = 1:45;
        
    case 'altCCA'
        curDataParams.classificationMethod = ['cca1' da_method_cur];
        vNcomponents = 1:45;

end

%% Dimension reduction
for cur_nComponents = vNcomponents
    curDataParams.nComponents = cur_nComponents;

    %% Wiggle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for w_n = 0:2
        curDataParams.wiggle_num = w_n;
        
        if ~strcmp(curDataParams.dataType,'Binary')
            temp_datType = curDataParams.dataType;
            temp_norm = curDataParams.normalisation;
            curDataParams.dataType = 'Binary';
            curDataParams.normalisation = false;

            [mFileNam,~] = matFileNamSelect('pSummary',curDataParams);
            load([mFileNam matExt],'p_list','p_suit','vbincentrs')

            curDataParams.dataType = temp_datType;
            curDataParams.normalisation = temp_norm;
        else
            [mFileNam,~] = matFileNamSelect('pSummary',curDataParams);
            load([mFileNam matExt],'p_list','p_suit','vbincentrs')
        end            
        
        [mFileNam,~] = matFileNamSelect('pClassification',curDataParams);
        load([mFileNam matExt],'class_out_LOO')
%         load([mFileNam matExt],'class_out_LOO','dvec','dvec_LOO')
        
%         stab = dvec'*dvec_LOO;
%         stab_av = mean(stab);
%         stab_me = median(stab);
%         stab_sd = sqrt(var(stab));
%         stab_min = min(stab);
%         stab_max = max(stab);
        
        p_list = p_list(p_suit);
        if length(p_list) == length(p_cols)
            if (sum(p_list ~= p_cols') == 0) 
                p_lnm = num2str(class_out_LOO,',%i');
                p_lnm = p_lnm(p_lnm ~= ' ');
            else
                p_lnm = repmat({'NA'},size(p_cols));
                for p_idx = 1:length(p_list)
                    p = p_list(p_idx);
                    p_lnm{p_cols==p} = num2str(class_out_LOO(p_list' == p));
                end
                p_lnm = char(strcat(',',p_lnm))';
                p_lnm = p_lnm(:)';
                p_lnm = p_lnm(p_lnm ~= ' ');
            end
        else
            p_lnm = repmat({'NA'},size(p_cols));
            for p_idx = 1:length(p_list)
                p = p_list(p_idx);
                p_lnm{p_cols==p} = num2str(class_out_LOO(p_list' == p));
            end
            p_lnm = char(strcat(',',p_lnm))';
            p_lnm = p_lnm(:)';
            p_lnm = p_lnm(p_lnm ~= ' ');
        end

%         if sum(isnan(stab)>0)
%             stab_str = 'NA,NA,NA,NA,NA';
%         else
%             stab_str = [num2str(stab_av) ',' ...
%                         num2str(stab_me) ',' ...
%                         num2str(stab_sd) ',' ...
%                         num2str(stab_min) ',' ...
%                         num2str(stab_max)];
%         end

        fprintf(fid,[vr_method_cur ',' ...
            da_method_cur ',' ...
            curDataParams.dataType ',' ...
            num2str(curDataParams.normalisation) ',' ...
            num2str(curDataParams.includeEmptyVals) ',' ...
            num2str(curDataParams.smoothParam) ',' ...
            num2str(curDataParams.use_cancer_annot) ',' ...
            num2str(curDataParams.nCal_tol) ',' ...
            num2str(curDataParams.nSpec_r_tol) ',' ...
            num2str(curDataParams.nSpec_bin_tol) ',' ...
            num2str(curDataParams.nComponents) ',' ...
            num2str(curDataParams.wiggle_num) ...
            p_lnm '\n']);
        
%         disp([stab_str ',' ...
%             curDataParams.folNam ',' ...
%             vr_method_cur ',' ...
%             da_method_cur ',' ...
%             curDataParams.dataType ',' ...
%             num2str(curDataParams.normalisation) ',' ...
%             num2str(curDataParams.includeEmptyVals) ',' ...
%             num2str(curDataParams.smoothParam) ',' ...
%             num2str(curDataParams.use_cancer_annot) ',' ...
%             num2str(curDataParams.nCal_tol) ',' ...
%             num2str(curDataParams.nSpec_r_tol) ',' ...
%             num2str(curDataParams.nSpec_bin_tol) ',' ...
%             num2str(curDataParams.nComponents) ',' ... 
%             num2str(curDataParams.wiggle_num) '\n']);
        
        
        
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
end
end

fclose(fid);













































