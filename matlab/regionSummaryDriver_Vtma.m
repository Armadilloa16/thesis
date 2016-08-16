%%
% Initialisation clean-up
clc
clear all
close all

% Add my code folders.
addpath('./functionMfiles/dataManagement')

[isOctave,matType,matExt] = checkIsOctave();
                    
vNcal = [0];
% vNcal = [0 3 4];

%%
% Calibrant Details.
calibrantMZ = [1296.685,1570.677,2147.199,2932.588];
calibrantName = {'Angiotensin I','[Glu]-Fibrinopeptide B', ... 
    'Dynorphin A','ACTH fragment (1-24)'};

%%
% Constant Proccessing parameters
curDataParams = struct();
curDataParams.binSize = 0.25;

disp('---------------------')
disp('  Summarise Regions  ')
disp('---------------------')

%% Dataset
vFolNams = {'VA1','VA2','VA3',...
            'VB1','VB2','VB3'};

for folnam_idx = 1:length(vFolNams)
curDataParams.folNam = vFolNams{folnam_idx};
disp(' ')
disp(curDataParams.folNam)

%% Data Type
vDatTypes = {'Binary','LogIntensity','Intensity','SN','Area'};
for datType_idx = 1:length(vDatTypes) 
curDataParams.dataType = vDatTypes{datType_idx};
disp(['  ' curDataParams.dataType ' Data'])
switch curDataParams.dataType
    case 'Binary'
        curDataParams.normalisation = false;
        %% Binary Smoothing
        for s = [0 0.15 0.25]
        curDataParams.smoothParam = s;
        disp(['    ' num2str(curDataParams.smoothParam) '-smooth'])

        %% Wiggled Binning
        curDataParams.wiggle = true;
        curDataParams.wiggle_den = 3;
        for w_n = 0:(curDataParams.wiggle_den-1)
        curDataParams.wiggle_num = w_n;
        disp(['      Wiggle ' num2str(w_n) ' / ' num2str(curDataParams.wiggle_den)])    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load core-patient map and cancer annotation information.
        [mFileNam,~] = matFileNamSelect('Annotation',curDataParams);
        load([mFileNam matExt],'annot_cancer')
        % Load Normalisation Results
        curDataParams.ppmNormTol = true;
        curDataParams.normTol = 10;
        curDataParams.dataType = 'Intensity';
        mFileNam = matFileNamSelect('Normalisation',curDataParams);
        load([mFileNam matExt],'nCal')
        curDataParams.dataType = vDatTypes{datType_idx};
        % Load Region Annotation
        mFileNam = matFileNamSelect('Raw',curDataParams);
        load([mFileNam matExt],'R','emptySpec','fExists')
        temp = (fExists + emptySpec);
        temp2 = zeros(1,sum(sum(temp>0)));
        temp(temp > 0) = 1:length(temp2);
        temp = temp(fExists);
        temp2(temp) = nCal;
        nCal = temp2;
        clear temp2
        if curDataParams.smoothParam > 0
            mFileNam = matFileNamSelect('Binned',curDataParams);
            load([mFileNam matExt],'fExists')
        end
        temp = (fExists + emptySpec)>0;
        non_empty_spec = fExists(temp);                %%%%% non_empty_spec
        r_number = R(non_empty_spec)';                 %%%%% r_number
        r_list = unique(r_number);                     %%%%% r_list
        nR = length(r_list);                           %%%%% nR
        nCal = nCal(non_empty_spec);                   %%%%% nCal
        annot_cancer = annot_cancer(non_empty_spec);   %%%%% annot_cancer
        mFileNam = matFileNamSelect('Binned',curDataParams);
        if curDataParams.smoothParam == 0
            load([mFileNam matExt],'mcdata')
            mbdata = mcdata > 0;
            clear mcdata
        else
            load([mFileNam matExt],'mbdata')
        end                                            %%%%% mbdata
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Annotated Cancer Spectra
        for use_cancer_annot = [false true]
        data_subset = true(1,size(mbdata,2));          %%%%% data_subset
        curDataParams.use_cancer_annot = use_cancer_annot;
        if use_cancer_annot
            data_subset = data_subset & annot_cancer;
            disp(['        ' 'Using only annotated cancer spectra'])
        else
            disp(['        ' 'Using all spectra'])
        end

        %% Calibrant QA filter
        for nCal_tol = vNcal
        curDataParams.nCal_tol = nCal_tol;
        disp(['          ' num2str(curDataParams.nCal_tol) ' calibrant QA filter'])
        data_subset = data_subset & (nCal >= nCal_tol);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Count how many spectra are in each region (after filtering for number of calibrants).
        r_list_count = r_list;
        for r = r_list
            r_list_count(r_list==r) = sum(r_number(data_subset) == r);
        end
        temp = r_list;
        r_list = r_list(r_list_count > 0);
        r_list_count = r_list_count(r_list_count > 0);
        % Make rdata summarising each region.
        rdata = zeros(size(mbdata,1),length(r_list));
        for r = r_list
            rdata(:,r_list==r) = mean(mbdata(:,(data_subset & (r_number == r))),2);
        end
        % Calculate number of peaks per bin
        nSpec = sum(mbdata(:,data_subset),2);
        % Save Results
        [mFileNam,regexpVars] = matFileNamSelect('rSummary',curDataParams);
        save(matType,[mFileNam matExt],'-regexp',regexpVars)
        r_list = temp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        end
        end
        end
        
        
        
        
        
    case {'LogIntensity','Intensity','SN','Area'}
        curDataParams.smoothParam = 0;
        %% Normalisation
        for n = [false true]
        curDataParams.normalisation = n;
        if ~n
            disp('    No normalisation')
        else
            disp('    With normalisation')
        end

        %% Wiggled Binning
        curDataParams.wiggle = true;
        curDataParams.wiggle_den = 3;
        for w_n = 0:(curDataParams.wiggle_den-1)
        curDataParams.wiggle_num = w_n;
        disp(['      Wiggle ' num2str(w_n) ' / ' num2str(curDataParams.wiggle_den)])    

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load core-patient map and cancer annotation information.
        [mFileNam,~] = matFileNamSelect('Annotation',curDataParams);
        load([mFileNam matExt],'annot_cancer')
        % Load Normalisation Results
        curDataParams.ppmNormTol = true;
        curDataParams.normTol = 10;
        mFileNam = matFileNamSelect('Normalisation',curDataParams);
        load([mFileNam matExt],'Shat','nCal','meanCalIntensity') 
                                                       %%%%% meanCalIntensity
                                                       %%%%% nCal
        % Load Region Annotation                       %%%%% Shat
        mFileNam = matFileNamSelect('Raw',curDataParams);
        load([mFileNam matExt],'R','emptySpec','fExists')
        temp = (fExists + emptySpec)>0;
        non_empty_spec = fExists(temp);                %%%%% non_empty_spec
        r_number = R(non_empty_spec)';                 %%%%% r_number
        r_list = unique(r_number);                     %%%%% r_list
        nR = length(r_list);                           %%%%% nR
        annot_cancer = annot_cancer(non_empty_spec);   %%%%% annot_cancer
        % Load Binned Intensity Data
        mFileNam = matFileNamSelect('Binned',curDataParams);
        load([mFileNam matExt],'mdata')                %%%%% mbdata
        % Normalisation.
        if curDataParams.normalisation
            mdata = mdata./repmat(exp(Shat),size(mdata,1),1);
            mdata = mdata./repmat(exp(meanCalIntensity),size(mdata,1),size(mdata,2));
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Annotated Cancer Spectra
        for use_cancer_annot = [false true]
        data_subset = true(1,size(mdata,2));          %%%%% data_subset
        curDataParams.use_cancer_annot = use_cancer_annot;
        if use_cancer_annot
            data_subset = data_subset & annot_cancer;
            disp(['        ' 'Using only annotated cancer spectra'])
        else
            disp(['        ' 'Using all spectra'])
        end

        %% Calibrant QA filter
        for nCal_tol = vNcal
        curDataParams.nCal_tol = nCal_tol;
        disp(['          ' num2str(curDataParams.nCal_tol) ' calibrant QA filter'])
        data_subset = data_subset & (nCal >= nCal_tol);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Count how many spectra are in each region (after filtering for number of calibrants).
        r_list_count = r_list;
        for r = r_list
            r_list_count(r_list==r) = sum(r_number(data_subset) == r);
        end
        temp = r_list;
        r_list = r_list(r_list_count > 0);
        % Make rdata summarising each region.
        rdata = zeros(size(mdata,1),length(r_list));
        switch curDataParams.dataType
            case 'LogIntensity'
                for r = r_list
                    rdata(:,r_list==r) = mean(log(mdata(:,(data_subset & (r_number == r)))+1),2);
                end
                
            case {'Intensity','SN','Area'}
                for r = r_list
                    rdata(:,r_list==r) = mean(mdata(:,(data_subset & (r_number == r))),2);
                end
        end
        % Save Results
        [mFileNam,regexpVars] = matFileNamSelect('rSummary',curDataParams);
        save(matType,[mFileNam matExt],'-regexp',regexpVars)
        r_list = temp;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        end
        end
        end
                
end
end
end
