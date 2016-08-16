clc
clear all
close all

% Add my code folders.
addpath('./functionMfiles/dataManagement')

[isOctave,matType,matExt] = checkIsOctave();

disp('---------------------')
disp('     Annotation      ')
disp('---------------------')

% Clinical Data
raw_annot = csvread('..\data\Etma_annotation\Etma_clinical_data_numeric.csv',1,0);
p_list = raw_annot(:,1);
p_lnm  = raw_annot(:,2);
p_suit = raw_annot(:,3);
p_lvsi = raw_annot(:,4);
p_grade= raw_annot(:,5);
p_size = raw_annot(:,6);
p_FIGO_old = raw_annot(:,7);
p_FIGO_old_lvls = {'1A' '1B' '1C' '2A' '3C'};
p_FIGO_new = raw_annot(:,8);
p_FIGO_new_lvls = {'1A' '1B' '3C1' '3C2'};
p_T_old = raw_annot(:,9);
p_T_old_lvls = {'1A' '1A' '1C' '2A'};
p_T_new = raw_annot(:,10);
p_T_new_lvls = {'1A' '1B'};
p_myo_inv = raw_annot(:,11);
p_myo_thi = raw_annot(:,12);
p_ser_dis = raw_annot(:,13);

[mFileNam,regexpVars] = matFileNamSelect('Etma Clinical Data');
save(matType,[mFileNam matExt],'-regexp',regexpVars)



% Preproccessing parameters
curDataParams = struct();

% Dataset
vFolNams = {'EA1','EA2',...
            'EB1','EB2'};

for folnam_idx = 1:length(vFolNams);
curDataParams.folNam = vFolNams{folnam_idx};
disp(curDataParams.folNam)

annot_cancer = csvread(['..\data\Etma_annotation\' curDataParams.folNam '_annotation.csv'],1,0);
annot_cancer = annot_cancer(annot_cancer(:,1),2)' == 1;

mFileNam = matFileNamSelect('Raw',curDataParams);
load([mFileNam matExt],'R','LXY')

R = unique(R);

switch curDataParams.folNam

    % TMA 1B
    case {'EA1','EA2','Etma1B1_1kHz','Etma1B2_1kHz'}
        p_number = [ 0  0  16 16 28 28 17 17 34 ...
                     34 18 18 39 39 25 25 36 0  ...
                     36 11 11 37 37 12 12 51 0  ...
                     51 13 13 46 46 26 26 40 40 ...
                     14 14 47 47 27 27 48 48 20 ...
                     20 49 49 21 21 56 56 24 24 ...
                     57 57 22 22 58 58 55 55 ];
        
    % TMA 2B
    case {'EB1','EB2','Etma2B1_1kHz','Etma2B2_1kHz'}
        p_number = [ 0  0  60 60 61 61 19 19 41 ...
                     41 8  8  1  1  35 35 32 0  ...
                     32 4  4  64 64 63 63 2  0  ...
                     2  30 30 33 33 9  9  65 65 ...
                     6  6  10 10 66 66 5  5  29 ...
                     29 52 52 7  7  53 53 54 54 ...
                     59 59 44 44 38 38 0 ];

    otherwise
        error('Manual class annotations have not been added for this dataset.')
        
end


switch curDataParams.folNam

    % UFX141002_EM-TMA1B_001_003_2000hz_PM
    % Previously 'EM_1B_1'
    case 'EA1'
        missing_cores = [ 20 41 52 53 ];
        ambiguous_regions = [];        
        
    % UFX150507_PM001-EMTMA1B_07__L17_2000hz_PM
    case 'EA2'
        p_number = fliplr(p_number);
        missing_cores = [ 9 11 22 23 25 49 55 56 57 60 ];
        ambiguous_regions = [ 10 ];        

    % UFX150409_PM001-EMTMA2B_05_TMA2B_L2_2000hzPM
    case 'EB1'
        missing_cores = [ 61 ] ;
        ambiguous_regions = [];        
   
    % UFX150512_PM001-EMTMA2B_07__L19_2000hz_PM
    case 'EB2'
        missing_cores = [ 17 27 29 30 37 38 54 ] ;
        ambiguous_regions = [];        
   
    % UFX140625_EM-TMA1B_001_002_1000HZ_PM
    % Previously 'EM_1B_1'
    case 'Etma1B1_1kHz'
        missing_cores = [ 3 7 8 9 41 42 52 ] ;
        ambiguous_regions = [];        
   
    % UFX150514_PM001-EMTMA1B_08__L3_1000hz_PM
    case 'Etma1B2_1kHz'
        missing_cores = [ 3 16 41 46 53 ] ;
        ambiguous_regions = [];        
   
    % UFX141202-EMTMA2B_PM001_04_1000hz_PM
    case 'Etma2B1_1kHz'
        missing_cores = [ 17 19 30 37 38  ] ;
        ambiguous_regions = [61 62];
   
    % UFX150505_PM001-EMTMA2B_06_TMA2B_L17_1000hz_PM
    case 'Etma2B2_1kHz'
        missing_cores = [ 17 20 21 27 29 30 31 37 38 61 ] ;
        ambiguous_regions = [];        
        
end

for r = ambiguous_regions
    p_number = [p_number(1:r-1) 0 p_number(r:end)];
end
temp = true(size(p_number));
temp(missing_cores) = false;
p_number = p_number(temp);
clear temp missing_cores

[mFileNam,regexpVars] = matFileNamSelect('Annotation',curDataParams);
save(matType,[mFileNam matExt],'-regexp',regexpVars)

dlmwrite(['./output/' curDataParams.folNam '_r_p_numbers.csv'],p_number')

end



