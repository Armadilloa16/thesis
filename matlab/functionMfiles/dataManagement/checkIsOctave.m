function [isOctave,matType,matExt] = checkIsOctave()

    v = version;
    if v(1) == '3'
        isOctave = true;
        matType = '-binary';
        matExt = '';
    elseif v(1) == '7'
        isOctave = false;
        matType = '-v7.3';
        matExt = '.mat';
    else
        disp('WARNING: Version of MATLAB/Octave not recognised? defaulting to MATLAB')
        disp(' ')
        isOctave = false;
        matType = '-v7.3';
        matExt = '.mat';
    end
    
end