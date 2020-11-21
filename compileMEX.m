function compileMEX

%   Locate directories
curdir = pwd;
basedir = which('SPSENSE_optset.m');
basedir = basedir(1:end-16);


if ismac
    % Code to run on Mac plaform
    
    cstr1 = 'mex COPTIMFLAGS="\$COPTIMFLAGS -O3" -largeArrayDims ';
    cstr2 = 'mex COPTIMFLAGS="\$COPTIMFLAGS -O3 -DMAXCORES=20" -largeArrayDims ';
    cstr3 = 'mex COPTIMFLAGS="\$COPTIMFLAGS -O3 -DMAXCORES=20" -lpthread -largeArrayDims ';
    
elseif isunix
    % Code to run on Linux plaform
    
    cstr1 = 'mex COPTIMFLAGS="\$COPTIMFLAGS -fopenmp -O3 -DMAXCORES=20" -lgomp -largeArrayDims ';

    %   Create compiler string for functions that can use more cores
    cstr2 = 'mex COPTIMFLAGS="\$COPTIMFLAGS -fopenmp -O3 -DMAXCORES=20" -lgomp -largeArrayDims ';

    %   Create compiler string for functions that use PThreads
    cstr3 = 'mex COPTIMFLAGS="\$COPTIMFLAGS -O3 -DMAXCORES=20" -lpthread -largeArrayDims ';
elseif ispc
    % Code to run on Windows platform
    error('Not specified for PC yet!')
else
    disp('Platform not supported')
end

%   Compile stuff in transform_functions/math_functions
cd([basedir 'Transform_Functions/Math_Functions']);
if exist('mult_add_MEX.c','file')
    fprintf('Compiling mult_add_MEX.c\n');
    eval([cstr1 'mult_add_MEX.c']);
end
if exist('absgradMEX.c','file')
    fprintf('Compiling absgradMEX.c\n');
    eval([cstr1 'absgradMEX.c']);
end
if exist('l1_normMEX.c','file')
    fprintf('Compiling l1_normMEX.c\n');
    eval([cstr1 'l1_normMEX.c']);
end
if exist('l2_norm_GR_MEX.c','file')
    fprintf('Compiling l2_norm_GR_MEX.c\n');
    eval([cstr1 'l2_norm_GR_MEX.c']);
end
if exist('l2_norm_k_MEX.c','file')
    fprintf('Compiling l2_norm_k_MEX.c\n');
    eval([cstr1 'l2_norm_k_MEX.c']);
end
if exist('convVS.c','file')
    fprintf('Compiling convVS.c\n');
    eval([cstr2 'convVS.c']);
end

%   Compile stuff in transform_functions/TV
cd([basedir 'Transform_Functions/Total_Variation']);
if exist('fFDMEX.c','file')
    fprintf('Compiling fFDMEX.c\n');
    eval([cstr1 'fFDMEX.c']);
end
if exist('iFDMEX.c','file')
    fprintf('Compiling iFDMEX.c\n');
    eval([cstr1 'iFDMEX.c']);
end

%   Compile stuff in transform_functions/wavelet
cd([basedir 'Transform_Functions/Wavelet']);
if exist('fwt1dMEX.c','file')
    fprintf('Compiling fwt1dMEX.c\n');
    eval([cstr3 'fwt1dMEX.c']);
end
if exist('iwt1dMEX.c','file')
    fprintf('Compiling iwt1dMEX.c\n');
    eval([cstr3 'iwt1dMEX.c']);
end

%   Compile stuff in transform_functions/GRAPPA
cd([basedir 'Transform_Functions/Parallel_Imaging']);
if exist('fGRAPPA_MEX.c','file')
    fprintf('Compiling fGRAPPA_MEX.c\n');
    eval([cstr1 'fGRAPPA_MEX.c']);
end
if exist('iGRAPPA_MEX.c','file')
    fprintf('Compiling iGRAPPA_MEX.c\n');
    eval([cstr1 'iGRAPPA_MEX.c']);
end
if exist('fSENSE_MEX.c','file')
    fprintf('Compiling fSENSE_MEX.c\n');
    eval([cstr1 'fSENSE_MEX.c']);
end
if exist('iSENSE_MEX.c','file')
    fprintf('Compiling iSENSE_MEX.c\n');
    eval([cstr1 'iSENSE_MEX.c']);
end

%   Return to original directory
cd(curdir);
fprintf('Done\n');

end
