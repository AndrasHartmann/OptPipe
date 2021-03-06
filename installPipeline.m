function installPipeline()
%function installPipeline()
% Function for installing the pipeline and setting all the path variables

cwd = cd;

%initializing COBRA
cobradir = [cwd '/external/opencobra'];

%addpath([cobradir '/external/libSBML-5.11.4-matlab/']);

%cd ([cobradir '/external/SBMLToolbox-4.1.0/toolbox']);
%InstallSBMLToolbox

addpath(cobradir);
initCobraToolbox

%external library for xls writing
addpath([cwd '/external/xlwrite']);
cd ([cwd '/external/xlwrite']);

Test_xlWrite

%external RobustKnock library 
addpath([cwd '/external/RobustKnock']);

%external Rankproduct library 
addpath([cwd '/external/Rankproduct']);

%methods & common_functions
addpath(genpath([cwd '/methods']));
addpath([cwd '/common_functions']);

%Checking parallel toolbox
toolboxName ='Parallel Computing Toolbox'
v = ver;
if any(strcmp(toolboxName, {v.Name}))&&license('test', 'distrib_computing_toolbox')
   %For older version
   %{
    if matlabpool('size') ~=0
        matlabpool close
    end    
    matlabpool
    %}
    poolobj = gcp('nocreate');
    delete(poolobj);
    parpool
end

cd(cwd);

disp('OptPipe installed successfully.')
