
clc;
clear;
path=mfilename('fullpath');
[pathstr,name] = fileparts(path);
%a=strfind(pathstr,'/dev');
%addpath(genpath(pathstr(1:a)));
%savepath;
a=strfind(pathstr,'/');
result_dir=strcat(pathstr(1:a(end)),'Results_glutamicum/fisetin');
mkdir(result_dir);
