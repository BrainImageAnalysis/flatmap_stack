


%%

%some of my tools that are needed
addpath includes


%%

% download, compile and set path to the sta-toolbox:
% https://bitbucket.org/skibbe/sta-toolbox/wiki/Home
addpath PATHTO/sta-toolbox/

sta_setPaths


%%
% path toe the gifti tools. Download them from here:
% https://www.artefact.tk/software/matlab/gifti/
addpath PATHTO/gifti-master/

%%

%STPT template
TC_avg_nii_fn = 'data/avg_TC_std.nii.gz'; 
%%
%Where to store debug data
debug_out_fd = 'PATHTO/test/';

%%
