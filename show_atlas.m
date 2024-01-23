% ADD all ABCD CMIG tools directories to MATLAB path
addpath(genpath('/home/d9smith/github/cmig_tools_internal'));

atlasVersion = 'ABCD2_cor10'; % TODO - update to ABCD3_cor10 ?
loadPrerenderedIfNeeded(atlasVersion) % loads atlas specific data
global PRI % saves atlas data as a global variable so does not need to be loaded every time you run showVol

% 4) Visualize using showVol
coords = [-7 -13 -4]; % Opens showVol to a specific location e.g. R NAcc
showVol(PRI.ABCD2.T1, PRI.ABCD2.CO, PRI.ABCD2.FOD, struct('roiatlas','ABCD2'), coords)


% While running, use these keyboard shortcuts to save images
% -- 'v': cycle to the next volume
% -- 'o': cycle orientations
% -- '!': save screenshot of main axis
