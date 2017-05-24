function [] = init_sbd()
%INIT_SBD   Initializes subdirectories and default config settings.
    
    % Add subdirectories to path
    fp = [fileparts(mfilename('fullpath')) '\'];
    addpath(fp);
    for d = {'utils', 'helpers'}
        addpath(genpath([fp d{1}]));
    end
end