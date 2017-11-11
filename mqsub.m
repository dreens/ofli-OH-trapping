%% This function submits another function to the queue
% Ideally I'd like this to handle data transfer on its own. Perhaps by
% having matlab write the result of the function to a mat file in a
% temporary working directory, and then retrieving it.
%
% The function called must live in the current directory when mqsub is
% called. The function probably needs all of its non-core matlab
% dependencies in this directory as well.
%
% The function must be referenced by a string containing its name as the 
% first argument to mqsub.
%
% The second argument to mqsub is an arraylist containing the arguments
% that will be passed to the function.
%
% Subsequent arguments must come in a standard matlab 'Name','Value' 
% pairing and control SLURM job submission flags, such as the number of 
% cores to use.
% 
% Warning: do not oversave the function during an mqsub call. A workaround
% is to make sure that a function has arguments that allow you to make any
% necessary adjustments by calling the function with different parameters.
%
% Warning: do not call tic or toc without arguments or use global variables
% while submitting jobs to be parallel processed. 
%
% Example Call to the function:
%  mqsub('simdecel',{},'Name','Brian_Compare','Cores',5,'Queue','batch')
function results = mqsub(varargin)

fname = varargin{1}; 
fargs = varargin{2};

if ~exist(fname,'file')
    error(['Function ' fname ' is not found. Change directories, son.']);
end

job = struct(...
    'Name','test_mqsub',...
    'Cores',8,...
    'Memory',1,...       % memory requested in GB
    'Partition','jila',...   % jila, long, slow, standby
    'Time','20:00:00',... %DD-HH:MM:SS
    'Data','../Cluster' ...
    );

for i=3:2:length(varargin)
    fn = varargin{i};
    if isfield(job,fn)
        job.(fn) = varargin{i+1};
    else
        error(['Field ' fn ' not recognized.'])
    end
end

startdirec = pwd;
t = datestr(now,'mmmm-dd-yyyy_HH-MM-SS');
dirname = [job.Data '/' job.Name '_' t];
mkdir(dirname)
cd(dirname);
datadirec = pwd;
save('fargs.mat','fargs');

fidm = fopen('mqsubwrap.m','w');
fprintf(fidm,[...
'%%%% Temporary Script for Qsubbing a matlab function           \n',...
'%% Dave Reens, 4/28/15                                         \n',...
'cd %s                                                          \n',...
'fprintf(''%%s: Opening Pool.\\n'',datestr(now,''HH-MM-SS''))   \n',...
'parpool(%d)                                                    \n',...
'args = open(''fargs.mat'');                                    \n',...
'args = args.fargs;                                             \n',...
'cd %s ;                                                        \n',...
'fprintf(''%%s: Running Func.\\n'',datestr(now,''HH-MM-SS''))   \n',...
'r = %s(args{:});                                               \n',...
'cd %s ;                                                        \n',...
'save(''results_%s.mat'',''r'');                                \n',...
],datadirec,job.Cores,startdirec,fname,datadirec,job.Name);
fclose(fidm);

fidj = fopen('mqsubjob','w');
fprintf(fidj,[...
'#!/bin/bash                                                    \n',...
'module load matlab                                             \n',...
'matlab -nosplash -nodisplay -nodesktop < mqsubwrap.m           \n',...
]);
fclose(fidj);
disp('text of job script at time of qsub:')
fileread('mqsubjob')
systemcall = sprintf('sbatch -J %s -N 1 -n %d -t %s -p %s --mem=%dG mqsubjob',...
    job.Name, job.Cores, job.Time, job.Partition, job.Memory);
disp('executing the following:')
disp(systemcall)
system(systemcall);
%delete('mqsubjob');
%delete('mqsubwrap')
cd(startdirec);
pause(0.05)
%system('squeue | grep dare');
results = 0;

end
