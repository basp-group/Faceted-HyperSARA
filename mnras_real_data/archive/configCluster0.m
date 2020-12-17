configCluster

c = parcluster('cirrus R2019a')

% Assign the project code for the job.  **[REQUIRED]**
c.AdditionalProperties.ProjectCode = 'ec110-guest';

% Specify the walltime (e.g. 5 hours).  **[REQUIRED]**
c.AdditionalProperties.WallTime = '96:00:00';

% Specify e-mail address to receive notifications about your job.
c.AdditionalProperties.EmailAddress = 'a.dabbech@hw.ac.uk';

c.AdditionalProperties.ProcsPerNode = 16
% Request a specific reservation to run your job.  It is better to
% use the queues rather than a reservation.
%%c.AdditionalProperties.Reservation = 'your-reservation';

% Set the job placement (e.g., pack, excl, scatter:excl).
% Usually the default of free is what you want.
c.AdditionalProperties.JobPlacement = 'pack'; %'scatter:excl';

% Request to run in a particular queue.  Usually the default (no
% specific queue requested) will route the job to the correct queue.
%%c.AdditionalProperties.QueueName = 'queue-name';

% If you are using GPUs, request up to 4 GPUs per node (this will
% override a requested queue name and will use the 'gpu' queue).
%%c.AdditionalProperties.GpusPerNode = 4;
c.NumWorkers = 50; % number of nodes
c.NumThreads = 1; % number of cores

c.saveProfile
