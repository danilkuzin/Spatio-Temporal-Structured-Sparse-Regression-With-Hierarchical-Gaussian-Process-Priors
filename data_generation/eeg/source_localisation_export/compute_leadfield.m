function [lf, grid_out] = compute_leadfield(cfg, data)

% FT_DIPOLEFITTING perform grid search and non-linear fit with one or multiple
% dipoles and try to find the location where the dipole model is best able
% to explain the measured EEG or MEG topography.
%
% This function will initially scan the whole brain with a single dipole on
% a regular coarse grid, and subsequently start at the most optimal location
% with a non-linear search. Alternatively you can specify the initial
% location of the dipole(s) and the non-linear search will start from there.
%
% Use as
%   [source] = ft_dipolefitting(cfg, data)
%
% The configuration has the following general fields
%   cfg.numdipoles  = number, default is 1
%   cfg.symmetry    = 'x', 'y' or 'z' symmetry for two dipoles, can be empty (default = [])
%   cfg.channel     = Nx1 cell-array with selection of channels (default = 'all'),
%                     see FT_CHANNELSELECTION for details
%   cfg.gridsearch  = 'yes' or 'no', perform global search for initial
%                     guess for the dipole parameters (default = 'yes')
%   cfg.nonlinear   = 'yes' or 'no', perform nonlinear search for optimal
%                     dipole parameters (default = 'yes')
%
% If you start with a grid search, the complete grid with dipole
% positions and optionally precomputed leadfields should be specified
%   cfg.grid            = structure, see FT_PREPARE_SOURCEMODEL or FT_PREPARE_LEADFIELD
% The positions of the dipoles can be specified as a regular 3-D
% grid that is aligned with the axes of the head coordinate system
%   cfg.grid.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.grid.resolution = number (e.g. 1 cm) for automatic grid generation
%   cfg.grid.inside     = N*1 vector with boolean value whether grid point is inside brain (optional)
%   cfg.grid.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
% If the source model destribes a triangulated cortical sheet, it is described as
%   cfg.grid.pos        = N*3 matrix with the vertex positions of the cortical sheet
%   cfg.grid.tri        = M*3 matrix that describes the triangles connecting the vertices
% Alternatively the position of a few dipoles at locations of interest can be
% specified, for example obtained from an anatomical or functional MRI
%   cfg.grid.pos        = N*3 matrix with position of each source
%
% If you do not start with a grid search, you have to give a starting location
% for the nonlinear search
%   cfg.dip.pos     = initial dipole position, matrix of Ndipoles x 3
%
% The conventional approach is to fit dipoles to event-related averages, which
% within FieldTrip can be obtained from the FT_TIMELOCKANALYSIS or from
% the FT_TIMELOCKGRANDAVERAGE function. This has the additional options
%   cfg.latency     = [begin end] in seconds or 'all' (default = 'all')
%   cfg.model       = 'moving' or 'regional'
% A moving dipole model has a different position (and orientation) for each
% timepoint, or for each component. A regional dipole model has the same
% position for each timepoint or component, and a different orientation.
%
% You can also fit dipoles to the spatial topographies of an independent
% component analysis, obtained from the FT_COMPONENTANALYSIS function.
% This has the additional options
%   cfg.component   = array with numbers (can be empty -> all)
%
% You can also fit dipoles to the spatial topographies that are present
% in the data in the frequency domain, which can be obtained using the
% FT_FREQANALYSIS function. This has the additional options
%   cfg.frequency   = single number (in Hz)
%
% Low level details of the fitting can be specified in the cfg.dipfit structure
%   cfg.dipfit.display  = level of display, can be 'off', 'iter', 'notify' or 'final' (default = 'iter')
%   cfg.dipfit.optimfun = function to use, can be 'fminsearch' or 'fminunc' (default is determined automatic)
%   cfg.dipfit.maxiter  = maximum number of function evaluations allowed (default depends on the optimfun)
%
% Optionally, you can modify the leadfields by reducing the rank, i.e. remove the weakest orientation
%   cfg.reducerank      = 'no', or number (default = 3 for EEG, 2 for MEG)
%
%
% The volume conduction model of the head should be specified as
%   cfg.headmodel     = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%
% The EEG or MEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions, see FT_DATATYPE_SENS
%   cfg.grad          = structure with gradiometer definition, see FT_DATATYPE_SENS
%   cfg.elecfile      = name of file containing the electrode positions, see FT_READ_SENS
%   cfg.gradfile      = name of file containing the gradiometer definition, see FT_READ_SENS
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
%   cfg.outputfile  =  ...
% If you specify one of these (or both) the input data will be read from a *.mat
% file on disk and/or the output data will be written to a *.mat file. These mat
% files should contain only a single variable, corresponding with the
% input/output structure.
%
% See also FT_SOURCEANALYSIS, FT_PREPARE_LEADFIELD, FT_PREPARE_HEADMODEL

% TODO change the output format, more suitable would be something like:
% dip.label
% dip.time
% dip.avg (instead of Vdata)
% dip.dip.pos
% dip.dip.mom
% dip.dip.model, or dip.dip.avg
% dip.dimord

% Undocumented local options:
%   cfg.dipfit.constr   = Source model constraints, depends on cfg.symmetry
% Optionally, you can include a noise covariance structure to sphere the data (is useful when using both
% magnetometers and gradiometers to fit your dipole)
%   cfg.dipfit.noisecov       = noise covariance matrix, see e.g. FT_TIMELOCK_ANALYSIS

% Copyright (C) 2004-2013, Robert Oostenveld
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$

% these are used by the ft_preamble/ft_postamble function and scripts
ft_revision = '$Id$';
ft_nargin   = nargin;
ft_nargout  = nargout;

% do the general setup of the function
ft_defaults
ft_preamble init
ft_preamble debug
ft_preamble loadvar data
ft_preamble provenance data
ft_preamble trackconfig

% the ft_abort variable is set to true or false in ft_preamble_init
if ft_abort
  return
end

% check if the input data is valid for this function
data = ft_checkdata(data, 'datatype', {'comp', 'timelock', 'freq'}, 'feedback', 'yes');

cfg = ft_checkconfig(cfg, 'renamed', {'hdmfile', 'headmodel'});
cfg = ft_checkconfig(cfg, 'renamed', {'vol',     'headmodel'});

% get the defaults
cfg.channel         = ft_getopt(cfg, 'channel', 'all');
cfg.component       = ft_getopt(cfg, 'component');        % for comp input
cfg.frequency       = ft_getopt(cfg, 'frequency');        % for freq input
cfg.latency         = ft_getopt(cfg, 'latency', 'all');   % for timeclock input
cfg.feedback        = ft_getopt(cfg, 'feedback', 'text');
cfg.gridsearch      = ft_getopt(cfg, 'gridsearch', 'yes');
cfg.nonlinear       = ft_getopt(cfg, 'nonlinear', 'yes');
cfg.symmetry        = ft_getopt(cfg, 'symmetry');
cfg.normalize       = ft_getopt(cfg, 'normalize');      % this is better not used in dipole fitting
cfg.normalizeparam  = ft_getopt(cfg, 'normalizeparam'); % this is better not used in dipole fitting
cfg.backproject     = ft_getopt(cfg, 'backproject');    % this is better not used in dipole fitting
cfg.reducerank      = ft_getopt(cfg, 'reducerank', []); % the default for this is handled below
cfg.dipfit          = ft_getopt(cfg, 'dipfit', []);   % the default for this is handled below

% put the low-level options pertaining to the dipole grid in their own field
cfg = ft_checkconfig(cfg, 'renamed', {'tightgrid', 'tight'}); % this is moved to cfg.grid.tight by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'renamed', {'sourceunits', 'unit'}); % this is moved to cfg.grid.unit by the subsequent createsubcfg
cfg = ft_checkconfig(cfg, 'createsubcfg',  {'grid'});

% the default for this depends on the data type
if ~isfield(cfg, 'model'),
  if ~isempty(cfg.component)
    % each component is fitted independently
    cfg.model = 'moving';
  elseif ~isempty(cfg.frequency)
    % fit the data with a dipole at one location
    cfg.model = 'regional';
  elseif ~isempty(cfg.latency)
    % fit the data with a dipole at one location
    cfg.model = 'regional';
  end
end

if ~isfield(cfg, 'numdipoles')
  if isfield(cfg, 'dip')
    cfg.numdipoles = size(cfg.dip(1).pos,1);
  else
    cfg.numdipoles = 1;
  end
end

% set up the symmetry constraints
if ~isempty(cfg.symmetry)
  if cfg.numdipoles~=2
    error('symmetry constraints are only supported for two-dipole models');
  elseif strcmp(cfg.symmetry, 'x')
    % this structure is passed onto the low-level ft_dipole_fit function
    cfg.dipfit.constr.reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    cfg.dipfit.constr.expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    cfg.dipfit.constr.mirror = [1 1 1 -1 1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 -x1 y1 z1]
  elseif strcmp(cfg.symmetry, 'y')
    % this structure is passed onto the low-level ft_dipole_fit function
    cfg.dipfit.constr.reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    cfg.dipfit.constr.expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    cfg.dipfit.constr.mirror = [1 1 1 1 -1 1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 -y1 z1]
  elseif strcmp(cfg.symmetry, 'z')
    % this structure is passed onto the low-level ft_dipole_fit function
    cfg.dipfit.constr.reduce = [1 2 3];         % select the parameters [x1 y1 z1]
    cfg.dipfit.constr.expand = [1 2 3 1 2 3];   % repeat them as [x1 y1 z1 x1 y1 z1]
    cfg.dipfit.constr.mirror = [1 1 1 1 1 -1];  % multiply each of them with 1 or -1, resulting in [x1 y1 z1 x1 y1 -z1]
  else
    error('unrecognized symmetry constraint');
  end
elseif ~isfield(cfg, 'dipfit') || ~isfield(cfg.dipfit, 'constr')
  % no symmetry constraints have been specified
  cfg.dipfit.constr = [];
end

if ft_getopt(cfg.dipfit.constr, 'sequential', false) && strcmp(cfg.model, 'moving')
  error('the moving dipole model does not combine with the sequential constraint')
  % see http://bugzilla.fieldtriptoolbox.org/show_bug.cgi?id=3119
end

if isfield(data, 'topolabel')
  % this looks like a component analysis
  iscomp = 1;
  % transform the data into a representation on which the timelocked dipole fit can perform its trick
  data = comp2timelock(cfg, data);
else
  iscomp = 0;
end

if isfield(data, 'freq')
  % this looks like a frequency analysis
  isfreq = 1;
  % transform the data into a representation on which the timelocked dipole fit can perform its trick
  data = freq2timelock(cfg, data);
else
  isfreq = 0;
end

% prepare the volume conduction model and the sensor array
% this updates the configuration with the appropriate fields
[headmodel, sens, cfg] = prepare_headmodel(cfg, data);

% set the default for reducing the rank of the leadfields
if isempty(cfg.reducerank)
  if ft_senstype(sens, 'eeg')
    cfg.reducerank = 'no';    % for EEG
  elseif ft_senstype(sens, 'meg') && ft_voltype(headmodel, 'infinite')
    cfg.reducerank = 'no';    % for MEG with a magnetic dipole, e.g. a HPI coil
  elseif ft_senstype(sens, 'meg')
    cfg.reducerank = 'yes';   % for MEG with a current dipole in a volume conductor
  end
end

% select the desired channels, the order should be the same as in the sensor structure
[selcfg, seldata] = match_str(cfg.channel, data.label);
% take the selected channels
Vdata = data.avg(selcfg, :);

% sphere the date using the noise covariance matrix supplied, if any
% this affects both the gridsearch and the nonlinear optimization
noisecov = ft_getopt(cfg.dipfit, 'noisecov');
if ~isempty(noisecov)
  [u, s] = svd(noisecov);
  tol = max(size(noisecov)) * eps(norm(s, inf));
  s = diag(s);
  r1 = sum(s > tol) + 1;
  s(1:(r1 - 1)) = 1 ./ sqrt(s(1:(r1 - 1)));
  s(r1:end)     = 0;
  sphere = diag(s) * u';
  % apply the sphering to the data
  Vdata = sphere * Vdata;
  % apply the sphering as a pre-multiplication to the sensor definition
  montage = [];
  montage.labelold = cfg.channel;
  montage.labelnew = cfg.channel;
  montage.tra = sphere;
  sens = ft_apply_montage(sens, montage, 'balancename', 'sphering');
end

if iscomp
  % select the desired component topographies
  Vdata = Vdata(:, cfg.component);
elseif isfreq
  % the desired frequencies have already been selected
  Vdata = Vdata(:, :);
else
  % select the desired latencies
  if ischar(cfg.latency) && strcmp(cfg.latency, 'all')
    cfg.latency = data.time([1 end]);
  end
  tbeg = nearest(data.time, cfg.latency(1));
  tend = nearest(data.time, cfg.latency(end));
  cfg.latency = [data.time(tbeg) data.time(tend)];
  Vdata = Vdata(:, tbeg:tend);
end

nchans = size(Vdata,1);
ntime  = size(Vdata,2);
Vmodel = zeros(nchans, ntime);
fprintf('selected %d channels\n', nchans);
fprintf('selected %d topographies\n', ntime);

if nchans<cfg.numdipoles*3
  warning('not enough channels to perform a dipole fit');
end

if ntime<1
  error('no spatial topography selected');
end

% check whether EEG is average referenced
if ft_senstype(sens, 'eeg')
  if any(rv(Vdata, avgref(Vdata))>0.001)
    warning('the EEG data is not average referenced, correcting this');
  end
  Vdata = avgref(Vdata);
end

% set to zeros if no initial dipole was specified
if ~isfield(cfg, 'dip')
  cfg.dip.pos = zeros(cfg.numdipoles, 3);
  cfg.dip.mom = zeros(3*cfg.numdipoles, 1);
end

% set to zeros if no initial dipole position was specified
if ~isfield(cfg.dip, 'pos')
  cfg.dip.pos = zeros(cfg.numdipoles, 3);
end

% set to zeros if no initial dipole moment was specified
if ~isfield(cfg.dip, 'mom')
  cfg.dip.mom = zeros(3*cfg.numdipoles, 1);
end

% check the specified dipole model
if numel(cfg.dip.pos)~=cfg.numdipoles*3 || numel(cfg.dip.mom)~=cfg.numdipoles*3
  error('inconsistent number of dipoles in configuration')
end

  % test whether we have a valid configuration for dipole scanning
  if cfg.numdipoles==1
    % this is ok
  elseif cfg.numdipoles==2 && ~isempty(cfg.dipfit.constr)
    % this is also ok
  elseif isfield(cfg.grid, 'pos') && size(cfg.grid.pos,2)==cfg.numdipoles*3
    % this is also ok
  else
    error('dipole scanning is only possible for a single dipole or a symmetric dipole pair');
  end

  % copy all options that are potentially used in ft_prepare_sourcemodel
  tmpcfg = keepfields(cfg, {'grid' 'mri' 'headshape' 'symmetry' 'smooth' 'threshold' 'spheremesh' 'inwardshift'});
  tmpcfg.headmodel = headmodel;
  if ft_senstype(sens, 'eeg')
    tmpcfg.elec = sens;
  else
    tmpcfg.grad = sens;
  end
  % construct the dipole grid on which the gridsearch will be done
  grid = ft_prepare_sourcemodel(tmpcfg);

  ngrid = size(grid.pos,1);

  switch cfg.model
    case 'regional'
      grid.error = nan(ngrid, 1);
    case 'moving'
      grid.error = nan(ngrid, ntime);
    otherwise
      error('unsupported cfg.model');
  end

  insideindx = find(grid.inside);
  thisindx = insideindx(:);
  lf = ft_compute_leadfield(grid.pos(thisindx,:), sens, headmodel, 'reducerank', cfg.reducerank, 'normalize', cfg.normalize, 'normalizeparam', cfg.normalizeparam, 'backproject', cfg.backproject);
  grid_out = grid.pos(thisindx,:);