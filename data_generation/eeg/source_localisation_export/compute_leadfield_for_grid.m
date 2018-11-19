function [ lf, grid ] = compute_leadfield_for_grid( EEG, select, xgrid, ygrid, zgrid )

if ~plugin_askinstall('Fieldtrip-lite', 'ft_sourceanalysis'), return; end;

if ~isfield(EEG, 'chanlocs')
  error('No electrodes present');
end

if ~isfield(EEG, 'icawinv')
  error('No ICA components to fit');
end

if ~isfield(EEG, 'dipfit')
  error('General dipolefit settings not specified');
end

if ~isfield(EEG.dipfit, 'vol') & ~isfield(EEG.dipfit, 'hdmfile')
  error('Dipolefit volume conductor model not specified');
end

dipfitdefs
    xgrid_initial  = eval( xgridstr );
    ygrid_initial  = eval( ygridstr );
    zgrid_initial  = eval( zgridstr );

  
  % perform batch fit with single dipole for all selected channels and components
  % warning off;
  %warning backtrace off;  
  comp = eeglab2fieldtrip(EEG, 'componentanalysis', 'dipfit');
  
  cfg = struct('component', select, 'xgrid', xgrid, 'ygrid', ygrid, 'zgrid', zgrid);
  cfg.model      = 'moving';
cfg.gridsearch = 'yes';
cfg.nonlinear  = 'no';
% add some additional settings from EEGLAB to the configuration
tmpchanlocs    = EEG.chanlocs;
cfg.channel    = { tmpchanlocs(EEG.dipfit.chansel).labels };
if isfield(EEG.dipfit, 'vol')
    cfg.vol        = EEG.dipfit.vol;
elseif isfield(EEG.dipfit, 'hdmfile')
    cfg.hdmfile    = EEG.dipfit.hdmfile;
else
    error('no head model in EEG.dipfit')
end
if isfield(EEG.dipfit, 'elecfile') & ~isempty(EEG.dipfit.elecfile)
    cfg.elecfile = EEG.dipfit.elecfile;
end
if isfield(EEG.dipfit, 'gradfile') & ~isempty(EEG.dipfit.gradfile)
    cfg.gradfile = EEG.dipfit.gradfile;
end

%Do some trick to force fieldtrip to use the multiple sphere model
if strcmpi(EEG.dipfit.coordformat, 'CTF')
   cfg = rmfield(cfg, 'channel');
   comp = rmfield(comp, 'elec');
   cfg.gradfile = EEG.dipfit.chanfile;
end


if ~isfield(cfg, 'component')
  % default is to scan all components
  cfg.component = 1:size(comp.topo,2);
end

% for each component scan the whole brain with dipoles using FIELDTRIPs
% dipolefitting function
[lf, grid] = compute_leadfield(cfg, comp);
X = lf;
save('lead_field_matrix_and_grid.mat', 'X', 'grid', 'EEG');

  %warning backtrace on;
