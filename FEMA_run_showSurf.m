function [vols beta_hat_spline] = FEMA_run_showSurf(vertvals, statname, fvals, clim, varargin)
% Function to plot vertexwize results of a FEMA analysis.
%% Inputs:

%
%% Output:

inputs = inputParser;
addParameter(inputs,'colormap',blueblackred());
addParameter(inputs,'legendPosition','south');
addParameter(inputs,'curvcontrast',[0.2 0.2]); % contrast of gyri/sulci
addParameter(inputs,'bgcol',[0 0 0]); % change to [1 1 1] for white background
addParameter(inputs,'ico',5); % change to [1 1 1] for white background
addParameter(inputs,'savepath',[]); % change to [1 1 1] for white background
parse(inputs,varargin{:})

load SurfView_surfs.mat % load surface templates

icnum=inputs.Results.ico+1; % index for ico number (icnum = ico + 1)
icnvert = size(icsurfs{icnum}.vertices,1); % indices of vertices for specified icosahedral order

vertvals_lh = vertvals(1:icnvert); % divide statistics by hemisphere for plotting
vertvals_rh = vertvals(icnvert+[1:icnvert]);

% Create figure
fh = figure('Units', 'centimeters', 'Position', [10 10 16 5], 'Color', inputs.Results.bgcol, 'InvertHardcopy', 'off');

% Define spacing for axes
hvgap = [0.02 0.02];
legendPosition = inputs.Results.legendPosition;
if strcmpi(legendPosition, 'south')
    lrgap = [0.02 0.02];
    btgap = [0.2 0.01];
else
    if strcmpi(legendPosition, 'east')
      lrgap = [0.02 0.17];
      btgap = [0.018 0.018];
    end
end

% Create axes
allH = tight_subplot(1,4, hvgap, btgap, lrgap);
hold(allH(:), 'on');

curvcontrast = inputs.Results.curvcontrast; 
bgcol = inputs.Results.bgcol; 
cm = inputs.Results.colormap;

axes(allH(1)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
axes(allH(2)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
axes(allH(3)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'right',[1 0],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;
axes(allH(4)); SurfView_show_new(surf_lh_pial,surf_rh_pial,vertvals_lh,vertvals_rh,fvals,cm,'left', [0 1],curvvec_lh,curvvec_rh,icsurfs{icnum},[],curvcontrast,bgcol); set(gca,'visible','off'); axis tight;

% Set colorbar
colormap(cm);
cb                    = colorbar('color', 'w');
cb.FontSize           = 8;
cb.Label.Interpreter  = 'none';
cb.Label.String       = statname;
cb.Label.FontSize     = 8;
cb.Label.FontWeight   = 'bold';   
cb.Box                = 'off';

if strcmpi(legendPosition, 'south')
    cb.Location = 'south';
    cb.Position(1)      = allH(1).Position(1);
    cb.Position(2)      = cb.Position(2) - btgap(1);
    cb.Position(3)      = 1- allH(1).Position(1) - hvgap(1);
else
    if strcmpi(legendPosition, 'east')
      cb.Location = 'eastoutside';
      cb.Position(1)      = allH(4).Position(1) + allH(4).Position(3) + 0.16;
      cb.Position(2)      = allH(3).Position(2) + .035;
      % cb.Position(4)      = allH(1).Position(4)*2 + hvgap(1);
      cb.Position(4)      = allH(1).Position(4) - .05;
    end
end
caxis(clim);

% save
savepath = inputs.Results.savepath;
if ~isempty(savepath)
    saveas(fh, savepath);
end