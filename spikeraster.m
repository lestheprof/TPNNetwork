function spikeraster(spikes, varargin)
% SPIKERASTER - Spike raster of multiple neurons
%
%    SPIKERASTER(spikes, options), plot spike trains given in variable
%    SPIKES, which is in format [neuron time].  Optional arguments:
% 
%    'Range', RANGE  Plot only neurons specified in the vector
%                    RANGE.  If neurons are specified that have no
%                    spikes a line will still be made for them in
%                    the raster.
%    'EndTime', ET   Plot up until ET
%    'Lines'         Make a line for each spike to sit on
%    'Xlabel'        Set the x-axis label.  The default is 'time'.
%
%    See also LINE, ISMEMBER
%    Original Author:     David Sterratt <David.C.Sterratt@ed.ac.uk>
%    Created:             Fri Apr 23 1999
%    $Revision: 1.4 $ $Date: 2000/07/04 11:13:42 $
%       Minor revision LSS 2024/08/13
%

% Defaults
if isempty(spikes) 
    disp('Spikeraster called with empty spike train') ;
    return
end
    
range = min(spikes(:,1)):max(spikes(:,1)); % Neurons to plot
endtime = [];				% The final time to plot.
                                        % Set to empty as flag for later
linesflag = 0;				% Line for each spike to
                                        % sit on?
% xlabelstr = 'time';			% xlabel
					
% Read in arguments
for i=1:(nargin-1)
  if ischar(varargin{i})
    switch lower(varargin{i})
     case 'range', range = varargin{i+1};
     case 'lines', linesflag = 1;
     case 'endtime', endtime = varargin{i+1};
     case 'xlabel', xlabelstr = varargin{i+1};
     case 'annotate', annotatestr = varargin{i+1};
    end
  end
end

% If endtime hasn't been specified in the arguments, set it to the
% time of the last spike of the neurons we want to look at (that is
% those specified by range).
if isempty(endtime)
  endtime = max(spikes(ismember(spikes(:,1),range),2));
end

% Prepare the axes
h = newplot;

% Full, Black lines
set(h,'LineStyleOrder', '-')   
set(h,'ColorOrder', [0 0 0])   

% Do the plotting one neuron at a time
for n = 1:length(range)
  s = spikes((spikes(:,1)==range(n))&(spikes(:,2)<=endtime),2);
  line([s'; s'], [(n-0.25)*ones(1,size(s,1)); ...
		  (n+0.25)*ones(1,size(s,1))]);
  % Make the plot the right length
  set(h, 'Xlim', [0 endtime])
  set(h, 'Ylim', [0.5 size(range,2)+0.5])
end

% Add lines for the spikes to sit on if required
if linesflag
  xline=get(gca,'XLim');
  for n=1:size(range,2)
    line(xline, [n n]);
  end
end
% add annotation if required
if exist('annotatestr','var')
    text(0.95 * endtime, size(range,2)+1, annotatestr) ;
end
% add x label if required
if exist('xlabelstr', 'var')
    xlabel(xlabelstr) ;
end
