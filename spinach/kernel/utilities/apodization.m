% Performs free induction decay apodization. Supports 1D, 2D and 3D FIDs with
% the following syntax:
%
%                  fid=apodization(fid,window_type,params)
%
% Arguments:
%
%     fid         - The free induction decay. The function expects a column
%                   vector in the case of 1D FID, a 2D matrix with the time
%                   origin located at the (1,1) corner point in the case of
%                   a 2D FID, and a 3D matrix with the time origin located
%                   at the (1,1,1) corner point in the case of a 3D FID.
%
%     window_type - Type of the window function. The following window func-
%                   tion types are supported: 
%
%                   'none-1d' divides the first point by 2 and does noth-
%                             ing else.
%
%                  'crisp-1d' divides the first point by 2 and multiplies
%                             the FID by a matched cos^8 half-bell.
%
%                    'exp-1d' divides the first point by 2 and multiplies
%                             the FID by a decaying exponential with the
%                             decay rate specified by the user.
%
%               'gaussian-1d' divides the first point by 2 and multiplies
%                             the FID by a decaying Gaussian with the de-
%                             cay rate specified by the user.
%
%                'cosbell-1d' divides the first point by 2 and multiplies
%                             the FID by a matched cosine half-bell.
%
%                   'echo-1d' multiplies the FID by sine square bell, used
%                             for DEER spin echoes.
%
%                 'kaiser-1d' divides the first point by 2 and multiplies
%                             the FID by a Kaiser function with the decay
%                             rate specified by the user.
%
%                'hamming-1d' divides the first point by 2 and multiplies
%                             the FID by a matched Hamming function.
%
%                   'none-2d' divides the corner point by 2 and does noth-
%                             ing else.
%
%                  'crisp-2d' divides the corner point by 2 and multipli-
%                             es the FID by a matched cos^4 half-bell.
%
%                    'exp-2d' divides the corner point by 2 and multipli-
%                             es the FID by a decaying exponential in both
%                             dimensions with the decay rate specified by
%                             the user.
%
%               'gaussian-2d' divides the corner point by 2 and multipli-
%                             the FID by a decaying Gaussian in both dim-
%                             ensions with the decay rate specified by 
%                             the user.
%
%                'cosbell-2d' divides the corner point by 2 and multipli-
%                             es the FID by a matched cosine half-bell in
%                             both dimensions.
%
%              'sqcosbell-2d' divides the first point by 2 and multiplies
%                             the FID by matched squared cosine half-bell.
%
%                   'echo-2d' multiplies one dimension by the square sine
%                             bell, used for DEER echo processing.
%
%                'cosbell-3d' divides the corner point by 2 and multipli-
%                             es the FID by a matched cosine half-bell in
%                             all three dimensions.
%
%          params - decay rate parameters for those window functions that
%                   require such parameters. If a function does not requi-
%                   re a parameter, the value of params is ignored. Defa-
%                   value of params is zero.
%
% i.kuprov@soton.ac.uk

function fid=apodization(fid,window_type,params)

% Set the defaults
if ~exist('params','var'), params=0; end

% Check consistency
grumble(fid,window_type,params);

% Apply the window function
switch window_type
    
    case 'none-1d'
        fid(1)=fid(1)/2;
    
    case 'crisp-1d'
        fid(1)=fid(1)/2;
        fid=fid.*cos(linspace(0,pi/2,numel(fid))').^8;
    
    case 'exp-1d'
        fid(1)=fid(1)/2;
        fid=fid.*exp(-params*linspace(0,1,numel(fid))');
    
    case 'gaussian-1d'
        fid(1)=fid(1)/2;
        fid=fid.*exp(-params*(linspace(0,1,numel(fid))').^2);
        
    case 'cosbell-1d'
        fid(1)=fid(1)/2;
        fid=fid.*cos(linspace(0,pi/2,numel(fid))');
        
    case 'sqcosbell-1d'
        fid(1)=fid(1)/2;
        fid=fid.*cos(linspace(0,pi/2,numel(fid))').^2;
        
    case 'kaiser-1d'
        fid(1)=fid(1)/2;
        fid=fid.*kaiser(numel(fid),params);
        
    case 'hamming-1d'
        fid(1)=fid(1)/2;
        fid=fid.*hamming(numel(fid));

    case 'echo-1d'
        fid=fid.*sin(linspace(0,pi,numel(fid))').^2;

    case 'echo-2d'
        decay_col=sin(linspace(0,pi,size(fid,1))').^2;
        decay_row=ones(size(fid,2),1);
        fid=fid.*kron(decay_col,decay_row');
        
    case 'none-2d'
        fid(1,1)=fid(1,1)/2;
        
    case 'crisp-2d'
        fid(1,1)=fid(1,1)/2;
        decay_col=cos(linspace(0,pi/2,size(fid,1))').^4;
        decay_row=cos(linspace(0,pi/2,size(fid,2))').^4;
        fid=fid.*kron(decay_col,decay_row');
    
    case 'exp-2d'
        fid(1,1)=fid(1,1)/2;
        decay_col=exp(-params*linspace(0,1,size(fid,1))');
        decay_row=exp(-params*linspace(0,1,size(fid,2))');
        fid=fid.*kron(decay_col,decay_row');
    
    case 'gaussian-2d'
        fid(1,1)=fid(1,1)/2;
        decay_col=exp(-params*(linspace(0,1,size(fid,1))').^2);
        decay_row=exp(-params*(linspace(0,1,size(fid,2))').^2);
        fid=fid.*kron(decay_col,decay_row');
        
    case 'cosbell-2d'
        fid(1,1)=fid(1,1)/2;
        decay_col=cos(linspace(0,pi/2,size(fid,1))');
        decay_row=cos(linspace(0,pi/2,size(fid,2))');
        fid=fid.*kron(decay_col,decay_row');
        
    case 'sqcosbell-2d'
        fid(1,1)=fid(1,1)/2;
        decay_col=cos(linspace(0,pi/2,size(fid,1))').^2;
        decay_row=cos(linspace(0,pi/2,size(fid,2))').^2;
        fid=fid.*kron(decay_col,decay_row');
        
    case 'cosbell-3d'
        fid(1,1,1)=fid(1,1,1)/2;
        [f1,f2,f3]=ndgrid(linspace(0,pi/2,size(fid,1)),...
                          linspace(0,pi/2,size(fid,2)),...
                          linspace(0,pi/2,size(fid,3)));
        fid=fid.*cos(f1).*cos(f2).*cos(f3);
        
    case 'sqcosbell-3d'
        fid(1,1,1)=fid(1,1,1)/2;
        [f1,f2,f3]=ndgrid(linspace(0,pi/2,size(fid,1)),...
                          linspace(0,pi/2,size(fid,2)),...
                          linspace(0,pi/2,size(fid,3)));
        fid=fid.*(cos(f1).^2).*(cos(f2).^2).*(cos(f3).^2);
        
    otherwise
        error(['function ' window_type ' not implemented.']);

end

end

% Consistency enforcement
function grumble(fid,window_type,decay_rate)
if (~isnumeric(fid))||(~isnumeric(decay_rate))
    error('fid and decay rate must be numeric.');
end
if (numel(decay_rate)~=1)||(~isreal(decay_rate))||(decay_rate<0)
    error('decay rate must be a positive real number.');
end
if ~ischar(window_type)
    error('window_type must be a character string.');
end
if strcmp(window_type((end-1):end),'1d')&&...
   ((numel(size(fid))~=2)||(size(fid,2)~=1)||(size(fid,1)<2))
    error('1D window functions require FID to be a column vector.');
end
if strcmp(window_type((end-1):end),'2d')&&...
   ((numel(size(fid))~=2)||any(size(fid)<2))
    error('2D window functions require FID to be a 2D matrix.');
end
if strcmp(window_type((end-1):end),'3d')&&...
   ((numel(size(fid))~=3)||any(size(fid)<2))
    error('3D window functions require FID to be a 3D matrix.');
end
end

% I have had my results for a long time, but I do not yet
% know how I am to arrive at them.
% 
% Carl Friedrich Gauss

