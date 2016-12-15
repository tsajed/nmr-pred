% Blue -> white -> red color map. White always corresponds
% to value zero. 
%   
% Cunjie Zhang
% Ilya Kuprov

function newmap=b2r(cmin,cmax)

% Check the input
if nargin~=2
   error('incorrect number of input arguments.')
end
if cmin>=cmax
   error('the first argument must be smaller than the second one.');
end

% Set basic colors
red_top     = [1 0 0];
white_middle= [1 1 1];
blue_bottom = [0 0 1];

% Set up color interpolation 
color_num = 255;   
color_input = [blue_bottom;  white_middle;  red_top];
oldsteps = linspace(-1, 1, length(color_input));
newsteps = linspace(-1, 1, color_num);  

% Interpolate colors
newmap_all = NaN(size(newsteps,2),3);
if (cmin<0)&&(cmax>0)
    
    if abs(cmin)<cmax
        
        for j=1:3
            newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
        end
        start_point = round((cmin+cmax)/2/cmax*color_num);
        newmap = squeeze(newmap_all(start_point:color_num,:));
        
    elseif abs(cmin)>=cmax
        
        for j=1:3
            newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
        end
        end_point = round((cmax-cmin)/2/abs(cmin)*color_num);
        newmap = squeeze(newmap_all(1:end_point,:));
        
    end
    
elseif cmin>=0
    
    for j=1:3
        newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
    end
    start_point = round((cmin+cmax)/2/cmax*color_num);
    newmap = squeeze(newmap_all(start_point:color_num,:));
    
elseif cmax <= 0
    
    for j=1:3
        newmap_all(:,j) = min(max(transpose(interp1(oldsteps, color_input(:,j), newsteps)), 0), 1);
    end
    end_point = round((cmax-cmin)/2/abs(cmin)*color_num);
    newmap = squeeze(newmap_all(1:end_point,:));
    
end

end
    
