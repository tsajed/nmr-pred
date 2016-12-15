% Merges means and standard deviations of the mean for multiple sample sets
% into a total mean and total standard deviation of the mean. Syntax:
%
%           [total_mean,total_stdm]=statmerge(means,stdms,npoints)
%
% Parameters:
%
%      means - a vector of mean values for each sample set
%
%      stdms - a vector of standard deviations of the mean
%              for each sample set
%
%    npoints - a vector specifying the number of samples in
%              each sample set
%
% i.kuprov@soton.ac.uk

function [total_mean,total_stdm]=statmerge(means,stdms,npoints)

% Compute the combined mean
total_mean=sum(means.*npoints)/sum(npoints);

% Convert stdms to stds
stds=stdms.*sqrt(npoints);

% Combine stds
total_std=sqrt(sum((npoints/sum(npoints)).*(stds.^2)));

% Convert std to stdm
total_stdm=total_std/sqrt(sum(npoints));

end

% Redemption, but not repentance.
%
% Michael Krug

