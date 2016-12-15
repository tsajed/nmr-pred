% A simple 1D bin packing algorithm.
%
% i.kuprov@soton.ac.uk

function bins=binpack(box_sizes,bin_size)

% Number the boxes
box_index=(1:numel(box_sizes))';

% Find boxes that are bigger than bins
big_boxes=(box_sizes>bin_size); 
bins=num2cell(box_index(big_boxes));
box_sizes(big_boxes)=[];
box_index(big_boxes)=[];

% Pack the rest of the boxes
while numel(box_index)>0
    
    % Get enough boxes to fill the bin
    current_boxes=(cumsum(box_sizes)<bin_size);

    % Stuff them into the bin
    bins{end+1}=box_index(current_boxes); %#ok<AGROW>
    box_sizes(current_boxes)=[];
    box_index(current_boxes)=[];
    
end

end

% "Ford!" he said, "there's an infinite number of monkeys
% outside who want to talk to us about this script for 
% Hamlet they've worked out."
%
% Douglas Adams, "The Hitchhiker's Guide to the Galaxy"

