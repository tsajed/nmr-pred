% Adds cell arrays element-by-element.
%
% i.kuprov@soton.ac.uk

function A=plus(A,B)

if iscell(A)&&iscell(B)&&all(size(A)==size(B))
    for n=1:numel(A), A{n}=A{n}+B{n}; end
else
    error('cell array sizes must match.');
end

end

% I came into the room, which was half dark, and presently spotted Lord
% Kelvin in the audience and realized that I was in for trouble at the last
% part of my speech dealing with the age of the Earth, where my views
% conflicted with his. To my relief, Kelvin fell fast asleep, but as I
% came to the important point, I saw the old bird sit up, open an eye and
% cock a baleful glance at me! Then a sudden inspiration came, and I said
% "Lord Kelvin had limited the age of the Earth, provided no new source of
% energy was discovered. That prophetic utterance refers to what we are
% now considering tonight, radium." Behold! The old boy beamed upon me.
% 
% Ernest Rutherford

