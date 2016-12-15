% Pads or cuts a string with spaces to make it n characters long. Syntax:
%
%             output_string=pad(input_string,new_length,side)
%
% where input_string and output_string are character strings and new_length
% is a positive integer number. The optional side argument ('left' or 'right')
% controls the side from which the string is padded.
%
% Konstantin Pervushin (p@trosy.com)
% Ilya Kuprov (i.kuprov@soton.ac.uk)

function output_string=pad(input_string,new_length,side)

% Check consistency
grumble(input_string,new_length);

% String padding
if (~exist('side','var'))||strcmp(side,'right')
    output_string=[input_string blanks(new_length)];
    output_string=output_string(1:new_length);
else
    output_string=[blanks(new_length) input_string];
    output_string=output_string((end-new_length):end);
end

end

% Consistency enforcement
function grumble(input_string,new_length)
if (~ischar(input_string))||(~isrow(input_string))
    error('input_string must be a character string.');
end
if (~isnumeric(new_length))||(~isreal(new_length))||...
   (~isscalar(new_length))||(mod(new_length,1)~=0)||(new_length<1)
    error('new_length mult be a positive integer.');
end
end

% By the grace of reality and the nature of life, man – every man – is an
% end in himself, he exists for his own sake, and the achievement of his
% own happiness is his highest moral purpose.
%
% Ayn Rand, "Atlas Shrugged"

