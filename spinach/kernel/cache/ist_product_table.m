% Structure coefficient tables for the associateive envelopes of su(mult)
% algebras. Syntax:
%
%    [product_table_left,product_table_right]=ist_product_table(mult)
%
% The input parameter is the multiplicity of the spin in question. Output
% contains the structure constants in the following conventions:
%
%    T{n}*T{m}=sum_over_k[product_table_left(n,m,k)*T{k}]
%    
%    T{m}*T{n}=sum_over_k[product_table_right(n,m,k)*T{k}]
%
% Disk caching is used and updated automatically.
%
% hannah.hogben@chem.ox.ac.uk
% i.kuprov@soton.ac.uk

function [product_table_left,product_table_right]=ist_product_table(mult)

% Check consistency
grumble(mult);

% Generate cache record name
table_file=['ist_product_table_' num2str(mult) '.mat'];

% Check the cache
if exist(table_file,'file')
    
    % Lift data from the cache if the file is already available
    load(table_file);
    
else
    
    % Get the irreducible spherical tensors
    T=irr_sph_ten(mult);
    
    % Preallocate the arrays
    product_table_left=zeros(mult^2,mult^2,mult^2);
    product_table_right=zeros(mult^2,mult^2,mult^2);
    
    % Get the structure coefficients
    for m=1:mult^2
        for k=1:mult^2
            normalization=sqrt(trace(T{k}*T{k}')*trace(T{m}*T{m}'));
            for n=1:mult^2
                product_table_left(n,m,k)=trace(T{n}*T{m}*T{k}')/normalization;
                product_table_right(n,m,k)=trace(T{m}*T{n}*T{k}')/normalization;
            end
        end
    end

    % Save the table to the cache
    save(table_file,'product_table_left','product_table_right');

end

end

% Consistency enforcement
function grumble(mult)
if (numel(mult)~=1)||(~isnumeric(mult))||(~isreal(mult))||...
   (mult<2)||(mod(mult,1)~=0)
    error('mult must be a real integer greater than 1.');
end
end

% According to Oxford Chemistry folklore, Peter Atkins has once asked the
% following question at an interview for a Lecturer post: "What is it that
% you have done that a technician would not do?". The candidate produced a
% reasonable answer to the effect that a technician would not have the re-
% quired skills. A better answer was suggested by a member of the Inter-
% view Board some time later: "And what is it that you, Atkins, have done
% that a journalist would not do?"

