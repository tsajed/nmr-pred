% Converts a quaternion rotation specification into a corresponding
% direction cosine matrix. Syntax:
%
%                           dcm=quat2anax(q)
%
% where q is a structure with four fields q.u, q.i, q.j, q.k giving 
% the four components of the quaternion.
%
% i.kuprov@soton.ac.uk

function dcm=quat2dcm(q)

% Check consistency
grumble(q);

% Normalize the quaternion
qnorm=norm([q.u q.i q.j q.k],2);
q.u=q.u/qnorm; q.i=q.i/qnorm;
q.j=q.j/qnorm; q.k=q.k/qnorm;

% Preallocate the answer
dcm=zeros(3,3);

% Compute the answer
dcm(1,1)=q.u^2+q.i^2-q.j^2-q.k^2; dcm(1,2)=2*(q.i*q.j+q.u*q.k);     dcm(1,3)=2*(q.i*q.k-q.u*q.j);
dcm(2,1)=2*(q.i*q.j-q.u*q.k);     dcm(2,2)=q.u^2-q.i^2+q.j^2-q.k^2; dcm(2,3)=2*(q.j*q.k+q.u*q.i);
dcm(3,1)=2*(q.i*q.k+q.u*q.j);     dcm(3,2)=2*(q.j*q.k-q.u*q.i);     dcm(3,3)=q.u^2-q.i^2-q.j^2+q.k^2;

end

% Consistency enforcement
function grumble(q)
if ~all(isfield(q,{'i','j','k','u'}))
    error('quaternion data structure must contain u, i, j, and k fields.');
end
if ~all(isreal([q.u q.i q.j q.k]))
    error('quaternion elements must be real.');
end
end

% Ah, there's nothing more exciting than science. You get all the
% fun of sitting still, being quiet, writing down numbers, paying
% attention... science has it all.
%
% Principal Skinner

