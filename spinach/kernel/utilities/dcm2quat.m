% Converts a direction cosine matrix representation of a rotation into
% the unit quaternion representation. Syntax:
%
%                           q=dcm2quat(dcm)
%
% Output: a structure with four fields q.u, q.i, q.j, q.k giving the
% four components of the quaternion. 
%
% i.kuprov@soton.ac.uk

function q=dcm2quat(dcm)

% Check consistency
grumble(dcm);

% Get the trace
tr_dcm=trace(dcm);

% Run the conversion
if tr_dcm>0
    % Most angles are straightforward
    A=sqrt(tr_dcm+1); q.u=0.5*A;
    q.i=(dcm(2,3)-dcm(3,2))/(2*A);
    q.j=(dcm(3,1)-dcm(1,3))/(2*A);
    q.k=(dcm(1,2)-dcm(2,1))/(2*A);
else
    % Zero trace case needs more care
    d=diag(dcm);
    if (d(2)>d(1))&&(d(2)>d(3))
        A=sqrt(d(2)-d(1)-d(3)+1.0); q.j=0.5*A;
        if A~=0, A=0.5/A; end
        q.u=(dcm(3,1)-dcm(1,3))*A;
        q.i=(dcm(1,2)+dcm(2,1))*A;
        q.k=(dcm(2,3)+dcm(3,2))*A;
    elseif d(3)>d(1)
        A=sqrt(d(3)-d(1)-d(2)+1.0); q.k=0.5*A;
        if A~=0, A = 0.5/A; end
        q.u=(dcm(1,2)-dcm(2,1))*A;
        q.i=(dcm(3,1)+dcm(1,3))*A;
        q.j=(dcm(2,3)+dcm(3,2))*A;
    else
        A=sqrt(d(1)-d(2)-d(3)+1.0); q.i=0.5*A;
        if A~=0, A=0.5/A; end
        q.u=(dcm(2,3)-dcm(3,2))*A;
        q.j=(dcm(1,2)+dcm(2,1))*A;
        q.k=(dcm(3,1)+dcm(1,3))*A;
    end
end

end

% Consistency enforcement
function grumble(dcm)
if (~isnumeric(dcm))||(~isreal(dcm))||(~all(size(dcm)==[3 3]))
    error('DCM must be a real 3x3 matrix.');
end
if norm(dcm'*dcm-eye(3))>1e-6
    warning('DCM is not orthogonal to 1e-6 tolerance - conversion accuracy not guaranteed.');
end
if norm(dcm'*dcm-eye(3))>1e-2
    error('DCM is not orthogonal to 1e-2 tolerance, cannot proceed with conversion.');
end
if abs(det(dcm)-1)>1e-6
    warning('DCM determinant is not unit to 1e-6 tolerance - conversion accuracy not guaranteed.');
end
if abs(det(dcm)-1)>1e-2
    error('DCM determinant is not unit to 1e-2 tolerance, cannot proceed with conversion.');
end
end

% I was playing in a tournament in Germany one year when a man approached
% me. Thinking he just wanted an autograph, I reached for my pen, when the
% man made a startling announcement... "I've solved chess!" I sensibly
% started to back away in case the man was dangerous as well as insane, but
% the man continued: "I'll bet you 50 marks that if you come back to my
% hotel room I can prove it to you." Well, 50 marks was 50 marks, so I
% humored the fellow and accompanied him to his room. Back at the room, we
% sat down at his chess board. "I've worked it all out, white mates in 12
% moves no matter what." I played with black perhaps a bit incautiously,
% but I found to my horror that white's pieces coordinated very strangely,
% and that I was going to be mated on the 12th move! I tried again, and I
% played a completely different opening that couldn't possibly result in
% such a position, but after a series of very queer-looking moves, once
% again I found my king surrounded, with mate to fall on the 12th move. I
% asked the man to wait while I ran downstairs and fetched Emmanuel Lasker,
% who was world champion before me. He was extremely skeptical, but agreed
% to at least come and play. Along the way we snagged Alekhine, who was
% then world champion, and the three of us ran back up to the room.
% 
% Lasker took no chances, but played as cautiously as could be, yet after
% a bizarre, pointless-looking series of maneuvers, found himself hemmed
% in a mating net from which there was no escape. Alekhine tried his hand,
% too, but all to no avail.
% 
% It was awful! Here we were, the finest players in the world, men who had
% devoted our very lives to the game, and it was all over! The tournaments,
% the matches, everything - chess had been solved, white wins.
% 
% We killed him, of course.
% 
% Jose Raul Capablanca, only half joking.

