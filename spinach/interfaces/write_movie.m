% Orbits the camera around a 3D plot and writes a movie.
%
% i.kuprov@outlook.com

function write_movie(file_name)

% Open the video writer object
writerObj=VideoWriter(file_name,'MPEG-4'); 
open(writerObj);

% Orbit the camera
for n=1:360
    
    % Grab the frame
    writeVideo(writerObj,getframe(gcf));
    pause(0.05); camorbit(1,0);
    
end

% Close the writer object
close(writerObj);

end

% Love is [...] friendship inspired by beauty.
%
% Marcus Tullius Cicero

