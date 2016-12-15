% Spectral differentiation kernel. To be used for accurate numerical
% differentiation of real signals in the following way:
%
% derivative=real(fft(ifft(signal).*fftdiff(1,length(signal),1)'));
%
% i.kuprov@soton.ac.uk

function kern=fftdiff(order,npoints,dx)

if mod(npoints,2)==1
    
    % Kernel for odd point counts
    kern=ifftshift((2i*pi*((1-npoints)/2:((npoints)/2))/(npoints*dx)).^order);
    
elseif mod(npoints,2)==0
    
    % Kernel for even point counts
    kern=ifftshift((2i*pi*(((-npoints)/2):((npoints-1)/2))/(npoints*dx)).^order);
    
else
    
    % Complain and bomb out
    error('npoints parameter must be an integer.');
    
end

end

% A few people laughed, a few people cried, most people were silent. I re-
% membered the line from the Hindu scripture, the Bhagavad-Gita; Vishnu is
% trying to persuade the Prince that he should do his duty and, to impress
% him, takes on his multi-armed form and says, 'Now I am become Death, the
% destroyer of worlds.' I suppose we all thought that, one way or another.
%
% J. Robert Oppenheimer, about the first atomic detonation

