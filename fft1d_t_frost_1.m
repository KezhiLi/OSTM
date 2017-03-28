% fft1d_t: implementation of the transpose of dense FFT-based sampling operator; 
% also refer to fft1d_f: implementation of the dense
% FFT-based sampling operator (e.g., partial FFT, scrambled FFT etc.); 

% Written by Thong Do, JHU, 2007
% modified by Lu Gan April. 2008

function x = fft1d_t_frost_1(b, N, OMEGA, p)
% Input Parameters:
% b: A signal vector with length of K;
% N: length of reconstructed signal; 
% OMEGA: a vector with K/2 frequency indices choosen uniform at random;
% p: vector that randomizes samples.
        %p is either a random permutation of [1:N] (for scrambling);
        %or a Bernoulli vector with 1, -1 entries (for random sign
        %flipping);
% Output: 
% x: a signal vector with length of N;
    % 1d DFT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA=sign(randn(1,129));
%A=ones(1,513);
B=zeros(1,256);
B(1:129)=AA(1:129);
for count=1:127;
    B(129+count)=AA(129-count);
end

para=B;

n=256;
nn=n/2;
for i=1:n;
     aa(i)=0;
     for j=1:n;
        aa(i)=aa(i)+cos((j-1)*(i-1)*pi/nn)*para(j);
     end
end
A=toeplitz(aa,aa)/256;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    K = length(b);
    fx = zeros(N,1);
    fx(OMEGA) = sqrt(2)*b(1:K/2) + i*sqrt(2)*b(K/2+1:K);

    x = zeros(N,1);
   if max(p)> 1
        % randomizing by permutation vector
        x(p) = sqrt(N)*real(ifft(fx));
   elseif max(p) == 1
        % randomizing by Bernoulli vector
        x = sqrt(N)*real(ifft(fx));
        x = x.*p;
   end
