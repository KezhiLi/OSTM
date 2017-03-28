%function test1d();

% Kezhi Li  12/03/2009

% parameters
N = 256; % signal length
% total_samp=100;% number of measurements
k = 20;% sparsity

% sparsifying transform
Psi = idct(eye(N));

% gaussian measurement matrix
Mtx1 = orth(randn(N,N));   

% partial Hadamard with random permuatation

Mtx2 = hadamard(N)/sqrt(N);

q2 = [1,randperm(N-1)+1];

p2 = randperm(N);

% %block-diagonal Hadamard with random permutation;
% blk_size = 32;
% Mtx3=kron(eye(N/blk_size),hadamard(blk_size)/sqrt(blk_size));
% q3 = randperm(N);
% p3 = randperm(N);

% %partial Hadamard with random sign reversal;
% Mtx4 = hadamard(N)/sqrt(N);
% p4 = 2*round(rand(N,1))-1;
% Mtx4 = Mtx4*diag(p4);
% q4 = randperm(N);
% AAa=sign(randn(1,blk_size/2+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
AA=sign(randn(1,N/2+1));

B=zeros(1,N);
B(1:N/2+1)=AA(1:N/2+1);
for count=1:N/2-1;
    B(N/2+1+count)=AA(N/2+1-count);
end

para=B;

n=N;
nn=n/2;
for i=1:n;
     aa(i)=0;
     for j=1:n;
        aa(i)=aa(i)+cos((j-1)*(i-1)*pi/nn)*para(j);
     end
end
A=toeplitz(aa,aa)/N;

Mtx5 = A;
q5 = [1,randperm(N-1)+1];

p5 = randperm(N);
%p5=1:N;
% Teoplitz random matrix

%temp1 = randn(1,2*N-1);
temp1=sign(randn(1,2*N-1))/sqrt(N);
Mtx6 = toeplitz(temp1(1:N),[temp1(1),temp1(N+1:2*N-1)]);

q6 = [1,randperm(N-1)+1];

p6 = randperm(N);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

trial_num = 500;
success_Gauss = zeros(13,1);
success_Had_RP = zeros(13,1);
%success_blkHad = zeros(13,1);
%success_blkOSTM = zeros(13,1);
success_OSTM = zeros(13,1);
success_Toep = zeros(13,1);


for j =1:13
    samp_num = 60+(j-1)*5;
    % Gaussian Sampling Operator;
    A1 = Mtx1(1:samp_num,:);
    % Full WHT Sampling Operator with random permutation;
    A2 = Mtx2(q2(1:samp_num),p2);
%      % Block WHT Sampling Operator with random permutation;
%     A3 = Mtx3(q3(1:samp_num),p3);
%     % OSTM32;
%     A4 = Mtx4(q4(1:samp_num),p4);
    % Orthogonal Symmetric Teoplitz Matrix;
    A5 = Mtx5(q5(1:samp_num),p5);
    % Toeplitz bernoulli matix;
    A6= Mtx6(q6(1:samp_num),p6);
    
    message=sprintf('Sampling No.=%d',samp_num);
    disp(message);
    
 %   simu_time=zeros(1,6);
    
    for i=1:trial_num

        % create a sparse signal in Psi domain
        alp = [randn(k,1); zeros(N-k,1)];

        p = randperm(N);
        alp = alp(p);
        %x = Psi*alp;
        x=alp;
        
        % observation
        y1 = A1*x;
        y2 = A2*x;
        y5 = A5*x;
        y6 = A6*x;
        
        % reconstruction using OMP alg.
        sigma = 0;
        
 %       st=cputime;
        alp1 = omp(x, A1, y1, k, sigma);
        %xr1 = Psi*alp1;
        xr1=alp1;
  %      simu_time(1)=cputime-st+simu_time(1);
        
  %      st=cputime;
        alp2 = omp(x, A2, y2, k, sigma);
        %xr2 = Psi*alp2;
        xr2=alp2;
  %      simu_time(2)=cputime-st+simu_time(2);
        
  %      st=cputime;
        alp3 = omp(x, A3, y3, k, sigma);
        %xr3 = Psi*alp3;
        xr3=alp3;
  %      simu_time(3)=cputime-st+simu_time(3);
        
  %      st=cputime;
        alp5 = omp(x, A5, y5, k, sigma);
        %xr5 = Psi*alp5;
        xr5=alp5;
 %       simu_time(5)=cputime-st+simu_time(5);
        
 %       st=cputime;
        alp6 = omp(x, A6, y6, k, sigma);
        %xr6 = Psi*alp6;
        xr6=alp6;
 %       simu_time(6)=cputime-st+simu_time(6);
        
% if snr >50, it is considered as perfect reconstruction
        if snr(x,xr1)>50
            success_Gauss(j) = success_Gauss(j)+1;
        end
        
         if snr(x,xr2)>50
            success_Had_RP(j) = success_Had_RP(j)+1;
         end
            
          if snr(x,xr5)>50
            success_OSTM(j) = success_OSTM(j)+1;
          end
          
         if snr(x,xr6)>50
            success_Toep(j) = success_Toep(j)+1;
         end
    end
    
    success_Gauss(j) = success_Gauss(j)/trial_num;
    success_Had_RP(j) = success_Had_RP(j)/trial_num;
    success_OSTM(j) = success_OSTM(j)/trial_num;
    success_Toep(j) = success_Toep(j)/trial_num;
end 


figure;hold on;
plot([50:5:110], success_Gauss);
plot([50:5:110], success_Had_RP, 'x');
%plot([50:5:110]/256, success_blkHad, '+');
%plot([50:5:110]/256, success_blkOSTM, '*');
plot([50:5:110], success_OSTM, 's');
plot([50:5:110], success_Toep, 'd');
%legend('i.i.d Gaussian Operator','WHT256+Random Permutation','WHT8+Random Permutation','WHT256+Ramdom Sign Flipping','OSTM','Toeplitz');
legend('i.i.d Gaussian Operator','Hadamard256','OSTM256','Toeplitz Bernoulli');
xlabel('No. of Samples');
ylabel('Successful rate');
hold off;

