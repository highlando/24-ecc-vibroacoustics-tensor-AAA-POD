% test_2D_plate_block_AAA.m
%
% Ion Victor Gosea, CSC Group, MPI Magdeburg
%
% The Block-AAA and SetValued-AAA algorithms applied to 391 snapshots of 
% sizes 3 x 3, 6 x 6, 12 x 12, or 21 x 21
%
% The required Matlab functions can be downloaded from the address below:
%
% https://github.com/nla-group/block_aaa
%
% Last modified: 08.11.2024
%


clear all; close all;

podtnsdata = load('../data/podtnsrdata.mat');

largeorsmall = 3;

ftnsr = podtnsdata.full;  % the full data

switch largeorsmall
  case 1
    dtnsr = podtnsdata.pod3x3.dttnsr_ot;
    dvcso = podtnsdata.pod3x3.svdvecs_o;
    dvcsot = podtnsdata.pod3x3.svdvecs_ot;
  case 2
    dtnsr = podtnsdata.pod6x6.dttnsr_ot;
    dvcso = podtnsdata.pod6x6.svdvecs_o;
    dvcsot = podtnsdata.pod6x6.svdvecs_ot;
  case 3
    dtnsr = podtnsdata.pod12x12.dttnsr_ot;
    dvcso = podtnsdata.pod12x12.svdvecs_o;
    dvcsot = podtnsdata.pod12x12.svdvecs_ot;
  case 4 
    dtnsr = ftnsr;  
end

pts = (10:400);

for ii = 1:length(pts)
    FF{ii} = dtnsr(:,:,ii);
end

size(FF)

% set tolerance here
%opts.tol = 1e-4;
opts.tol = 1e-6;
% opts.tol = 1e-8;
%opts.tol = 1e-10;

% set maximum numbe rof iterations (100 to make sure it stops at the
% tolerancevalue)
opts.maxit = 100;

opts.return = 'best';

% this is Block-AAA
tic;[R1,rmse1,out1] = block_aaa(FF,pts,opts); toc

ord1 = length(out1.zk)

% this is SetValued-AAA
tic;[R2,rmse2,zk2] = set_val_aaa(FF,pts,opts.maxit,opts.tol);toc

ord2 = length(zk2)


disp(['EXAMPLE 1 - best RMSE achieved: ' num2str(rmse(pts,R1,FF)) ])
figure, semilogy(0:length(rmse1)-1,rmse1), hold on
semilogy([0,opts.maxit-1],opts.tol*[1,1],'k--');
semilogy(0:length(rmse2)-1,rmse2,'-.r');
xlabel('degree'), ylabel('RMSE')
legend('Block-AAA convergence','Set-valued AAA convergence')

%Here are the interpolation points (frequencies) that the algorithm
%Block-AAA selects
selected_pts = out1.zk.';

for ii = 1:length(pts)
    err1(ii)= norm(FF{ii}-R1(pts(ii)));
    err2(ii)= norm(FF{ii}-R2(pts(ii)));
end

figure
semilogy(err1,'b*','markersize',14); hold on;
semilogy(err2,'ro','markersize',14);
legend('Block-AAA','SetValued-AAA');
ylabel('Error (2-norm of the deviation)');
xlabel('Frequency');

AAAcoretensor = zeros(size(dtnsr));
for ii = 1:length(pts)
  AAAcoretensor(:, :, ii) = R1(pts(ii));
end

prjdtnsro = nmodeproduct(AAAcoretensor, dvcsot, 2);
prjdtnsr = nmodeproduct(prjdtnsro, dvcso, 1);
% compare to full data
norm(ftnsr(:)-prjdtnsr(:), 'fro')
