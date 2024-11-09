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


close all;


if exist('clsp')
  largeorsmall = clsp;
else
  clear all;
  largeorsmall = 3;
end

podtnsdata = load('../data/podtnsrdata.mat');
ftnsr = podtnsdata.full;  % the full data

switch largeorsmall
  case 1
    dtnsr = podtnsdata.pod3x3.dttnsr_ot;
    dvcso = podtnsdata.pod3x3.svdvecs_o;
    dvcsot = podtnsdata.pod3x3.svdvecs_ot;
    nx = 3;
  case 2
    dtnsr = podtnsdata.pod6x6.dttnsr_ot;
    dvcso = podtnsdata.pod6x6.svdvecs_o;
    dvcsot = podtnsdata.pod6x6.svdvecs_ot;
    nx = 6;
  case 3
    dtnsr = podtnsdata.pod12x12.dttnsr_ot;
    dvcso = podtnsdata.pod12x12.svdvecs_o;
    dvcsot = podtnsdata.pod12x12.svdvecs_ot;
    nx = 12;
  case 4 
    dtnsr = ftnsr;  
    nx = 21;
end

pts = (10:400);

for ii = 1:length(pts)
    FF{ii} = dtnsr(:,:,ii);
end

size(FF)

% set tolerance here
% opts.tol = 1e-6;
% opts.tol = 1e-8;
%opts.tol = 1e-10;

tolvec = [1e-4, 1e-6, 1e-8, 1e-10];

for ctol = tolvec
  opts.tol = ctol;

  % set maximum numbe rof iterations (100 to make sure it stops at the
  % tolerancevalue)
  opts.maxit = 100;

  opts.return = 'best';

  % this is Block-AAA
  tic;[R1,rmse1,out1] = block_aaa(FF,pts,opts); toc

  ord1 = length(out1.zk);

  % this is SetValued-AAA
  tic;[R2,rmse2,zk2] = set_val_aaa(FF,pts,opts.maxit,opts.tol);toc

  ord2 = length(zk2);

  %Here are the interpolation points (frequencies) that the algorithm
  %Block-AAA selects
  selected_pts = out1.zk.';

  for ii = 1:length(pts)
  %     err1(ii)= norm(FF{ii}-R1(pts(ii)));
  %     err2(ii)= norm(FF{ii}-R2(pts(ii)));
      Fap1 = dvcso * R1(pts(ii)) * dvcsot';
      Fap2 = dvcso * R2(pts(ii)) * dvcsot';
      Fex = ftnsr(:, :, ii);
      err1(ii)= norm(Fex-Fap1, 'fro');
      err2(ii)= norm(Fex-Fap2, 'fro');
  end

  err1 = zeros(length(pts), 1);
  err2 = zeros(length(pts), 1);
  for ii = 1:length(pts)
      app1 = dvcso * R1(pts(ii)) * dvcsot';
      app2 = dvcso * R2(pts(ii)) * dvcsot';
      err1(ii) = norm(ftnsr(:, :, ii)-app1, 'fro');
      err2(ii) = norm(ftnsr(:, :, ii)-app2, 'fro');
  end
  fprintf('blk: nx: %d, tol: %.1e, m: %d, err: %.2e\n', ... 
          nx, ctol, ord1-1, norm(err1))
  fprintf('set: nx: %d, tol: %.1e, m: %d, err: %.2e\n', ...
          nx, ctol, ord2-1, norm(err2))

end

figure
semilogy(err1,'b*','markersize',14); hold on;
semilogy(err2,'ro','markersize',14);
legend('Block-AAA','SetValued-AAA');
ylabel('Error (2-norm of the deviation)');
xlabel('Frequency');

% disp(['EXAMPLE 1 - best RMSE achieved: ' num2str(rmse(pts,R1,FF)) ])
figure, semilogy(0:length(rmse1)-1,rmse1), hold on
semilogy([0,opts.maxit-1],opts.tol*[1,1],'k--');
semilogy(0:length(rmse2)-1,rmse2,'-.r');
xlabel('degree'), ylabel('RMSE')
legend('Block-AAA convergence','Set-valued AAA convergence')
