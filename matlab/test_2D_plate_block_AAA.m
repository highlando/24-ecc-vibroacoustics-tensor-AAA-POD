% test_2D_plate_block_AAA.m
%
% Ion Victor Gosea, CSC Group, MPI Magdeburg
%
% The Block-AAA algorithm applied to 391 snapshots of sizes 3 x 3, 6 x 6,
% 12 x 12, or 21 x 21
%
% The required Matlab functions can be downloaded from the address below:
%
% https://github.com/nla-group/block_aaa
%
% Last modified: 22.10.2024
%

set(0,'DefaultFigurePosition', [100 100 1000 400]);
set(0,'defaultlinelinewidth',3)
set(0,'defaultlinemarkersize',20)
set(0,'defaultaxesfontsize',24)

clear all; close all;

podtnsdata = load('../data/podtnsrdata.mat');

largeorsmall = 4;

ftnsr = podtnsdata.full;  % the full data

for largeorsmall = 1:4
  switch largeorsmall
    case 1
      dtnsr = podtnsdata.pod3x3.dttnsr_ot;
      dvcso = podtnsdata.pod3x3.svdvecs_o;
      dvcsot = podtnsdata.pod3x3.svdvecs_ot;
      clgnd = '3x3 HOSVD'
      nxy = 3
    case 2
      dtnsr = podtnsdata.pod6x6.dttnsr_ot;
      dvcso = podtnsdata.pod6x6.svdvecs_o;
      dvcsot = podtnsdata.pod6x6.svdvecs_ot;
      clgnd = '6x6 HOSVD'
      nxy = 6
    case 3
      dtnsr = podtnsdata.pod12x12.dttnsr_ot;
      dvcso = podtnsdata.pod12x12.svdvecs_o;
      dvcsot = podtnsdata.pod12x12.svdvecs_ot;
      clgnd = '12x12 HOSVD'
      nxy = 12
    case 4 
      dtnsr = ftnsr;  
      clgnd = 'full data'
      nxy = 21
  end

  pts = (10:400);

  for ii = 1:length(pts)
      FF{ii} = dtnsr(:,:,ii);
  end

  size(FF)

  opts.tol = 1e-8;
  opts.maxit = 100;
  opts.return = 'best';

  tic;[R1,rmse1,out1] = block_aaa(FF,pts,opts); toc

  disp(['EXAMPLE 1 - best RMSE achieved: ' num2str(rmse(pts,R1,FF)) ])
  figure, semilogy(0:length(rmse1)-1,rmse1), hold on
  semilogy([0,opts.maxit-1],opts.tol*[1,1],'k--')
  xlabel('degree'), ylabel('RMSE')
  title('block-AAA convergence')

  %Here are the interpolation points (frequencies) that the algorithms
  %selects
  selected_pts = out1.zk.';

  aaatensor = zeros(nxy, nxy, 391);
  for w = pts
    aaatensor(:, :, w-9) = R1(w);
  end
  if not (largeorsmall == 4)
    prjdtnsro = nmodeproduct(aaatensor, dvcsot, 2);
    prjdtnsr = nmodeproduct(prjdtnsro, dvcso, 1);
    aaatensor = prjdtnsr;
  end

  figure(10101)
  semilogy(abs(ftnsr(:)-aaatensor(:))/norm(ftnsr(:)), 'DisplayName', clgnd)
  hold('on')
end
