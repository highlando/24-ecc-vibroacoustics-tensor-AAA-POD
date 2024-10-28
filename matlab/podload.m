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
end

% inflate the (reduced) tensors to full size
prjdtnsro = nmodeproduct(dtnsr, dvcsot, 2);
prjdtnsr = nmodeproduct(prjdtnsro, dvcso, 1);

% compare to full data
norm(ftnsr(:)-prjdtnsr(:), 'fro')
