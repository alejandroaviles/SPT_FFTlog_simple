(* ::Package:: *)

Imatrix[nu1_,nu2_]:=1/(8 Pi^1.5) ((Gamma[3/2-nu1]Gamma[3/2-nu2]Gamma[nu1+nu2-3/2])/(Gamma[nu1]Gamma[nu2]Gamma[3-nu1-nu2]));
M22[nu1_,nu2_]:=Imatrix[nu1,nu2](3/2.-nu1-nu2)(1/2.-nu1-nu2) ( nu1 nu2 (98. (nu1+nu2)^2-14 (nu1+nu2) +36 )  -91 (nu1+nu2)^2 + 3 (nu1+nu2) +58  ) /  (  196  nu1(1+nu1)(1/2. - nu1)nu2(1+nu2)(1/2-nu2) );

M13[nu1_]:=(1+9. nu1)/4. Tan[nu1 * Pi]/(28.* Pi (nu1+1)nu1(nu1-1)(nu1-2)(nu1-3));


(* ::Section:: *)
(*M22 and M13*)


cmTetamTM[inputpkT_,kmin_,kmax_,Nfftlog_,biasnu_]:=Module[{kT,etaT,PrecmT,cmT,toFFTT,int,cmetamM,inputpk},
int=Log[kmax/kmin]/(Nfftlog-1);
kT=Table[kmin Exp[(i-1)int],{i,1,Nfftlog}];
inputpk=Interpolation[inputpkT];
toFFTT=Table[inputpk[kT[[i]]](kT[[i]]/kmin)^(-biasnu),{i,1,Nfftlog}] ;

etaT=Table[2 Pi (m-Nfftlog/2-1) (Nfftlog-1)/(Log[kmax/kmin](Nfftlog)),{m,1,Nfftlog+1}];
PrecmT=Fourier[toFFTT,FourierParameters->{-1, -1}];

cmT=Table[If[(j-Nfftlog/2)<1,kmin^(-biasnu- I etaT[[j]])Conjugate[PrecmT[[-j+Nfftlog/2+2]]],kmin^(-biasnu-I etaT[[j]])PrecmT[[j-Nfftlog/2]]],{j,1,Nfftlog+1}];
cmT[[1]]=cmT[[1]]/2;
cmT[[Length[cmT]]]=cmT[[Length[cmT]]]/2;
cmetamM=Table[{cmT[[j]],etaT[[j]]},{j,1,Nfftlog+1}];
cmetamM
];



M22M[kmin_,kmax_,Nfftlog_,biasnu_]:=Module[{nuT,M22matrix},
M22matrix=ConstantArray[0,{Nfftlog+1,Nfftlog+1}];
nuT=Table[-0.5(biasnu + I 2. Pi (m-Nfftlog/2.-1.) (Nfftlog-1.)/(Log[kmax/kmin](Nfftlog)) ),{m,1,Nfftlog+1}];
(*nu = -1/2 (biasnu + I eta) *)
Do[
Do[
M22matrix[[i,j]]=M22[nuT[[i]],nuT[[j]]];
,{j,1,Nfftlog+1}];
,{i,1,Nfftlog+1}];
M22matrix
];

M13M[kmin_,kmax_,Nfftlog_,biasnu_]:=Module[{cmTetamT,nuT,nu2T,M13vector},
M13vector=ConstantArray[0,{Nfftlog+1}];
nuT=Table[-0.5(biasnu + I 2. Pi (m-Nfftlog/2.-1.) (Nfftlog-1.)/(Log[kmax/kmin](Nfftlog)) ),{m,1,Nfftlog+1}];
(*nu = -1/2 (biasnu + I eta) *)
Do[M13vector[[i]]=M13[nuT[[i]]];,{i,1,Nfftlog+1}];
M13vector
];




P22M[kTout_,cmT_,etamT_,biasnu_,M22matrix_]:=Module[{k,vec,P22ofk,P22out},
P22out=ConstantArray[0,{Length@kTout}];
Do[
k=kTout[[i]];
vec=Table[cmT[[j]]k^(biasnu+ I etamT[[j]]),{j,1,Length@cmT}];
P22out[[i]]={k,Re[k^3 vec.M22matrix.vec]*Exp[-(k/(6*0.6711))^6]};
,{i,1,Length@kTout}];
P22out
]

P13UVM[kTout_,cmT_,etamT_,biasnu_,M13vector_,inputpkT_]:=Module[{k,inputpk,vec,sigma2v,P13ofk,P13out},
P13out=ConstantArray[0,{Length@kTout}];
(*sigma2v=Sum[inputpkT[[i,2]](inputpkT[[i,1]]-inputpkT[[i-1,1]]),{i,2,Length@inputpkT}]/(6. Pi^2);*)
inputpk=Interpolation[inputpkT];
sigma2v=NIntegrate[inputpk[p],{p,inputpkT[[1,1]],Last[inputpkT[[1]]]}]/(6. Pi^2);
Do[
k=kTout[[i]];
vec=Table[cmT[[j]]k^(biasnu+ I etamT[[j]]),{j,1,Length@cmT}];
P13out[[i]]={k,(Re[k^3 inputpk[k]vec.M13vector]-61/105  inputpk[k] k^2 sigma2v)*Exp[(k/(6*0.6711))^6](*exp(-pow(kdisc[index_j]/cutoff, 6.)*)};
,{i,1,Length@kTout}];
P13out
];

P13IRM[kTout_,cmT_,etamT_,biasnu_,M13vector_,inputpkT_]:=Module[{k,inputpk,vec,sigma2v,P13ofk,P13out},
P13out=ConstantArray[0,{Length@kTout}];
(*sigma2v=Sum[inputpkT[[i,2]](inputpkT[[i,1]]-inputpkT[[i-1,1]]),{i,2,Length@inputpkT}]/(6. Pi^2);*)
inputpk=Interpolation[inputpkT];
Do[
k=kTout[[i]];
vec=Table[cmT[[j]]k^(biasnu + I etamT[[j]]),{j,1,Length@cmT}];
P13out[[i]]={k,Re[k^3 inputpk[k]vec.M13vector](*-  inputpk[k] k^2 sigma2v*)};
,{i,1,Length@kTout}];
P13out
];


