(* ::Package:: *)

(* ::Title:: *)
(*Module to get pk_nw *)


(* ::Text:: *)
(*Alejandro Aviles*)
(*avilescervantes@gmail.com*)
(**)
(*call *)
(**)
(*pknwM[kmin, kmax, Nk, inputPkLinear, hubble] *)
(**)
(*to get the non-wiggle power spectrum with Nk points spaced in log-intervals between [kmin,kmax].*)
(**)
(*You have to add the hubble = H0 /(km/s/Mpc)*)
(**)
(*inputPklinear is a two-columns table { k [h/Mpc], Pk [Mpc/h]^3 }    *)
(**)
(**)
(**)


pknwM[outputkmin_,outputkmax_,outputsize_,inputPklTable_,h_]:=Module[
{ksT,ksmin,ksmax,Nks,logkpkT,PSL,FSTlogkpkT,FSTlogkpkEvenT,FSTlogkpkOddT,mcutmin,mcutmax,nFSTlogkpkEvenTcuttedT,nFSTlogkpkOddcutted,nFSTlogkpkOddTcuttedT,nFSTlogkpkEvencutted,preT,FSTofFSTlogkpkNWT,PNW1T,PNW1,DeltaAppf,PNW,outkT,FSTtype,intorder,result},

ksmin=7*10^(-5)/h;
ksmax=7./h;
Nks=2^16;

ksT=Table[ksmin +(ii-1) (ksmax-ksmin)/(Nks-1),{ii,1,Nks}];
PSL=Interpolation[inputPklTable];
logkpkT=Table[Log[ksT[[i]] PSL[ksT[[i]]]],{i,1,Length@ksT}];

(* Add here: check sizes of table,  abort if dont match *)
FSTtype=1;
intorder=4;

FSTlogkpkT=FourierDST[logkpkT,FSTtype];
FSTlogkpkEvenT=Table[FSTlogkpkT[[m]],{m,2,Length@FSTlogkpkT,2}];
FSTlogkpkOddT=Table[FSTlogkpkT[[m]],{m,1,Length@FSTlogkpkT,2}];


(* cut tables:*)
mcutmin=120;
mcutmax=240;


nFSTlogkpkEvenTcuttedT=Join[Table[{m,FSTlogkpkEvenT[[m]]},{m,2,mcutmin-2}],Table[{m,FSTlogkpkEvenT[[m]]},{m,mcutmax+2,Length@FSTlogkpkEvenT}]];
nFSTlogkpkOddTcuttedT=Join[Table[{m,FSTlogkpkOddT[[m]]},{m,1,mcutmin-1}],Table[{m,FSTlogkpkOddT[[m]]},{m,mcutmax+1,Length@FSTlogkpkOddT}]];

nFSTlogkpkOddcutted=Interpolation[nFSTlogkpkOddTcuttedT,InterpolationOrder->intorder(*,Method->"Spline"*)];
nFSTlogkpkEvencutted=Interpolation[nFSTlogkpkEvenTcuttedT,InterpolationOrder->intorder(*,Method->"Spline"*)];

preT=ConstantArray[0,Length@FSTlogkpkT];
Do[
If[mcutmin<i<mcutmax,preT[[2*i]]={2*i,nFSTlogkpkEvencutted[i+1]};preT[[2*i-1]]={2*i-1,nFSTlogkpkOddcutted[i]}];
If[mcutmin>=i ||i>=mcutmax,preT[[2*i]]={2*i,FSTlogkpkT[[2*i]]};preT[[2*i-1]]={2*i-1,FSTlogkpkT[[2*i-1]]}];
,{i,1,Length@FSTlogkpkT/2}];

(* Inverse Sine transformation*)
FSTofFSTlogkpkNWT=FourierDST[preT[[All,2]],FSTtype];

PNW1T=Table[{ksT[[m]],Exp[FSTofFSTlogkpkNWT[[m]]]/ksT[[m]]},{m,1,Length@ksT}];
PNW1=Interpolation[PNW1T];

DeltaAppf[k_]:=k (PSL[ksT[[8]]]-PNW1[ksT[[8]]])/PNW1[ksT[[8]]]/ksT[[8]];

PNW[k_]:=Piecewise[{{PSL[k]/(DeltaAppf[k]+1),k<0.001},{PNW1[k],0.001<=k<=Last[ksT]},{PSL[k],k> Last[ksT]}}];

outkT=Table[outputkmin Exp[(i-1)Log[outputkmax/outputkmin]/(outputsize-1)],{i,1,outputsize}];
result=Table[{outkT[[i]],PNW[outkT[[i]]]},{i,1,Length@outkT}];
result
];
