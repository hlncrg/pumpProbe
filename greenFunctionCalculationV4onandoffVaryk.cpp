#include<iostream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<time.h>
#include<stdio.h>
#include<sstream>
#include<stdlib.h>
#include<complex>

using namespace std;

extern "C" void zgetrf_(int *m,int *n, complex<double> A[],int *lda, int tau[],int *info);
extern "C" void zgetri_(int *m,complex<double> A[], int *n, int tau[],complex<double> work[], int* p,int *info);

void Initialize();
void Finish();
void CalcBareGreens();
void CalcLocalGreens();
void CalcPhononGreens();
void CalcSelfEnergy();
void CalcGreensFunction();
void CalcPhotoemission();
double heavySide(double time,double timep);
complex<double> Energy(complex<double> kx,double ky);
complex<double> energyIntegral(int time,int timep,complex<double> kx,double ky);
complex<double> fermi(complex<double> energy,double mu);
complex<double> thetaPump(complex<double> time);
void invertc(double** A,double** Ainv,int N);
string toa(int num);
string toa(double num);

complex<double> **BareGreens;
complex<double> **LocalGreens;
complex<double> **PhononGreens;
complex<double> **SelfEnergy;
complex<double> *timer;
complex<double> **GreensFunction;
complex<double> **spectral;

double mu,thopping, tphopping,gsquared,temperature,e0;
int timeMax,timeMin,timeBeta,timeTotal;
double temperatureslice, timeslice, kmax,barekxMin,barekxMax,kindex,numbarekx;
double omegaPhonon,sigmaPump,timeOn,timeOff,sigmaProbe;
int omegaProbe,tmeas;
double omegaMax,omegaMin;
double PI=2.0*asin(1);
complex<double> eye(0,1);
double oneminusheavy;
double Nboson;
string outputFileName;



complex<double>* Ainv;
int* TAU;
int INFO;
complex<double>* WORK;
int lwork;




int main(){

Initialize();
CalcLocalGreens();
CalcPhononGreens();
CalcSelfEnergy();

for(kindex=barekxMin;kindex<barekxMax;kindex+=(barekxMax-barekxMin)/numbarekx){
CalcBareGreens();
CalcGreensFunction();
CalcPhotoemission();}
//calculate green functions in momentum space.

Finish();
return 0;}


void Initialize(){
/*
outputFileName="out.txt";
timeslice=1.0/8.0;
timeMin=-60;timeMax=90;timeBeta=20;temperature=.01;timeOn=timeslice*double(timeMin)/2.0;
timeTotal=2*(timeMax-timeMin)+timeBeta+1;
sigmaPump=250.0*timeslice*8.0;
omegaPhonon=0.1;
e0=7;mu=-0.255;
gsquared=.01;thopping=.25;tphopping=.3*thopping;
kmax=2;
omegaProbe=8;
sigmaProbe=625.0;
*/

//read in parameters.
cin>>outputFileName;
cin>>timeslice;
cin>>timeMin;
cin>>timeMax;
cin>>timeBeta;
cin>>temperature;
cin>>timeOn;
cin>>timeOff;
cin>>sigmaPump;
cin>>omegaPhonon;
cin>>e0;
cin>>mu;
cin>>gsquared;
cin>>thopping;
cin>>tphopping;
cin>>kmax;
cin>>omegaProbe;
cin>>omegaMin;
cin>>omegaMax;
cin>>sigmaProbe;
cin>>barekxMin;
cin>>barekxMax;
cin>>numbarekx;
cin>>tmeas;

outputFileName+=
"timeslice"+toa(timeslice)+
"timeMin"+toa(timeMin)+
"timeMax"+toa(timeMax)+
"timeBeta"+toa(timeBeta)+
"temperature"+toa(temperature)+
"timeOn"+toa(timeOn)+
"timeOff"+toa(timeOff)+
"sigmaPump"+toa(sigmaPump)+
"omegaPhonon"+toa(omegaPhonon)+
"e0"+toa(e0)+
"mu"+toa(mu)+
"gsquared"+toa(gsquared)+
"thopping"+toa(thopping)+
"tphopping"+toa(tphopping)+
"kmax"+toa(kmax)+
"omegaProbe"+toa(omegaProbe)+
"omegaMin"+toa(omegaMin)+
"omegaMax"+toa(omegaMax)+
"sigmaProbe"+toa(sigmaProbe)+
"barekxMin"+toa(barekxMin)+
"barekxMax"+toa(barekxMax)+
"tmeas"+toa(tmeas)+
"numbarekx"+toa(numbarekx)+
".txt";
//create name of output file


timeTotal=2*(timeMax-timeMin)+timeBeta+1;
lwork=timeTotal*timeTotal;
Nboson=1.0/(exp(omegaPhonon/temperature)-1.0);
//Nboson=0;
BareGreens=new complex<double>* [timeTotal];
for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
BareGreens[timeIndex]=new complex<double> [timeTotal];}

for(int timeIndex1=0;timeIndex1<timeTotal;timeIndex1++){
for(int timeIndex2=0;timeIndex2<timeTotal;timeIndex2++){
BareGreens[timeIndex1][timeIndex2]=0.0+0.0*eye;}}

LocalGreens=new complex<double>* [timeTotal];
for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
LocalGreens[timeIndex]=new complex<double> [timeTotal];}

for(int timeIndex1=0;timeIndex1<timeTotal;timeIndex1++){
for(int timeIndex2=0;timeIndex2<timeTotal;timeIndex2++){
LocalGreens[timeIndex1][timeIndex2]=0.0+0.0*eye;}}

PhononGreens=new complex<double>* [timeTotal];
for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
PhononGreens[timeIndex]=new complex<double> [timeTotal];}

for(int timeIndex1=0;timeIndex1<timeTotal;timeIndex1++){
for(int timeIndex2=0;timeIndex2<timeTotal;timeIndex2++){
PhononGreens[timeIndex1][timeIndex2]=0.0+0.0*eye;}}

SelfEnergy=new complex<double>* [timeTotal];
for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
SelfEnergy[timeIndex]=new complex<double> [timeTotal];}

for(int timeIndex1=0;timeIndex1<timeTotal;timeIndex1++){
for(int timeIndex2=0;timeIndex2<timeTotal;timeIndex2++){
SelfEnergy[timeIndex1][timeIndex2]=0.0+0.0*eye;}}

GreensFunction=new complex<double>* [timeTotal];
for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
GreensFunction[timeIndex]=new complex<double> [timeTotal];}

for(int timeIndex1=0;timeIndex1<timeTotal;timeIndex1++){
for(int timeIndex2=0;timeIndex2<timeTotal;timeIndex2++){
GreensFunction[timeIndex1][timeIndex2]=0.0+0.0*eye;}}

spectral=new complex<double>* [timeMax-timeMin];
for(int timeIndex=0;timeIndex<timeMax-timeMin;timeIndex++){
spectral[timeIndex]=new complex<double> [omegaProbe];}

timer=new complex<double> [timeTotal];
for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
if(timeIndex<(timeMax-timeMin)){timer[timeIndex]=timeslice*(timeMin+timeIndex);}
else{
if(timeIndex<2*(timeMax-timeMin)){timer[timeIndex]=timeslice*(2*timeMax-timeMin-timeIndex);}
else{timer[timeIndex]=timeslice*timeMin-eye/temperature*double((timeIndex-2*(timeMax-timeMin))/double(timeBeta));}
}}


Ainv=new complex<double> [timeTotal*timeTotal];
TAU=new int [timeTotal];
WORK=new complex<double> [timeTotal*timeTotal];


return;}
/************************************************/
void Finish(){

for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
delete [] BareGreens[timeIndex];}
delete [] BareGreens;

for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
delete [] LocalGreens[timeIndex];}
delete [] LocalGreens;

for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
delete [] PhononGreens[timeIndex];}
delete [] PhononGreens;

for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
delete [] SelfEnergy[timeIndex];}
delete [] SelfEnergy;

for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
delete [] GreensFunction[timeIndex];}
delete [] GreensFunction;

for(int timeIndex=0;timeIndex<timeMax-timeMin;timeIndex++){
delete [] spectral[timeIndex];}
delete [] spectral;

delete [] timer;


delete [] Ainv;
delete [] TAU;
delete [] WORK;


return;}
/************************************************/
void CalcBareGreens(){

for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
for(int timepIndex=0;timepIndex<timeTotal;timepIndex++){
cout<<"kx"<<kindex<<endl;
complex<double> timeValue=timer[timeIndex];
complex<double> timepValue=timer[timepIndex];
BareGreens[timeIndex][timepIndex]=eye*(fermi(Energy(PI*kindex,PI*kindex),mu)-heavySide(timeIndex,timepIndex))*
				exp(eye*mu*(timeValue-timepValue))*
				exp(-eye*energyIntegral(timeIndex,timepIndex,PI*kindex,PI*kindex));
}}


return;}
/************************************************/
void CalcLocalGreens(){

complex<double> timeValue,timepValue;
for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
for(int timepIndex=0;timepIndex<timeTotal;timepIndex++){
LocalGreens[timeIndex][timepIndex]=0.0+eye*0.0;
}}

for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
for(int timepIndex=0;timepIndex<timeTotal;timepIndex++){
for(int kxIndex=int(-kmax);kxIndex<kmax;kxIndex++){
for(int kyIndex=int(-kmax);kyIndex<kmax;kyIndex++){
timeValue=timer[timeIndex];
timepValue=timer[timepIndex];
LocalGreens[timeIndex][timepIndex]+=eye*(fermi(Energy(PI*kxIndex/kmax,PI*kyIndex/kmax),mu)-heavySide(timeIndex,timepIndex))*
                                exp(eye*mu*(timeValue-timepValue))*
				exp(-eye*energyIntegral(timeIndex,timepIndex,PI*kxIndex/double(kmax),PI*kyIndex/double(kmax)))
				/(4*kmax*kmax);
}}
}}


return;}
/************************************************/
void CalcPhononGreens(){

complex<double> timeValue,timepValue;

for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
for(int timepIndex=0;timepIndex<timeTotal;timepIndex++){
oneminusheavy=1-heavySide(timeIndex,timepIndex);
timeValue=timer[timeIndex];
timepValue=timer[timepIndex];
PhononGreens[timeIndex][timepIndex]=-eye*((Nboson+oneminusheavy)*exp(eye*omegaPhonon*(timeValue-timepValue)))
			 	   -eye*((Nboson+heavySide(timeIndex,timepIndex))*exp(-eye*omegaPhonon*(timeValue-timepValue)));

}}


return;}
/************************************************/
void CalcSelfEnergy(){
for(int timeIndex=0;timeIndex<timeTotal-1;timeIndex++){
for(int timepIndex=0;timepIndex<timeTotal-1;timepIndex++){
SelfEnergy[timeIndex][timepIndex]=eye*gsquared*LocalGreens[timeIndex][timepIndex]*PhononGreens[timeIndex][timepIndex]*(timer[timeIndex+1]-timer[timeIndex])*(timer[timepIndex+1]-timer[timepIndex]);
}
SelfEnergy[timeIndex][timeTotal-1]=eye*gsquared*LocalGreens[timeIndex][timeTotal-1]*PhononGreens[timeIndex][timeTotal-1]*(timer[timeIndex+1]-timer[timeIndex])*(timer[1]-timer[0]);
SelfEnergy[timeTotal-1][timeIndex]=eye*gsquared*LocalGreens[timeTotal-1][timeIndex]*PhononGreens[timeTotal-1][timeIndex]*(timer[timeIndex+1]-timer[timeIndex])*(timer[1]-timer[0]);
}

SelfEnergy[timeTotal-1][timeTotal-1]=eye*gsquared*LocalGreens[timeTotal-1][timeTotal-1]*PhononGreens[timeTotal-1][timeTotal-1]*(timer[1]-timer[0])*(timer[1]-timer[0]);

return;}
/************************************************/
void CalcGreensFunction(){

int bigIndex=0;
for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
for(int timepIndex=0;timepIndex<timeTotal;timepIndex++){
Ainv[bigIndex++]=BareGreens[timeIndex][timepIndex];}}

zgetrf_(&timeTotal,&timeTotal,Ainv,&timeTotal,TAU,&INFO);
zgetri_(&timeTotal,Ainv,&timeTotal,TAU,WORK,&lwork,&INFO);

bigIndex=0;
for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
for(int timepIndex=0;timepIndex<timeTotal;timepIndex++){
Ainv[bigIndex++]-=SelfEnergy[timeIndex][timepIndex];
}}

zgetrf_(&timeTotal,&timeTotal,Ainv,&timeTotal,TAU,&INFO);
zgetri_(&timeTotal,Ainv,&timeTotal,TAU,WORK,&lwork,&INFO);

bigIndex=0;
for(int timeIndex=0;timeIndex<timeTotal;timeIndex++){
for(int timepIndex=0;timepIndex<timeTotal;timepIndex++){
GreensFunction[timeIndex][timepIndex]=Ainv[bigIndex++];
}}


return;}
/************************************************/
void CalcPhotoemission(){
double spec;
spec=0.0;

ofstream outputFile;
outputFile.open(outputFileName.c_str(),fstream::app);

for(double omegaIndex=omegaMin;omegaIndex<omegaMax;omegaIndex+=(omegaMax-omegaMin)/double(omegaProbe)){spec=0.0;
for(int timeIndex=0;timeIndex<timeMax-timeMin;timeIndex++){
for(int timepIndex=timeMax-timeMin;timepIndex<2*(timeMax-timeMin);timepIndex++){
spec+=imag(1/PI*GreensFunction[timeIndex][timepIndex]*
      exp(-pow(timer[timeIndex]-timer[tmeas-timeMin],2)/(2.0*pow(sigmaProbe,2)))/(sigmaProbe*sqrt(2*PI))*
      exp(-pow(timer[timepIndex]-timer[tmeas-timeMin],2)/(2.0*pow(sigmaProbe,2)))/(sigmaProbe*sqrt(2*PI))*
                                    exp(eye*omegaIndex*(timer[timeIndex]-timer[timepIndex]))*
                                    timeslice*timeslice);


}}
outputFile<<setprecision(15)<<kindex*PI<<" "<<omegaIndex<<" "<<spec<<endl;
}outputFile<<endl;



outputFile.close();
return;}
/************************************************/
double heavySide(double time,double timep){
if(time>timep)return 1;
return 0;}
/************************************************/
complex<double> Energy(complex<double> kx,double ky){
return -2.0*thopping*(cos(kx)+cos(ky))+4.0*tphopping*cos(kx)*cos(ky);}
/************************************************/
complex<double> energyIntegral(int time,int timep,complex<double> kx,double ky){
complex<double> sum;

if(time>timep){for(int timeIndex=timep;timeIndex<time;timeIndex++){
sum=sum+Energy(kx+e0*thetaPump(timer[timeIndex])*timer[timeIndex],ky)*(timer[timeIndex+1]-timer[timeIndex]);
}}  
else{for(int timeIndex=time;timeIndex<timep;timeIndex++){
sum=sum-Energy(kx+e0*thetaPump(timer[timeIndex])*timer[timeIndex],ky)*(timer[timeIndex+1]-timer[timeIndex]);
}}

return sum;}
/************************************************/
complex<double> fermi(complex<double> Energy,double mu){
return 1.0/(1.0+exp((Energy-mu)/temperature));}
/************************************************/
complex<double> thetaPump(complex<double> time){
if(real(time)>=timeOn && real(time)<=timeOff){
return 1.0;}
return 0.0;}
/************************************************/
string toa(int num){
stringstream converter;
converter << num;
return converter.str();}

string toa(double num){
stringstream converter;
converter << num;
return converter.str();}


