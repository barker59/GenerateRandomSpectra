function [Hm0,Tm01,Tm02,Tp,Te,fp,f,Syy] = Random_Spectrum_Oceanlyz(d)
%%Performs Spectral Analysis using Wave Spectra Fun from Oceanlyz
%%(http://www.arashkarimpour.com/download.html)
%Created by Aaron Barker - 20/11/2015

%d=waterdepthsample %Data import
module=1;
%INPUT---------------------------------------------------------------------

%measurement properties
burst=1;               %number of burst in the input file
duration=2048;         %duration time that data collected in each burst in second
nfft=2^11;             %NFFT for Fast Fourier Transform
fs=1;                 %sampling frequency that data collected at in (Hz)
heightfrombed=0.05;    %Pressure sensor height from bed in (m)

%input values
fmin=0.04;             %minimum frequency for cut off the lower part of spectra (Hz)
fmax=5;                %maximum frequency for cut off the upper part of spectra (Hz)
fminpcorr=0.15;        %minimum frequency that automated calculated fmaxpcorr can have if autofmaxpcorr=1 (Hz)
fmaxpcorr=0.55;        %maximum frequency for applying pressure attenuation factor (Hz)
ftailcorrection=0.9;   %frequency that diagnostic tail apply after that (typically set at 2.5fm, fm=1/Tm01) (Hz)
tailpower=-5;          %power that diagnostic tail apply based on that (-3 for shallow water to -5 for deep water)
fminswell=0.1;         %minimum frequency that swell can have (it is used for Tpswell calculation) (Hz)
fmaxswell=0.25;        %maximum frequency that swell can have, It is about 0.2 in Gulf of Mexico (Hz)

%switching module on or off
pressureattenuation=2; %define if to apply pressure attenuation factor or not (0:off, 1:on, without correction after fmaxpcorr, 2:on, with constant correction after fmaxpcorr)
autofmaxpcorr=1;       %define if to calculate fmaxpcorr and ftailcorrection based on water depth or not (0:off, 1:on)
mincutoff=1;           %define if to cut off the spectra below fmin (0: cutoff off, 1: cutoff on)
maxcutoff=0;           %define if to cut off the spectra beyond fmax (0: cutoff off, 1: cutoff on)
tailcorrection=0;      %define if to apply tail correction or not (0: not apply, 1: diagnostic tail, 2: TMA Spectrum tail)
ploton=0;              %define if to plot spectra or not (0: not plot, 1: plot)
saveon=1;              %define if to save data or not (0: not save, 1: save)

%FUNCTION------------------------------------------------------------------
% Calling calculation function

%[wave]=CalcWaveFun(inputfilelocation,inputfilename,module,burst,duration,nfft,fs,heightfrombed,fmin,fmax,fminpcorr,fmaxpcorr,ftailcorrection,tailpower,fminswell,fmaxswell,pressureattenuation,autofmaxpcorr,mincutoff,maxcutoff,tailcorrection,ploton);
sample=fs*duration; %number of sample in 1 burst
fmaxpcorr1=fmaxpcorr; %storing user defined fmaxpcorr
pressureattenuation=0; %define if pressure attenuation factor apply (0:off, 1:on)

%loading data
%currentpath=cd(inputfilelocation);
%d=importdata(inputfilename);
i=burst;
j1=(i-1)*sample+1
j2=i*sample
input=d(j1:j2,1);
h=mean(input(:,1)); % calculating mean water depth
if h<=0
    warning('Mean water depth is Zero or negative, Oceanlyz continues with mean water depth=0.001 m.');
    h=0.001;
end
	
%Spectral Analysis Function
[wave.Hm0(i,1),wave.Tm01(i,1),wave.Tm02(i,1),wave.Tp(i,1),wave.Te(i,1),wave.fp(i,1),wave.f(i,:),wave.Syy(i,:)]=WaveSpectraFun(input,fs,duration,nfft,heightfrombed,fmin,fmax,ftailcorrection,tailpower,mincutoff,maxcutoff,tailcorrection,ploton);
Hm0=wave.Hm0;
Tm01=wave.Tm01;
Tm02=wave.Tm02;
Tp=wave.Tp;
Te=wave.Te;
fp=wave.fp;
f=wave.f;
Syy=wave.Syy;
 