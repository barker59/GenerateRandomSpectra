%++++++++++++IN THE NAME OF GOD THE COMPASSIONATE THE MERCIFUL++++++++++++%
% Code to calculate wave properties using wave surface elevation power    % 
% spectral density                                                        %
% Ver 1.3                                                                 %
%                                                     by: Arash Karimpour %
%                                                   www.arashkarimpour.com%
%+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++%

%INPUT---------------------------------------------------------------------
% calculating wave properties from power spectral density
%
% input=importdata('d.mat');      %load water depth (d)/surface elevation (Eta) data and rename it "input" in (m)
% fs=10;                          %sampling frequency that data collected at in (Hz)
% duration=1024;                  %duration time that data collected in input in each burst in second
% nfft=2^10;                      %NFFT for Fast Fourier Transform
% heightfrombed=0.0;              %sensor height from bed
% fmin=0.04;                      %minimum frequency for cut off the lower part of spectra
% fmax=1;                         %maximum frequency for cut off the upper part of spectra
% ftailcorrection=1;              %frequency that diagnostic tail apply after that (typically set at 2.5fm, fm=1/Tm01)
% tailpower=-4;                   %power that diagnostic tail apply based on that (-3 for shallow water to -5 for deep water)
% mincutoff=1;                    %define if to cut off the spectra below fmin or not (0: cutoff off, 1: cutoff on)
% maxcutoff=1;                    %define if to cut off the spectra beyond fmax or not (0: cutoff off, 1: cutoff on)
% tailcorrection=0;               %define if apply tail correction or not (0: not apply, 1: diagnostic tail, 2: TMA Spectrum tail)
% ploton=1;                       %define if to plot results or not (0: not plot, 1: plot)

%OUTPUT--------------------------------------------------------------------

% Hm0                             : Zero-Moment Wave Height (m)
% Tm01                            : Wave Period from m01 (second), Mean Wave Period
% Tm02                            : Wave Period from m02 (second), Mean Zero Crossing Period
% Tp                              : Peak Wave Period (second)
% fp                              : Peak Wave Frequency (Hz)
% f                               : Frequency (Hz)
% Syy                             : Wave Surface Elevation Power Spectrum (m^2s)

%--------------------------------------------------------------------------

function [Hm0,Tm01,Tm02,Tp,Te,fp,f,Syy]=WaveSpectraFun(input,fs,duration,nfft,heightfrombed,fmin,fmax,ftailcorrection,tailpower,mincutoff,maxcutoff,tailcorrection,ploton)

%--------------------------------------------------------------------------
%deterending

input1=detrend(input,'linear');

%--------------------------------------------------------------------------
%calculating Fast Fourier transform and power density

fmax(fmax>fs/2)=fix(fs/2);
%nfft = 2^(nextpow2(length(input1)));

%--------------------------------------------------------------------------
clear j1 j2 i

%calculating power density 

clear f w HS hpsd Syy

windowelem = 256; % number of elements in each window
overlapelem = 128; %number of overlap element

[Syy,f] = pwelch(input1,hanning(windowelem),overlapelem,nfft,fs); %Wave power spectrum and Frequency

w=2*pi*f; %angular velocity
deltaf=f(2,1)-f(1,1);

%--------------------------------------------------------------------------
% applying tail correction

if tailcorrection==1 % applying diagnostic frequency tail based on SWAN after fmax
    
    fmaxloc=find(f>=ftailcorrection);
    for i=1:length(Syy(:,1))
        if f(i,1)>ftailcorrection
            Syy(i,1)=Syy(fmaxloc(1,1),1)*(f(i,1)/ftailcorrection).^tailpower; %Adding diagnostic tail
        end
    end
    Syy(Syy<0)=0; %Syy can not be negative
    
end

if tailcorrection==2 % applying TMA Spectrum tail after fmax
    
    omega=2*pi.*f.*sqrt(h/9.81);
    
    THETA=ones(length(omega(:,1)),1);
    THETA(omega<=1)=omega(omega<=1).^2./2;
    THETA(omega>1)=1-0.5*(2-omega(omega>1)).^2;
    
    fmaxloc=find(f>=ftailcorrection);
    [THETAmax THETAmaxloc]=max(THETA(:,1));
    for i=1:length(Syy(:,1))
        if f(i,1)>ftailcorrection
            if f(i,1)>=f(THETAmaxloc,1)
                Syy(i,1)=Syy(fmaxloc(1,1),1)*THETA(i,1)*(f(i,1)/ftailcorrection).^tailpower; %Adding TMA Spectrum tail
            elseif f(i,1)<f(THETAmaxloc,1)
                Syy(i,1)=Syy(fmaxloc(1,1),1)*(f(i,1)/ftailcorrection).^tailpower; %Adding JONSWAP Spectrum tail
            end
        end
    end
    Syy(Syy<0)=0; %Syy can not be negative
    
end

%--------------------------------------------------------------------------
%cut off spectra based on fmin and fmax 

if mincutoff==1
    Syy(f<fmin)=0;
end

if maxcutoff==1
    Syy(f>fmax)=0;
end

%--------------------------------------------------------------------------

%Calculating spectral moments
clear m0 m1 m2 m_1
m_1=nansum(Syy.*f.^-1*deltaf);
setappdata(0,'m_1',m_1);
m0=sum(Syy.*f.^0*deltaf);
m1=sum(Syy.*f.^1*deltaf);
m2=sum(Syy.*f.^2*deltaf);

% calculating wave properties
Hm0=4*sqrt(m0); %Zero-Moment wave height
Tm01=m0/m1; %mean period
Tm02=(m0/m2)^0.5; %zero crossing period

% calculation peak period
[Syymax loc4]=max(Syy(:,1));
Tp=1/f(loc4,1); %peak period
Te=m_1/m0;

% calculating peak frequency from weighted integral (Young, 1995)
fp=(sum(Syy.^5.*f.^1*deltaf))./(sum(Syy.^5.*f.^0*deltaf)); %peak frequency
%--------------------------------------------------------------------------
%plotting

if ploton==1
    
    val=[m0 m1 m2 Hm0 Tm01 Tm02 Tp Te fp];
    name={'m0','m1','m2','Hm0','Tm01','Tm02','Tp','Te','fp'};
    for i=1:length(val)
        fprintf('%14s   %g\n',name{i},val(i));
    end
    
    %plotting
    loglog(f(f~=0),Syy(f~=0))
    hold on
    
    title('Power Spectral Density')
    xlabel('Frequency(Hz)')
    ylabel('Spectral Density(m^2s)')
    
end
%--------------------------------------------------------------------------

