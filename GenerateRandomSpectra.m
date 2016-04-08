%% 'Generate Random Spectra'
% Written by Aaron Barker, 6/4/2016

% Creates spectra for analysis by students of CE4013 Harbour and Coastal
% Engineering. 
% Output consists of a '.csv' file of frequency and spectral
% ordinates, '.png' image of the spectrum (for student),'.png' image of
% spectrum with results of analysis overlaid (for lecturer) and a Summary
% 'SpectrumSummary(Lecturer).csv' file for the lecturer containing
% frequency and spectral ordinates alongside the correct results.

%% Setup
no_of_spectra = 11;
spectraGenerationMethod = 'Empirical Spectral Formulation';
SaveFileName='SpectrumSummary(Lecturer).csv'
myfile = fopen(SaveFileName ,'w'); %Clear the file before use
fclose(myfile)

%% Generate waves in the time domain and convert to spectrum
if strcmp(spectraGenerationMethod,'Time Domain') == 1
    frequency_step = 0.0156;
    end_frequency = frequency_step *2048;
    frequency = linspace(0,end_frequency,end_frequency/frequency_step);

    WavePeriods = [10,6,7,8,12]
    WaveMagnitudes = [6,2,8,2,2]
    k=1
    for t=1:1:10000
        for j=1:length(WavePeriods)
            wave(k,j) = (WaveMagnitudes(j) * sin(2*pi*(1/WavePeriods(j))*t));
        end
        k=k+1;
    end

    [Hm0,Tm01,Tm02,Tp,Te,fp,f,Syy] = Random_Spectrum_Oceanlyz(wave(:,1))
    %plot(f,Syy)
    clear('wave')
end

%% Generate in the frequency domain via empirical spectra
if strcmp(spectraGenerationMethod,'Empirical Spectral Formulation') == 1
    %Wave magnitudes and periods can be replaced with random number
    %generation, linspace, or other function driven generation
    WavePeriods = [10,9,8,7,7.12,5.54,12,15,13.43,12.1,5];
    WaveMagnitudes = [4,3.2,2.7,2,2,1.23,3,5,5,4,12];
    SecondaryWaveMagnitudes = WaveMagnitudes;
    SecondaryWavePeriods = WavePeriods*0.5;
    SpectraType = {'bs','bs','bs','bs','bs','bs','bs','bs','bs','bs','bs'};
    noModes = {2,1,1,1,2,1,1,1,1,2,1}; %'1' creates a unimodal spectrum, '2' results in a bimodal spectrum with second period shifted by 50%
    InfoGrav = {0,0,1,0,0,0,1,0,0,0,0}; %'1' Flag adds a constant InfaGravity wave to the spectrum
    for i=1:no_of_spectra
        
        %Generate Frequency ordinates
        frequency_step = 0.001;
        end_frequency = 0.6;
        frequency = linspace(0,end_frequency,2048);

        %Generate Spectra based on choice of spectral type
        if strcmp(SpectraType{i},'bs') == 1
            %Generate an empirical Bretschneider Spectrum
            if noModes{i} == 1
                SpectraName = 'Bretschneider'
                [SpectralOrdinates,Spectra_freq] = generateBretschneiderSpectrum(WaveMagnitudes(i),WavePeriods(i),frequency);
                SpectralOrdinates(isnan(SpectralOrdinates))=0;
                Spectra_df = Spectra_freq(3,1)-Spectra_freq(2,1);
            else
                SpectraName = 'Multi-Modal Bretschneider'
                [SpectralOrdinates1,Spectra_freq] = generateBretschneiderSpectrum(WaveMagnitudes(i),WavePeriods(i),frequency);
                [SpectralOrdinates2,Spectra_freq] = generateBretschneiderSpectrum(SecondaryWaveMagnitudes(i),SecondaryWavePeriods(i),frequency);
                SpectralOrdinates = SpectralOrdinates1 + SpectralOrdinates2;
                SpectralOrdinates(isnan(SpectralOrdinates))=0;
                Spectra_df = Spectra_freq(3,1)-Spectra_freq(2,1);
            end
            if InfoGrav{i} == 1
               %Add an InfaGravity Wave
               % Hard coded to keep external files to a minimum for this
               % quick project.
               SpectralOrdinates(20:57,1)= [0.0000	0.0000	0.0065	0.0330	0.0916	0.1944	0.3465	0.5431	0.7788	1.0481	1.3442	1.6491	1.9391	2.1905	2.3794	2.4899	2.5272	2.5005	2.4187	2.2909	2.1257	1.9314	1.7164	1.4889	1.2573	1.0296	0.8140	0.6185	0.4512	0.3165	0.2115	0.1326	0.0761	0.0385	0.0159	0.0045	0.0005	0.0000]

            end
        elseif strcmp(SpectraType{i},'JS') == 1
            %Generate am empirical Jonswap Spectrum
            %Work TBC
            SpectraName = 'Jonswap';
        else
            %Generate am empirical Pierson-Moskowitz Spectrum
            %Work TBC
            SpectraName = 'Pierson Moskowitz';    
        end

        %Plot the spectrum
        for plotrun =1:2
            figure('units','normalized','outerposition',[0 0 1 1])
            p = plot(Spectra_freq,SpectralOrdinates)
            xlim([0,0.6])

            if noModes{i} == 1
                plotname = ['Spectrum (',num2str(i),') - ',SpectraName,' -  with ',num2str(WaveMagnitudes(i)),'m Hm0 and ',num2str(WavePeriods(i)),'s T02']
            else
                plotname = {['Spectrum (',num2str(i),') - ',SpectraName,' -  with ',num2str(WaveMagnitudes(i)),'m Hm0 and ',num2str(WavePeriods(i)),'s T02'],[' PLUS secondary spectrum with ',num2str(SecondaryWaveMagnitudes(i)),'m Hm0 and ',num2str(SecondaryWavePeriods(i)),'s T02']}
            end
            if plotrun ==1
                title(['Spectrum (',num2str(i),')'],'fontweight','bold')
            else
                title(plotname,'fontweight','bold')
            end
            ylabel('Spectral Energy (m^2/Hz)','fontweight','bold')
            xlabel('Frequency (Hz)','fontweight','bold')

            % Check that the spectra was appropriately created using spectral moment
            % analysis to determine key parameters

            l=length(SpectralOrdinates)
            for k = 1:l;
                m0_temp(k) = (Spectra_df*SpectralOrdinates(k)*Spectra_freq(k,1)^0);
            end
            Spectra_m0 = sum(m0_temp)

            for k = 1:l;
                m_1_temp(k) = (Spectra_df*SpectralOrdinates(k)*Spectra_freq(k,1)^-1);
            end
            Spectra_m_1 = nansum(m_1_temp)


            for k = 1:l;
                m1_temp(k) = (Spectra_df*SpectralOrdinates(k)*Spectra_freq(k,1)^1);
            end
            Spectra_m1 = nansum(m1_temp)

            for k = 1:l;
                m2_temp(k) = (Spectra_df*SpectralOrdinates(k)*Spectra_freq(k,1)^2);
            end
            Spectra_m2 = nansum(m2_temp)

            for k = 1:l;
                m_2_temp(k) = (Spectra_df*SpectralOrdinates(k)*Spectra_freq(k,1)^-2);
            end
            Spectra_m_2 = nansum(m_2_temp)

            for k = 1:l;
                m4_temp(k) = (Spectra_df*SpectralOrdinates(k)*Spectra_freq(k,1)^4);
            end
            Spectra_m4 = nansum(m4_temp)



            % Spectral Parameters from spectral moment analysis
            Spectra_Hm0(i) = 4*sqrt(Spectra_m0);
            Spectra_T02(i) = sqrt(Spectra_m0/Spectra_m2);
            Spectra_Te(i) = (Spectra_m_1/Spectra_m0);
            Spectra_Tpc(i) = (Spectra_m_2*Spectra_m1/Spectra_m0^2);
            Spectra_Tc(i) = (Spectra_m2/Spectra_m4)^.5;

            dim = [0.5 0.4 0.25 0.2];
            str = {'Results of spectral calculation:',['Hm0: ',num2str(Spectra_Hm0(i)),' m'],...
                ['T02: ',num2str(Spectra_T02(i)),' s'],['Te: ',num2str(Spectra_Te(i)),' s'],['Tp (Calculated): ',num2str(Spectra_Tpc(i)),' s'],['Tc: ',num2str(Spectra_Tc(i)),' s']};

            if plotrun == 2
                % Annotate Images with correct answers for lecturer
                annotation('textbox',dim,'String',str,'FitBoxToText','off');
            end
            if plotrun ==1
                saveas(gcf,['Spectrum ',num2str(i),'_Student.png'])
            else
                saveas(gcf,['Spectrum ',num2str(i),'.png'])
            end
        end
        
        %Output the Spectral files for students

        SpectralOrdinatesOutput{i} = SpectralOrdinates;
        SpectralFrequencyOutput{i} = Spectra_freq;
        SpectralOutputMatrix = [SpectralFrequencyOutput{i},SpectralOrdinatesOutput{i}]
        SpectralOutputMatrixSave{i} = SpectralOutputMatrix;
        csvfilename = ['Spectrum',num2str(i),'.csv']
        csvwrite(csvfilename,SpectralOutputMatrix)

        %Output the Spectral Summary file for Lecturer

        DataToWrite=SpectralOutputMatrix';
        myfile = fopen(SaveFileName ,'a');
        fprintf(myfile,['Spectrum Number: ',num2str(i),'\n']);
        fprintf(myfile,[SpectraName]);
        fprintf(myfile,',\n');
        fprintf(myfile,['Infagravity Wave: ',num2str(InfoGrav{i})]);
        fprintf(myfile,',\n');
        fprintf(myfile,['Hm0: ',num2str(Spectra_Hm0(i)),' m']);
        fprintf(myfile,',\n');
        fprintf(myfile,['T02: ',num2str(Spectra_T02(i)),' s']);
        fprintf(myfile,',\n');
        fprintf(myfile,['Te: ',num2str(Spectra_Te(i)),' s']);
        fprintf(myfile,',\n');
        fprintf(myfile,['Tp (Calculated): ',num2str(Spectra_Tpc(i)),' s']);
        fprintf(myfile,',\n');
        fprintf(myfile,['Tc: ',num2str(Spectra_Tc(i)),' s']);
        fprintf(myfile,',\n');
        dlmwrite(SaveFileName,DataToWrite, '-append','roffset',2,'delimiter',',');
        fprintf(myfile,',\n');
        fclose(myfile);
        %close all
    end
    close all
end

