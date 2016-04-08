%Generate Spectra Via Bretschneider, Pierson Moskowitz etc.

function [BS,BS_freq] = generateBretschneiderSpectrum(H,T,freq)
B = (0.751/T)^4;
A = B*(H^2)/4;
%BS_freq=0.005:0.001:.516;
BS_freq = freq';
BS = (A.*BS_freq(:,1).^(-5)).*exp(-B.*(BS_freq(:,1).^(-4)));
BS_df=BS_freq(3,1)-BS_freq(2,1);
end
