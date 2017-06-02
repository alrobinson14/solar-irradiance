% Alley Robinson 
% April 26, 2017 - May 7, 2017
% LASP Interview Assignment
% This program will be able to calculate irradiance measurement (W/m^2), plot the data, and compare to
% the reference spectrum provided by LASP. 

format long g; 

% INPUTS: read the CSV files provided
  
referenceSpectrum = csvread('referenceSpectrum.txt',1,0);
wavelengthRef   = referenceSpectrum(:,1);
irradianceRef   = referenceSpectrum(:,2);

detectorTemps = csvread('detectorTemp.txt',1,0);
microsecondsSinceGpsEpochDT = detectorTemps(:,1);
tempInC                     = detectorTemps(:,2);
   
distanceAndDoppler = csvread('distanceAndDoppler.txt',1,0);
microsecondsSinceGpsEpochDD   = distanceAndDoppler(:,1);
sunObserverDistanceCorrection = distanceAndDoppler(:,2);
sunObserverDopplerFactor      = distanceAndDoppler(:,3);
  
instrumentTelemetry = csvread('instrumentTelemetry.txt',1,0);
microsecondsSinceGpsEpochTELE = instrumentTelemetry(:,1);
gratPos                       = instrumentTelemetry(:,2);
counts                        = instrumentTelemetry(:,3);
  
integrationTimetxt = csvread('integrationTime.txt',1,0);
microsecondsSinceGpsEpochINT = integrationTimetxt(:,1);
intTime                      = integrationTimetxt(:,2);

plans = readtable('plans.txt');
planName  = plans(:,1);
startTime = plans(:,2);
endTime   = plans(:,3);


% Step 1: Get the Wavelength from the Grating equation
wavelength = wavelengthCalculator(gratPos);


% Step 2: Find out the counts/second/area 
% Step 2(a): Interpolate the integration time in order to successfully get a count rate
integrationTime = interpolateTime(microsecondsSinceGpsEpochINT,...
           microsecondsSinceGpsEpochTELE, intTime);
           
% Step 2(b): Calculate the count rate and the photonsPerSecondPerCm2
[cr, photonsPerSecondPerCm2] = eventsCalculator(counts, integrationTime);
     

% Step 3: Calculate the irradiance
irradiance = irradianceCalculator(wavelength, photonsPerSecondPerCm2);
% 
% % QuickScan
% quickScan = wavelengthCalculator(gratPos(165:17932));
% 
% [qsCr, qsPhotonsPerSecPerCm2] = eventsCalculator(counts(165:17932), integrationTime(165:4607));
% 
% qsIrr = irradianceCalculator(quickScan, qsPhotonsPerSecPerCm2);
% 
% % Constant Wavelength
% constantWavelength = wavelengthCalculator(gratPos(19328:28193));
% 
% [cwCr, cwPhotonsPerSecPerCm2] = eventsCalculator(counts(19328:28193), integrationTime(6003:10436));
%  
% cwIrr = irradianceCalculator(constantWavelength, cwPhotonsPerSecPerCm2);
% 
% % Down Scan
% downScan = wavelengthCalculator(gratPos(29598:32126));
% 
% [dsCr, dsPhotonsPerSecPerCm2] = eventsCalculator(counts(29598:32126), intTime(11841:16264));
% 
% dsIrr = irradianceCalculator(downScan, dsPhotonsPerSecPerCm2);
% 
% % Dark
% dark = wavelengthCalculator(gratPos(33540:37954));
% 
% [drkIntTime, drkCr, drkPhotonsPerSecPerCm2] = eventsCalculator(counts(33540:37954), intTime(17678:22092));
% 
% drkIrr = irradianceCalculator(dark, drkPhotonsPerSecPerCm2);
% 
% 
% % Up Scan
% upScan = wavelengthCalculator(gratPos(39377:41892));
% 
% [usCr, usPhotonsPerSecPerCm2] = eventsCalculator(counts(39377:41892), intTime(23515:27918));
% 
% usIrr = irradianceCalculator(upScan, usPhotonsPerSecPerCm2);
% 
% 
[wavelength, index] = unique(wavelength);
wlFinal = interp1(wavelength, irradiance(index), wavelengthRef, 'linear');

% OUTPUTS: plotting reference data and measured data
% figure('Name', 'Irradiance Data', 'NumberTitle', 'off');
% % Estimated fix for the outliers
% s = size(wavelength);
% for i = 1:s
%   if wavelength(i) > 185
%        wavelength(i) = NaN;
%   end
% end    
% plot(wavelengthRef, irradianceRef, 'r', 'DisplayName','Reference Data');
% xlim([178 182]);
% title('Irradiance Data');
% xlabel('Wavelength (nm)');
% ylabel('Solar Irradiance (watts/m^2)');
% hold on;
% plot(wavelength, irradiance, 'm', 'DisplayName', 'Measured Data');
% hold off;
% legend show;
% 
% 
% figure('Name', 'Quick Scan', 'NumberTitle', 'off');
%  
% 
% s = size(quickScan);
% for i = 1:s
%   if quickScan(i) > 185
%        quickScan(i) = NaN;
%   end
% end
%   
% 
% plot(quickScan, qsIrr, 'g');
% title('Irradiance at the Quick Scan');
% xlabel('Wavelength (nm)');
% ylabel('Solar Irradiance (watts/m^2)');
% 
% figure('Name', 'Constant Wavelength', 'NumberTitle', 'off');
% 
% s = size(constantWavelength);
% for i = 1:s
%   if constantWavelength(i) > 185
%        constantWavelength(i) = NaN;
%   end
% end
% 
% plot(constantWavelength, cwIrr, 'b');
% title('Irradiance at a Constant Wavelength');
% xlabel('Wavelength (nm)');
% ylabel('Solar Irradiance (watts/m^2)');
% 
% figure('Name', 'Down Scan', 'NumberTitle', 'off');
% 
% s = size(downScan);
% for i = 1:s
%   if downScan(i) > 185
%        downScan(i) = NaN;
%   end
% end
% 
% plot(downScan, dsIrr, 'y');
% title('Irradiance at the Down Scan');
% xlabel('Wavelength (nm)');
% ylabel('Solar Irradiance (watts/m^2)');
% 
% figure('Name', 'Dark', 'NumberTitle', 'off');
% 
% plot(dark, drkIrr, 'k');
% title('Irradiance at the Dark Scan');
% xlabel('Wavelength (nm)');
% ylabel('Solar Irradiance (watts/m^2)');
% 
% figure('Name', 'Up Scan', 'NumberTitle', 'off');
% 
% s = size(upScan);
% for i = 1:s
%   if upScan(i) > 185
%        upScan(i) = NaN;
%   end
% end
% 
% plot(upScan, usIrr, 'c');
% title('Irradiance at the Up Scan');
% xlabel('Wavelength (nm)');
% ylabel('Solar Irradiance (watts/m^2)');


% PROCESSING - process the data, store that in separate values
function integrationTime = interpolateTime(microsecondsSinceGpsEpochINT,...
           microsecondsSinceGpsEpochTELE, intTime)
      
      conversionFactor = 1000;   
       
      intTimeConversion = conversionFactor ./ intTime;     
      integrationTime = interp1(microsecondsSinceGpsEpochINT, intTimeConversion,...
           microsecondsSinceGpsEpochTELE, 'nearest');
       
end              

% This function will calculate the ang1 values and the wavelength in nm from the grating equation
function wavelength = wavelengthCalculator(gratPos)
    
     % Initialize variables to what has been given
     offset     = 239532.38;
     stepSize   = 2.4237772022101214e-6;
     phiGInRads = 0.08503244115716374;
     d          = 277.77777777777777;
    
     ang1 = (offset - gratPos) * stepSize;
     wavelength = 2 * d * sin(ang1) * cos(phiGInRads / 2.0); %[nm]  
    
end

% This function will calculate the counds/seconds/area and photon events by using the provided equation
function [cr, photonsPerSecondPerCm2] = eventsCalculator(counts, integrationTime)
        
       apArea = 0.01; 
    
       cr = counts ./ integrationTime;        %[counts / sec]
       photonsPerSecondPerCm2 = cr ./ apArea; %[photons/sec/cm^2] 
       
end

% This function will calculate solar irradiance
function irradiance = irradianceCalculator(wavelength, photonsPerSecondPerCm2)
      
      % Initializing variables provided
       conversionFactor = 0.0000000010; 
       h = 6.62606957e-34; %[m^2 * kg / s]
       c = 299792458.0; %[m/s]
      
       wavelengthInM = wavelength * conversionFactor;
       energyPerPhoton = h * c ./ wavelengthInM;
       
      irradiance = photonsPerSecondPerCm2 .* 100 .* 100 .* energyPerPhoton; %[watts/m^2]
end

