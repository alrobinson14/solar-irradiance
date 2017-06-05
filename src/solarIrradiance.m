% Alley Robinson 
% April 26, 2017 - May 7, 2017
% LASP Interview Assignment
% This program will be able to calculate irradiance measurement (W/m^2), plot the data, and compare to
% the reference spectrum provided by LASP. 

format long g; % for accuracy purposes

% INPUTS: read the CSV files provided
% Read in the file, and split the columns into its own vector

  
referenceSpectrumData = csvread('referenceSpectrum.txt',1,0);
wavelengthRef   = referenceSpectrumData(:,1);
irradianceRef   = referenceSpectrumData(:,2);

detectorTempsData = csvread('detectorTemp.txt',1,0);
microsecondsSinceGpsEpochDT = detectorTempsData(:,1);
tempInC                     = detectorTempsData(:,2);
   
distanceAndDopplerData = csvread('distanceAndDoppler.txt',1,0);
microsecondsSinceGpsEpochDD   = distanceAndDopplerData(:,1);
sunObserverDistanceCorrection = distanceAndDopplerData(:,2);
sunObserverDopplerFactor      = distanceAndDopplerData(:,3);
  
instrumentTelemetryData = csvread('instrumentTelemetry.txt',1,0);
microsecondsSinceGpsEpochTELE = instrumentTelemetryData(:,1);
gratPos                       = instrumentTelemetryData(:,2);
counts                        = instrumentTelemetryData(:,3);
  
integrationTimeData = csvread('integrationTime.txt',1,0);
microsecondsSinceGpsEpochINT = integrationTimeData(:,1);
intTime                      = integrationTimeData(:,2);

% Plans.txt has strings in the first column, so csvread will not work here
plansData = readtable('plans.txt');
planName  = plansData(:,1);
startTime = plansData(:,2);
endTime   = plansData(:,3);


% PROCESSING: calculations split into individual functions
measured_WL = wavelengthCalculator(gratPos);

[cr, photonsPerSecondPerCm2] = eventsCalculator(microsecondsSinceGpsEpochINT,...
        microsecondsSinceGpsEpochTELE, intTime, counts);
      

irradiance = irradianceCalculator(measured_WL, photonsPerSecondPerCm2);


% % QuickScan
% quickScan = wavelengthCalculator(gratPos(165:17932));
% 
% [qsIntTime, qsCr, qsPhotonsPerSecPerCm2] = eventsCalculator(microsecondsSinceGpsEpochINT(165:4607),...
%     microsecondsSinceGpsEpochTELE(165:17932), intTime(165:4607), counts(165:17932));
% 
% qsIrr = irradianceCalculator(quickScan, qsPhotonsPerSecPerCm2);
% 
% % Constant Wavelength
% constantWavelength = wavelengthCalculator(gratPos(19328:28193));
% 
% [cwIntTime, cwCr, cwPhotonsPerSecPerCm2] = eventsCalculator(microsecondsSinceGpsEpochINT(6003:10436),...
%      microsecondsSinceGpsEpochTELE(19328:28193), intTime(6003:10436), counts(19328:28193));
%  
% cwIrr = irradianceCalculator(constantWavelength, cwPhotonsPerSecPerCm2);
% 
% % Down Scan
% downScan = wavelengthCalculator(gratPos(29598:32126));
% 
% [dsIntTime, dsCr, dsPhotonsPerSecPerCm2] = eventsCalculator(microsecondsSinceGpsEpochINT(11841:16264),...
%      microsecondsSinceGpsEpochTELE(29598:32126), intTime(11841:16264), counts(29598:32126));
% 
% dsIrr = irradianceCalculator(downScan, dsPhotonsPerSecPerCm2);
% 
% % Dark
% dark = wavelengthCalculator(gratPos(33540:37954));
% 
% [drkIntTime, drkCr, drkPhotonsPerSecPerCm2] = eventsCalculator(microsecondsSinceGpsEpochINT(17678:22092),...
%      microsecondsSinceGpsEpochTELE(33540:37954), intTime(17678:22092), counts(33540:37954));
% 
% drkIrr = irradianceCalculator(dark, drkPhotonsPerSecPerCm2);
% 
% 
% % Up Scan
% upScan = wavelengthCalculator(gratPos(39377:41892));
% 
% [usIntTime, usCr, usPhotonsPerSecPerCm2] = eventsCalculator(microsecondsSinceGpsEpochINT(23515:27918),...
%      microsecondsSinceGpsEpochTELE(39377:41892), intTime(23515:27918), counts(39377:41892));
% 
% usIrr = irradianceCalculator(upScan, usPhotonsPerSecPerCm2);
% 
% % OUTPUTS: plotting reference data and measured data
 figure('Name', 'Irradiance Data', 'NumberTitle', 'off');
% % Estimated fix for the outliers
%  s = size(measured_WL);
%  for i = 1:s
%    if measured_WL(i) > 185
%         measured_WL(i) = NaN;
%    end
%  end    
%  plot(wavelengthRef, irradianceRef, 'k', 'Linewidth', 1.00, 'DisplayName','Reference Data');
%  xlim([178 182]);
%  title('Irradiance Data');
%  xlabel('Wavelength (nm)');
%  ylabel('Solar Irradiance (watts/m^2)');
%  hold on;
%  plot(measured_WL, irradiance, 'm', 'DisplayName', 'Measured Data');
%  hold off;
%  legend show;
% %  
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
% 

% PROCESSING - process the data, store that in separate values

% This function will calculate the ang1 values and the wavelength in nm from the grating equation
function measured_WL = wavelengthCalculator(gratPos)
    
    % Initialize variables to what has been given
    offset     = 239532.38;
    stepSize   = 2.4237772022101214e-6;
    phiGInRads = 0.08503244115716374;
    d          = 277.77777777777777;
    
    ang1 = (offset - gratPos) * stepSize;
    measured_WL = 2 * d * sin(ang1) * cos(phiGInRads / 2.0); %[nm]  
end

% This function will calculate the counds/seconds/area and photon events by using the provided equation
function [cr, photonsPerSecondPerCm2] = eventsCalculator(microsecondsSinceGpsEpochINT,...
        microsecondsSinceGpsEpochTELE, intTime, counts)
        
     % Initializing variables provided 
       apArea = 0.01; %[cm^2]
       conversionFactor = 1000;
     
       intTimeConversion = intTime ./ conversionFactor; 
       
       % Now, we need to use interpolation to get the data to align
       integrationTime = interp1(microsecondsSinceGpsEpochINT, intTimeConversion,...
           microsecondsSinceGpsEpochTELE, 'nearest');
      
       cr = counts ./ integrationTime; %[counts / sec]
      
       photonsPerSecondPerCm2 = cr ./ apArea; %[photons/sec/cm^2] 
end

% This function will calculate solar irradiance
function irradiance = irradianceCalculator(measured_WL, photonsPerSecondPerCm2)
      
      % Initializing variables provided
       conversionFactor = 0.0000000010; 
       h = 6.62606957e-34; %[m^2 * kg / s]
       c = 299792458.0; %[m/s]
      
       wl_meters = measured_WL * conversionFactor;
       energyPerPhoton = h * c ./ wl_meters;
       
      irradiance = photonsPerSecondPerCm2 .* 100 .* 100 .* energyPerPhoton; %[watts/m^2]
end
