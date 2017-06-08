% Alley Robinson 
% April 26, 2017 - May 7, 2017
% June 1, 2017 - June 8, 2017 - enhancements and bug fixes
% LASP Interview Assignment
% This program will be able to calculate irradiance measurement (W/m^2) and graph all results

format long g; % For accuracy purposes

% INPUTS: read the CSV files provided
referenceSpectrumData = csvread('referenceSpectrum.txt',1,0);
wavelength_ref   = referenceSpectrumData(:,1);
irradiance_ref   = referenceSpectrumData(:,2);

detectorTempsData = csvread('detectorTemp.txt',1,0);
microsecondsSinceGpsEpoch_DT = detectorTempsData(:,1);
tempInC                      = detectorTempsData(:,2);
   
distanceAndDopplerData = csvread('distanceAndDoppler.txt',1,0);
microsecondsSinceGpsEpoch_DD   = distanceAndDopplerData(:,1);
sunObserverDistanceCorrection  = distanceAndDopplerData(:,2);
sunObserverDopplerFactor       = distanceAndDopplerData(:,3);
  
instrumentTelemetryData = csvread('instrumentTelemetry.txt',1,0);
microsecondsSinceGpsEpoch_TELE = instrumentTelemetryData(:,1);
gratPos                        = instrumentTelemetryData(:,2);
counts                         = instrumentTelemetryData(:,3);
  
integrationTimeData = csvread('integrationTime.txt',1,0);
microsecondsSinceGpsEpoch_INT = integrationTimeData(:,1);
intTime                       = integrationTimeData(:,2);

% csvread will not work due to the planName column being a string
plansData = readtable('plans.txt');
planName  = plansData(:,1);
startTime = plansData(:,2);
endTime   = plansData(:,3);


% PROCESSING: calculations split into individual functions
    measured_WL = wavelengthCalculator(gratPos);
    
    [cr, photonsPerSecondPerCm2] = eventsCalculator(microsecondsSinceGpsEpoch_INT,...
        microsecondsSinceGpsEpoch_TELE, intTime, counts);
          
    irradiance = irradianceCalculator(measured_WL, photonsPerSecondPerCm2);
    
    % Yi = interpt1(X, Y, Xi, 'method')
    final_irradiance_ref = interp1(wavelength_ref, irradiance_ref, measured_WL, 'linear');
    

% Irradiance at the different scans using logical indexing

    % QuickScan
    quickScan = wavelengthCalculator(gratPos(165:17932));
        
    [qs_cr, qs_photonsPerSecPerCm2] = eventsCalculator(microsecondsSinceGpsEpoch_INT(165:4607),...
        microsecondsSinceGpsEpoch_TELE(165:17932), intTime(165:4607), counts(165:17932));
        
    qs_irradiance = irradianceCalculator(quickScan, qs_photonsPerSecPerCm2);
        
    % Constant Wavelength
    constant_WL = wavelengthCalculator(gratPos(19328:28193));
        
    [cw_cr, cw_photonsPerSecPerCm2] = eventsCalculator(microsecondsSinceGpsEpoch_INT(6003:10436),...
        microsecondsSinceGpsEpoch_TELE(19328:28193), intTime(6003:10436), counts(19328:28193));
         
    cw_irradiance = irradianceCalculator(constant_WL, cw_photonsPerSecPerCm2);
        
    % Down Scan
    downScan = wavelengthCalculator(gratPos(29598:32126));
        
    [ds_cr, ds_photonsPerSecPerCm2] = eventsCalculator(microsecondsSinceGpsEpoch_INT(11841:16264),...
        microsecondsSinceGpsEpoch_TELE(29598:32126), intTime(11841:16264), counts(29598:32126));
        
    ds_irradiance = irradianceCalculator(downScan, ds_photonsPerSecPerCm2);
        
    % Dark
    dark = wavelengthCalculator(gratPos(33540:37954));
        
    [drk_cr, drk_photonsPerSecPerCm2] = eventsCalculator(microsecondsSinceGpsEpoch_INT(17678:22092),...
        microsecondsSinceGpsEpoch_TELE(33540:37954), intTime(17678:22092), counts(33540:37954));
        
    drk_irradiance = irradianceCalculator(dark, drk_photonsPerSecPerCm2);
        
    % Up Scan
    upScan = wavelengthCalculator(gratPos(39377:41892));
        
    [us_cr, us_photonsPerSecPerCm2] = eventsCalculator(microsecondsSinceGpsEpoch_INT(23515:27918),...
        microsecondsSinceGpsEpoch_TELE(39377:41892), intTime(23515:27918), counts(39377:41892));
        
    us_irradiance = irradianceCalculator(upScan, us_photonsPerSecPerCm2);

% OUTPUTS: plotting reference data and measured data
    figure('Name', 'Irradiance Data', 'NumberTitle', 'off');
        for i = 1:size(measured_WL)
          if measured_WL(i) > 185
              measured_WL(i) = NaN;
          end
        end   
    plotReference(wavelength_ref, irradiance_ref);
    title('Irradiance Data');
    hold on;
    plot(measured_WL, final_irradiance_ref, 'm', 'DisplayName', 'Measured Data');
    hold off;
    legend show;
    
    % Plotting reference data and measured data at the different scans
    % Quick Scan
    figure('Name', 'Quick Scan', 'NumberTitle', 'off');
        for i = 1:size(quickScan)
          if quickScan(i) > 185
              quickScan(i) = NaN;
          end
        end       
    plotReference(wavelength_ref, irradiance_ref);
    title('Irradiance at the Quick Scan');
    hold on;
    plot(quickScan, qs_irradiance, 'g', 'DisplayName', 'Measured Data');
    hold off;
    legend show;
    
    % Constant Wavelength
    figure('Name', 'Constant Wavelength', 'NumberTitle', 'off');
        for i = 1:size(constant_WL)
          if constant_WL(i) > 185
              constant_WL(i) = NaN;
          end
        end   
    plotReference(wavelength_ref, irradiance_ref);
    title('Irradiance at a Constant Wavelength');
    hold on;
    plot(constant_WL, cw_irradiance, 'b', 'DisplayName', 'Measured Data');
    hold off;
    legend show;
    
    % Down Scan
    figure('Name', 'Down Scan', 'NumberTitle', 'off');
        for i = 1:size(downScan)
          if downScan(i) > 185
               downScan(i) = NaN;
          end
        end
    plotReference(wavelength_ref, irradiance_ref);
    title('Irradiance at the Down Scan');
    hold on;
    plot(downScan, ds_irradiance, 'c', 'DisplayName', 'Measured Data');
    hold off;
    legend show;
    
    % Dark Scan
    figure('Name', 'Dark', 'NumberTitle', 'off');
    plot(dark, drk_irradiance, 'c', 'DisplayName', 'Measured Data');
    title('Irradiance at the Dark Scan');
    xlabel('Wavelength (nm)');
    ylabel('Solar Irradiance (watts/m^2)');
    legend show;
    
    % Up Scan
    figure('Name', 'Up Scan', 'NumberTitle', 'off');
        for i = 1:size(upScan)
          if upScan(i) > 185
               upScan(i) = NaN;
          end
        end
    plotReference(wavelength_ref, irradiance_ref);
    title('Irradiance at the Up Scan');
    hold on;
    plot(upScan, us_irradiance, 'r', 'DisplayName', 'Measured Data');
    hold off;
    legend show;

% PROCESSING:
function measured_WL = wavelengthCalculator(gratPos)
       % wavelengthCalculator: calculates the wavelength in nm from the grating equation
       
     offset     = 239532.38;
     stepSize   = 2.4237772022101214e-6;
     phiGInRads = 0.08503244115716374;
     d          = 277.77777777777777;
        
     ang1        = (offset - gratPos) * stepSize;
     measured_WL = 2 * d * sin(ang1) * cos(phiGInRads / 2.0); %[nm]  
end

function [cr, photonsPerSecondPerCm2] = eventsCalculator(microsecondsSinceGpsEpoch_INT,...
        microsecondsSinceGpsEpoch_TELE, intTime, counts)
    % eventsCalculator: calculates the count rate and photon events 
      
    apArea           = 0.01; %[cm^2]
    conversionFactor = 1000;
     
    intTimeConversion = intTime ./ conversionFactor; 
       
    integrationTime = interp1(microsecondsSinceGpsEpoch_INT, intTimeConversion,...
        microsecondsSinceGpsEpoch_TELE, 'nearest'); % Interpolation of time
      
    cr                     = counts ./ integrationTime; %[counts / sec]
    photonsPerSecondPerCm2 = cr ./ apArea; %[photons/sec/cm^2] 
end


function irradiance = irradianceCalculator(measured_WL, photonsPerSecondPerCm2)
    % irradianceCalculator: calculates solar irradiance

    conversionFactor = 0.0000000010; 
    h                = 6.62606957e-34; %[m^2 * kg / s]
    c                = 299792458.0; %[m/s]
      
    wl_meters = measured_WL * conversionFactor;
    energyPerPhoton = h * c ./ wl_meters;
       
    irradiance = photonsPerSecondPerCm2 .* 100 .* 100 .* energyPerPhoton; %[watts/m^2]
end

function plotReference(wavelength_ref, irradiance_ref) 
     % plotReference: plots the reference data
         
     plot(wavelength_ref, irradiance_ref, 'k', 'Linewidth', 2.00, 'DisplayName','Reference Data');
     xlim([179 181]);
     set(gca,'xtick', 179:0.1:181);
     xlabel('Wavelength (nm)');
     ylabel('Solar Irradiance (watts/m^2)');
end  