% Carter Lybbert 2022. University of Utah. carterlybbert46@gmail.com

%% Use this script to generate effect site concentration curves for each treatment

clc;
%get participant number
Participant_Demographics = [0, 0, 0, 0];
participant_number = input('Enter the participant number using two digits: ');
participant_number_str = num2str(participant_number);
Participant_Demographics_Filename = sprintf('Participant_Demographics_PROP%d.mat',participant_number);

try %define where everything will be saved
    save_directory = sprintf('C:\\Users\\u6027582\\OneDrive\\Masters Materials\\Matlab\\PROP%d_Matlab\\', participant_number);
    cd(save_directory)
    addpath( 'C:\Users\u6027582\OneDrive\Masters Materials\Matlab\');
    addpath(save_directory);
    
catch
    save_directory = sprintf('C:\\OneDrive\\Masters Materials\\Matlab\\PROP%d_Matlab\\', participant_number);
    cd(save_directory)
    addpath( 'C:\OneDrive\Masters Materials\Matlab\');
    addpath(save_directory);
end

try % get participant demographics
    Participant_Demographics = csvread(['PROP',participant_number_str,'_Demographics.csv']);
catch
    %if Participant_Demographics == [0 0 0 0]
    for n = 1:4
        if n==1
            Participant_Demographics(n) = input('Enter the sex of particpant as a number (1 being male and 0 being female): ');
        elseif n==2
            Participant_Demographics(n) = input('Enter the height of particpant in cm to one decimal place: ');
        elseif n==3
            Participant_Demographics(n) = input('Enter the age of the particpant in years from conception: ');
        elseif n==4
            Participant_Demographics(n) = input('Enter the weight of the participant in kg to one decimal point: ');
        end
    end
    csvwrite(['PROP',participant_number_str,'_Demographics.csv'], Participant_Demographics); % Save the modeled parameters as a CSV in an Excel Spreadsheet
end

input('\nOnly continue once you have entered the .spa file of the BSR into an excel spreadsheet \n Press enter once this is completed'); %see SOP document for how to do this
number_of_treatments = input('Enter the treatment numbers of the past treatments that you want to include in this analysis as a vector\n in the following format [1 2 3 5] : ');
PROP_Ke0 = []; % create blank vectors to store the PK/PD parameters
PROP_Hill = [];
PROP_EC50 = [];
PROP_RMSE = [];
step = 1; % a stepping variable

% define fit function parameters... look up help fit.
fitType = @(Hill, EC50, x) 100./(1+(EC50./x).^Hill); % This defines the equation used by the fit function, the equation of our regression. This is the Hill equation or sigmoid Emax model
fitOptions = fitoptions( 'Method', 'NonlinearLeastSquares', ...
    'Lower', [1 1], ... % lower guess for hill and Ec50
    'Upper', [25 25], ... %upper guess for hill and Ec50
    'Start', [7 6], ... %starting point. 7 and 6 are a decent starting point.
    'Robust', 'LAR');

for treatment_selector = 1:length(number_of_treatments) % loop for each treatment
    Infusion_dose = csvread(['PROP',participant_number_str,'_T',num2str(number_of_treatments(treatment_selector)),'_Dosing.csv']); % second by second dosing in mcg
    BSR = csvread(['PROP',participant_number_str,'_T',num2str(number_of_treatments(treatment_selector)),'_BSR.csv']); % Load second-by-second BSR of the treatment associated with loaded dosing.
    length_BSR = length(BSR); % find the time duration of the BSR vector
    BSR = BSR(1:step:end); % I'm not sure that this does anything. Maybe delete.
    ke0_test_vector = [0.00:0.05:3]; % start with 0.08 and work your way up to 2, testing each value with the model
    RMSEVector = [];
    PDVector = [];
    y_error_saver = [];
    
    
    for m = 1:length(ke0_test_vector) %guessing and checking various Ke0 values to find the one associated with the lowest Root Mean Square Error. m is the index of the selected Ke0
        ke0 = ke0_test_vector(m); %start with lowest index of first value in Ke0 test vector.
        [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, ke0); % Use eleveld model function to find plasma and effect site concentration of propofol
        effect_conc = Ce_effect_concentration(1:step:length_BSR)'; % shorten Ce_effect concentration to the length of BSR and transpose
        [fitobject, gof, output] = fit(effect_conc, BSR, fitType, fitOptions); % fit type and fit object are defined above. This performs a nonlinear leastsquares regression of BSR to effect site concentration using the Hill equation a the regression equation
        PD = coeffvalues(fitobject); %This saves the optimal Hill and EC50 values from the regression that produced the best fit
        RMSE = gof.rmse; %safe RMSE of the regression.
        RMSEVector = [RMSEVector RMSE];
        PDVector = [PDVector; PD]; %store the PD factors that produced the lowest values.
        y_error = mean((output.residuals).^2);
        y_error_saver = [y_error_saver y_error];
    end
    
    [min_RMSE, min_RMSE_index] = min(RMSEVector); %find the index of the lowest RMSE
    Oke0 = ke0_test_vector(min_RMSE_index); % lowest RMSE in Ke0 test vector.
    ke0_test_vector = [Oke0-0.05:0.025:Oke0+0.05]; % zoom in and perform another guess and check where you try to find lowest RMSE of various Ke0 values.
    RMSEVector = [];
    PDVector = []; %do it again but increase the resolution
    y_error_saver = [];
    
    
    for m = 1:length(ke0_test_vector)
        ke0 = ke0_test_vector(m);
        [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, ke0);
        effect_conc = Ce_effect_concentration(1:step:length_BSR)';
        [fitobject, gof, output] = fit(effect_conc, BSR, fitType, fitOptions);
        PD = coeffvalues(fitobject);
        RMSE = gof.rmse;
        RMSEVector = [RMSEVector RMSE];
        PDVector = [PDVector; PD];
        y_error = mean((output.residuals).^2);
        y_error_saver = [y_error_saver y_error];
    end
    
    [min_RMSE, min_RMSE_index] = min(RMSEVector);
    Oke0 = ke0_test_vector(min_RMSE_index);
    ke0_test_vector = [Oke0-0.005:0.001:Oke0+0.005];
    RMSEVector = [];
    PDVector = [];
    y_error_saver = [];
    
    for m = 1:length(ke0_test_vector) % increaseing the resolution again to find lowest RMSE Ke0
        ke0 = ke0_test_vector(m);
        [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, ke0);
        effect_conc = Ce_effect_concentration(1:step:length_BSR)';
        [fitobject, gof, output] = fit(effect_conc, BSR, fitType, fitOptions);
        PD = coeffvalues(fitobject);
        RMSE = gof.rmse;
        RMSEVector = [RMSEVector RMSE];
        PDVector = [PDVector; PD];
        y_error = mean((output.residuals).^2);
        y_error_saver = [y_error_saver y_error];
    end
    
    
    [min_RMSE, min_RMSE_index] = min(y_error_saver); %find the best vector with the least error and find the associated Hill, EC50 and Ke0.
    Oke0 = ke0_test_vector(min_RMSE_index);
    
    % [min_RMSE, min_RMSE_index] = min(RMSEVector); %find the best vector with the least error and find the associated Hill, EC50 and Ke0.
    % Oke0 = ke0_test_vector(min_RMSE_index);
    PROP_Ke0 = [PROP_Ke0 Oke0];
    PROP_Hill = [PROP_Hill PDVector(min_RMSE_index,1)];
    PROP_EC50 = [PROP_EC50 PDVector(min_RMSE_index, 2)];
    PROP_RMSE = [PROP_RMSE RMSEVector(min_RMSE_index)];
    Pharm_parameters_saver = [PROP_Ke0; PROP_Hill; PROP_EC50; PROP_RMSE]';
    [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, Oke0);
    cd(save_directory)
    csvwrite(['PROP',participant_number_str,'_T',num2str(number_of_treatments(treatment_selector)),'_Ce.csv'], Ce_effect_concentration); % Save effect site concentration as a CSV file.
end

% Plot the BSR as a function of effect site concentration.
% plot(fitobject, x, BSR);
% xlabel('Predicted Ce (mcg/mL)');
% ylabel('BSR');
% legend('Observed BSR', 'RMSE Fit', 'Location', "northwest")
% title('PROP PD Model (RMSE Optimization)')
% Pharm_parameters_saver = [PROP_Ke0; PROP_Hill; PROP_EC50; PROP_RMSE]'
%display pharm parameters saver is a more user friendly format
