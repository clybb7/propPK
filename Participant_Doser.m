% Create Vector for Total Drug Admin(t). Adapted from createInfusion.m by Jason Huang
% This is a tool to help you transcribe the paper document tracking dosing during the procedure and convert it into an excel .csv file for further analysis
% Written by Carter Lybbert, 2022. carterlybbert46@gmail.com 

%% Select the correct save directory and load correct participant information
clc; clear
Participant_Demographics = [0, 0, 0, 0];
participant_number = input('Enter the participant number: ');
participant_number_str = num2str(participant_number);
Participant_Demographics_Print = sprintf('Participant_Demographics_PROP%d.mat',participant_number);
first_treatment = input('Is this the first high dose propofol treatment that this patient has undergone? Enter 1 for yes, 0 for no: ');
save_directory = sprintf('C:\\Users\\u6027582\\OneDrive\\Masters Materials\\Matlab\\');
cd(save_directory); addpath( 'C:\Users\u6027582\OneDrive\Masters Materials\Matlab\');
addpath(save_directory);

try
    save_directory = sprintf('C:\\Users\\u6027582\\OneDrive\\Masters Materials\\Matlab\\PROP%d_Matlab\\', participant_number);
    cd(save_directory)
    addpath( 'C:\Users\u6027582\OneDrive\Masters Materials\Matlab\');
    addpath(save_directory);
catch
    if exist(['PROP',participant_number_str,'_Matlab'],'dir') == 0
        mkdir(['PROP',participant_number_str,'_Matlab']);
    end
    save_directory = sprintf('C:\\Users\\u6027582\\OneDrive\\Masters Materials\\Matlab\\PROP%d_Matlab\\', participant_number);
    cd(save_directory)
    addpath( 'C:\Users\u6027582\OneDrive\Masters Materials\Matlab\');
    addpath(save_directory);
end

try
    Participant_Demographics = csvread(['PROP',participant_number_str,'_Demographics.csv']);
catch
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

if first_treatment == 0 %% Enter the dosing of past treatments if desired
    participant_entered = 1;
    enter_number_dosing_treatments = input('Enter the number of past treatment dosing schemes you would like to enter and save: ');
    
    if enter_number_dosing_treatments > 0
        for dosing_entry = 1:enter_number_dosing_treatments
            treatment_number = input('Enter the treatment number of this participant: ');
            treatment_number_str = num2str(treatment_number);
            Infusion_Saver = zeros(1,3600);
            pump_weight = input('Enter the weight of the participant that was entered into the infusion pump to one decimal place: ');
            number_of_boluses = input('Enter the number of boluses given during the treatment: ');
            number_diff_infusion_rates = input('Enter the number of different infusion rates given during the treatment: ');
            
            for bolus_selector = 1:number_of_boluses
                if bolus_selector == 1
                    Bolus_mg = input('Enter the mass in mg of the first bolus: ');
                    Bolus1_saver = num2str(Bolus_mg);
                    Bolus_duration = input('Enter the duration of the first bolus in seconds:');
                    Infusion_Saver(1:Bolus_duration) = Infusion_Saver(1:Bolus_duration) + (Bolus_mg*1000/Bolus_duration); %why are we multiplying bolus mass by 1000?
                    
                elseif bolus_selector == 2
                    Bolus2_mg = input('Enter the mass in mg of the second bolus: ');
                    Bolus2_duration = input('Enter the duration of the second bolus in seconds:');
                    Bolus2_time = (input('Enter the time in seconds after the beginning of treatment that the second bolus was administered: ')+1);
                    Infusion_Saver(Bolus2_time:(Bolus2_time+Bolus2_duration)) = Infusion_Saver(Bolus2_time:(Bolus2_time+Bolus2_duration)) + (Bolus2_mg*1000/Bolus2_duration);
                    
                elseif bolus_selector == 3
                    Bolus3_mg = input('Enter the mass in mg of the third bolus: ');
                    Bolus3_duration = input('Enter the duration of the third bolus in seconds: ');
                    Bolus3_time = (input('Enter the time in seconds after the beginning of treatment that the third bolus was administered: ')+1);
                    Infusion_Saver(Bolus3_time:(Bolus3_time+Bolus3_duration)) = Infusion_Saver(Bolus3_time:(Bolus3_time+Bolus3_duration)) + (Bolus3_mg*1000/Bolus3_duration);
                    
                elseif bolus_selector == 4
                    Bolus4_mg = input('Enter the mass in mg of the fourth bolus: ');
                    Bolus4_duration = input('Enter the duration of the fourth bolus in seconds: ');
                    Bolus4_time = (input('Enter the time in seconds after the beginning of treatment that the fourth bolus was administered: ')+1);
                    Infusion_Saver(Bolus4_time:(Bolus4_time+Bolus4_duration)) = Infusion_Saver(Bolus4_time:(Bolus4_time+Bolus4_duration)) + (Bolus4_mg*1000/Bolus4_duration);
                    
                elseif bolus_selector == 5
                    Bolus5_mg = input('Enter the mass in mg of the fifth bolus: ');
                    Bolus5_duration = input('Enter the duration of the fifth bolus in seconds: ');
                    Bolus5_time = (input('Enter the time in seconds after the beginning of treatment that the fifth bolus was administered: ')+1);
                    Infusion_Saver(Bolus5_time:(Bolus5_time+Bolus5_duration)) = Infusion_Saver(Bolus5_time:(Bolus5_time+Bolus5_duration)) + (Bolus5_mg*1000/Bolus5_duration);
                end
            end
           
            for infusion_selector = 1:number_diff_infusion_rates
                
                if infusion_selector ==1
                    infusion_rate = input('Enter the infusion rate of the first infusion in mcg/kg/min: ');
                    infusion_rate1_saver = num2str(infusion_rate);
                    infusion_start = input('Enter the time in seconds after the beginning of treatment that the first infusion was started: ');
                    infusion_end = input('Enter the end of the first infusion in seconds as recorded from the beginning of treatment: ');
                    Infusion_Saver(infusion_start:infusion_end) = Infusion_Saver(infusion_start:infusion_end) + infusion_rate*pump_weight/60;
                    
                elseif infusion_selector == 2
                    infusion_rate = input('Enter the infusion rate of the second infusion in mcg/kg/min: ');
                    infusion_start = infusion_end+1;
                    infusion_end = input('Enter the end of the second infusion in seconds from the beginning of treatment: ');
                    Infusion_Saver(infusion_start:infusion_end) = Infusion_Saver(infusion_start:infusion_end) + infusion_rate*pump_weight/60;
                    
                elseif infusion_selector == 3
                    infusion_rate = input('Enter the infusion rate of the third infusion in mcg/kg/min: ');
                    infusion_start = infusion_end+1;
                    infusion_end = input('Enter the end of the third infusion in seconds from the beginning of treatment: ');
                    Infusion_Saver(infusion_start:infusion_end) = Infusion_Saver(infusion_start:infusion_end) + infusion_rate*pump_weight/60;
                    
                elseif infusion_selector == 4
                    infusion_rate = input('Enter the infusion rate of the fourth infusion in mcg/kg/min: ');
                    infusion_start = infusion_end+1;
                    infusion_end = input('Enter the end of the fourth infusion in seconds from the beginning of treatment: ');
                    Infusion_Saver(infusion_start:infusion_end) = Infusion_Saver(infusion_start:infusion_end) + infusion_rate*pump_weight/60;
                    
                elseif infusion_selector == 5
                    infusion_rate = input('Enter the infusion rate of the fifth infusion in mcg/kg/min: ');
                    infusion_start = infusion_end+1;
                    infusion_end = input('Enter the end of the fifth infusion in seconds from the beginning of treatment: ');
                    Infusion_Saver(infusion_start:infusion_end) = Infusion_Saver(infusion_start:infusion_end) + infusion_rate*pump_weight/60;
                    
                elseif infusion_selector == 6
                    infusion_rate = input('Enter the infusion rate of the sixth infusion in mcg/kg/min: ');
                    infusion_start = infusion_end+1;
                    infusion_end = input('Enter the end of the sixth infusion in seconds from the beginning of treatment: ');
                    Infusion_Saver(infusion_start:infusion_end) = Infusion_Saver(infusion_start:infusion_end) + infusion_rate*pump_weight/60;
                    
                elseif infusion_selector == 7
                    infusion_rate = input('Enter the infusion rate of the seventh infusion in mcg/kg/min: ');
                    infusion_start = infusion_end+1;
                    infusion_end = input('Enter the end of the seventh infusion in seconds from the beginning of treatment: ');
                    Infusion_Saver(infusion_start:infusion_end) = Infusion_Saver(infusion_start:infusion_end) + infusion_rate*pump_weight/60;
                    
                elseif infusion_selector == 8
                    infusion_rate = input('Enter the infusion rate of the eighth infusion in mcg/kg/min: ');
                    infusion_start = infusion_end+1;
                    infusion_end = input('Enter the end of the eighth infusion in seconds from the beginning of treatment: ');
                    Infusion_Saver(infusion_start:infusion_end) = Infusion_Saver(infusion_start:infusion_end) + infusion_rate*pump_weight/60;
                    
                elseif infusion_selector == 9
                    infusion_rate = input('Enter the infusion rate of the ninth infusion in mcg/kg/min: ');
                    infusion_start = infusion_end+1;
                    infusion_end = input('Enter the end of the ninth infusion in seconds from the beginning of treatment: ');
                    Infusion_Saver(infusion_start:infusion_end) = Infusion_Saver(infusion_start:infusion_end) + infusion_rate*pump_weight/60;
                    
                elseif infusion_selector == 10
                    infusion_rate = input('Enter the infusion rate of the tenth infusion in mcg/kg/min: ');
                    infusion_start = infusion_end+1;
                    infusion_end = input('Enter the end of the tenth infusion in seconds from the beginning of treatment: ');
                    Infusion_Saver(infusion_start:infusion_end) = Infusion_Saver(infusion_start:infusion_end) + infusion_rate*pump_weight/60;
                end
            end
            clf; plot(Infusion_Saver) % 4.Preview of Combined Drug Administered during the treatment that you are saving as .csv
            hold on;
            xticks([0, 3*60:3*60:60*60])
            xticklabels({'0','3','6', '9', '12', '15', '18', '21', '24', '27', '30', '33', '36', '39','42','45','48','51','54','57','60'})
            xlim([0 30*60]); xlabel 'Time (Minutes)'; ylabel 'Infusion + Bolus (mcg/s)'
            title 'Total Drug Admin'
            csvwrite(['PROP',participant_number_str,'_T',treatment_number_str,'_Dosing.csv'], Infusion_Saver); %specify file name to export;
        end
    end
    
    % After this, import the .spa file into excel, then crop out the H column and save seperate BSR file as .csv
    % Then after you have both the Dosing and BSR data...proceed to determine ke0, Hill, EC50
%% Read in dosing from previous treatments if they exist
    input('\nOnly continue once you have entered the .spa file of the BSR into an excel spreadsheet \n Press enter once this is completed');
    number_of_treatments = [];
    try
        Pharm_parameters_saver = csvread(['PROP',participant_number_str,'_Compiled_PD_Model.csv']); %Load PK/PD parameters if previously determined
        use_previous_param = 1;
    catch
        use_previous_param = 0;
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(1),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [1];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(2),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [number_of_treatments 2];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(3),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [number_of_treatments 3];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(4),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [number_of_treatments 4];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(5),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [number_of_treatments 5];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(6),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [number_of_treatments 6];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(7),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [number_of_treatments 7];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(8),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [number_of_treatments 8];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(9),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [number_of_treatments 9];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(10),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [number_of_treatments 10];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(11),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [number_of_treatments 11];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(12),'_Dosing.csv']); % second by second dosing in mcg
        number_of_treatments = [number_of_treatments 12];
    catch
    end
    
    %% if PK/PD parameters have not been calculated for the past treatments of this participant, do that now and save them an a .CSV file
    if  use_previous_param == 0;
        PROP_Ke0 = []; % create blank vectors to store the PK/PD parameters
        PROP_Hill = [];
        PROP_Ce50 = [];
        PROP_RMSE = [];
        fitType = @(Hill, EC50, x) 100./(1+(EC50./x).^Hill); % This is the tissue response variation of the hill equation. Used to define our non linear least squares regression
        fitOptions = fitoptions( 'Method', 'NonlinearLeastSquares',...
            'Lower', [1 1], ... % lower guess for hill and Ec50.
            'Upper', [20 10], ... %upper guess for hill and Ec50.
            'StartPoint', [7 6], ... %starting point. 7 and 6 are a decent starting point.
            'Robust', 'LAR'); %least absolute residuals is robust linear least squares fitting method.
        
        for treatment_selector = 1:length(number_of_treatments) % n is big loop for each treatment
            Infusion_dose = csvread(['PROP',participant_number_str,'_T',num2str(number_of_treatments(treatment_selector)),'_Dosing.csv']); % second by second dosing in mcg
            BSR = csvread(['PROP',participant_number_str,'_T',num2str(number_of_treatments(treatment_selector)),'_BSR.csv']); % Is this second by second BSR?
            length_BSR = length(BSR); % is this both the time in seconds and number of samples of BSR?
            BSR = BSR(1:1:end); % burst suppression ratio in percentage throughout the treatment.
            ke0_test_vector = [0.08:0.005:0.20]; % start with 0.08 and work your way up to .2, testing each value with the model
            RMSEVector = [];
            PDVector = []; % gives us a vector range that we are checking
            y_error_saver = [];
            
            for m = 1:length(ke0_test_vector) %guessing and checking Ke0. Testing tthe Ke0 vector to find Ce. m is the length of the vector
                ke0 = ke0_test_vector(m);
                [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, ke0);
                effect_conc = Ce_effect_concentration(1:1:length_BSR)'; % shorten Ce_effect concentration to the length of BSR and transpose???
                [fitobject, gof, output] = fit(effect_conc, BSR, fitType, fitOptions); 
                PD = coeffvalues(fitobject) 
                RMSE = gof.rmse; %root mean square error
                RMSEVector = [RMSEVector RMSE];
                PDVector = [PDVector; PD]; %store what we got as output
                y_error = mean(abs(output.residuals)); %mean of residual absolute error
                y_error_saver = [y_error_saver y_error];
            end
            
            [min_RMSE, min_RMSE_index] = min(RMSEVector);
            Oke0 = ke0_test_vector(min_RMSE_index); %min_RMSE_inded is index of lowest RMSE in Ke0 test vector.
            ke0_test_vector = [Oke0-0.005:0.0025:Oke0+0.005]; %do it again but increase the resolution
            RMSEVector = [];
            PDVector = []; 
            y_error_saver = [];
            
            for m = 1:length(ke0_test_vector)
                ke0 = ke0_test_vector(m);
                [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, ke0);
                effect_conc = Ce_effect_concentration(1:1:length_BSR)';
                [fitobject, gof, output] = fit(effect_conc, BSR, fitType, fitOptions);
                PD = coeffvalues(fitobject);
                RMSE = gof.rmse;
                RMSEVector = [RMSEVector RMSE];
                PDVector = [PDVector; PD];
                y_error = mean(abs(output.residuals)); %mean of residual absolute error
                y_error_saver = [y_error_saver y_error];
            end
            
            [min_RMSE, min_RMSE_index] = min(y_error_saver); %find the best vector with the least error and find the associated Hill, EC50 and Ke0.
            Oke0 = ke0_test_vector(min_RMSE_index);
            ke0_test_vector = [Oke0-0.0005:0.0001:Oke0+0.0005];  % increaseing the resolution again to find Ke0
            RMSEVector = [];
            PDVector = []; 
            y_error_saver = [];
            
            for m = 1:length(ke0_test_vector)
                ke0 = ke0_test_vector(m);
                [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, ke0);
                effect_conc = Ce_effect_concentration(1:1:length_BSR)';
                [fitobject, gof, output] = fit(effect_conc, BSR, fitType, fitOptions);
                PD = coeffvalues(fitobject);
                RMSE = gof.rmse;
                RMSEVector = [RMSEVector RMSE];
                PDVector = [PDVector; PD];
                y_error = mean(abs(output.residuals)); %mean of residual absolute error
                y_error_saver = [y_error_saver y_error];
            end
            
            [min_RMSE, min_RMSE_index] = min(y_error_saver); %find the lowest error Ke0 value
            Oke0 = ke0_test_vector(min_RMSE_index);
            PROP_Ke0 = [PROP_Ke0 Oke0];
            PROP_Hill = [ PDVector(min_RMSE_index,1),PROP_Hill];
            PROP_Ce50 = [PDVector(min_RMSE_index, 2),PROP_Ce50];
            PROP_RMSE = [RMSEVector(min_RMSE_index),PROP_RMSE];
            Pharm_parameters_saver = [PROP_Ke0; PROP_Hill; PROP_Ce50; PROP_RMSE]';
            [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, Oke0);
        end
        
        plot(fitobject, effect_conc, BSR);
        xlabel('Predicted Ce (mcg/mL)');
        ylabel('BSR');
        legend('Observed BSR', 'RMSE Fit', 'Location', "northwest")
        title('PROP PD Model (RMSE Optimization)')
        Pharm_parameters_saver = [PROP_Ke0; PROP_Hill; PROP_Ce50; PROP_RMSE]'
        PD_Model_Filename = sprintf(['PROP',participant_number_str,'_Compiled_PD_Model.csv'], Pharm_parameters_saver);
        
        if  exist(PD_Model_Filename) ~= 0 %#ok<EXIST>
            Pharm_parameters_saver_read = csvread(PD_Model_Filename);
            Pharm_parameters_saver_read = [Pharm_parameters_saver_read; Pharm_parameters_saver ];
            csvwrite(['PROP',participant_number_str,'_Compiled_PD_Model.csv'], Pharm_parameters_saver_read); % Save the modeled parameters as a CSV in an Excel Spreadsheet
        else
            csvwrite(['PROP',participant_number_str,'_Compiled_PD_Model.csv'], Pharm_parameters_saver); % Save the modeled parameters as a CSV in an Excel Spreadsheet
        end
    end
end
%%Generate a dosing recommendation for the participant under
%%consideration. This is done by trying a bunch of different dosing combos
%%and seeing which one gives the desired predicted BSR

if first_treatment == 0
    generate_dosing = input('Would you like to also generate dosing recommendations today? Enter 1 for yes, 0 for no: ');
elseif first_treatment == 1
    generate_dosing = 1;
end

if generate_dosing == 1 %If this is their first treatment, use population Hill and Ce50 and Eleveld Ke0 
    if first_treatment == 1
        PROP_Ke0 = ((Participant_Demographics(4)/70)^(-0.25))*0.146; %Eleveld model Ke0
        PROP_Hill = 6.979;
        PROP_Ce50 = 6.275;
    else %Otherwise, use the most recent treatment PK/PD parameters. 
        Pharm_parameters_saver = csvread(['PROP',participant_number_str,'_Compiled_PD_Model.csv']);
        PROP_Ke0 = Pharm_parameters_saver(end,1);
        PROP_Hill = Pharm_parameters_saver(end,2);
        PROP_Ce50 = Pharm_parameters_saver(end,3);
    end
    minutes=30; 
    Infusion_Saver =  zeros(1,60*minutes); %create an empty vector to save simulated dosing
    Bolus_duration = 30; %standardize the bolus duration to 30 seconds
    infusion_rate = 100; %initialize the infusion rate at a rate smaller than any participant would ever need
    infusion_rate1_saver = num2str(infusion_rate);
    infusion_start = Bolus_duration+1;
    infusion_end = 860;
    Infusion_Saver(infusion_start:infusion_end) = Infusion_Saver(infusion_start:infusion_end) + infusion_rate*Participant_Demographics(4)/60; 
    
    [Cp, Ce] = EleveldPKFun([Participant_Demographics(1) Participant_Demographics(2) Participant_Demographics(3) Participant_Demographics(4)], Infusion_Saver, PROP_Ke0(end));
    fun = @(z,xdata) (100*xdata.^z(1))  ./ (z(2).^z(1)+xdata.^z(1)); %This is the hill equation 
    BSR = fun([PROP_Hill(end) PROP_Ce50(end)],Ce); 
    effect_conc=find(BSR>65);
    x_empty = isempty(effect_conc);
    bolus_stepper = 100;
    infusion_rate_stepper = 800;
    Bolus_mg = bolus_stepper;
    Bolus1_saver = num2str(Bolus_mg);
    Infusion_Saver(1:Bolus_duration) = (bolus_stepper*1000/Bolus_duration);
    figure;
    y = 1;
    while y == 1
        if x_empty == false
            if  effect_conc(1) > 240   %increase the size of the bolus if the threshold isn't reached within four minutes
                bolus_stepper = bolus_stepper + 5;
                disp('Bolus Increased, threshold not reached in time')
            elseif effect_conc(1) < 240
                y = 0;
            end
        elseif x_empty == true
            bolus_stepper = bolus_stepper + 5;
            disp('Bolus Increased, threshold not reached in time')
        end
        Infusion_Saver(1:Bolus_duration) = (bolus_stepper*1000/Bolus_duration);
        [Cp, Ce] = EleveldPKFun([Participant_Demographics(1) Participant_Demographics(2) Participant_Demographics(3) Participant_Demographics(4)], Infusion_Saver, PROP_Ke0(1));
        fun = @(z,xdata) (100*xdata.^z(1))  ./ (z(2).^z(1)+xdata.^z(1)); 
        BSR = fun([PROP_Hill(1) PROP_Ce50(1)],Ce); %insert [Hill, EC50] from EstimateModelRSME
        effect_conc=find(BSR>65);
        x_empty = isempty(effect_conc);
        
        Bolus_Saver = num2str(bolus_stepper);
        infusion_rate_saver = num2str(infusion_rate_stepper);
        plot(BSR(1:end),'b','LineWidth',1.5)
        title( [Bolus_Saver, ' mg bolus for 30 sec, ', infusion_rate_saver, ' mcg/kg/min Infusion @ ', num2str(Participant_Demographics(4)),' Kg'])
        xlim([0 30*60])
        xticks([0, 3*60, 6*60, 9*60, 12*60, 15*60, 18*60, 21*60, 24*60, 27*60, 30*60])
        xticklabels({'0','3', '6', '9', '12', '15', '18', '21', '24', '27', '30'})
        ylabel('Burst Suppresion Ratio (BSR)')
        xlabel('Time (Minutes)'); ylim([0 100])
        yline(70,'g','LineWidth',1.5)
        yline(90,'g','LineWidth',1.5)
        pause(0.01)
    end
    pause(1)
    Bolus_Saver = num2str(bolus_stepper);
    effect_conc = find(BSR>65); 
    x_empty = isempty(effect_conc);
    
    while x_empty == 0
        
        if isempty(effect_conc) == 0
            bolus_stepper = bolus_stepper - 5;
            disp('Bolus decreased')
        else
            break
        end
        Infusion_Saver(1:Bolus_duration) = (bolus_stepper*1000/Bolus_duration);
        [Cp, Ce] = EleveldPKFun([Participant_Demographics(1) Participant_Demographics(2) Participant_Demographics(3) Participant_Demographics(4)], Infusion_Saver, PROP_Ke0(1));
        fun = @(z,xdata) (100*xdata.^z(1))  ./ (z(2).^z(1)+xdata.^z(1)); 
        BSR = fun([PROP_Hill(1) PROP_Ce50(1)],Ce); 
        effect_conc=find(BSR>65);
        x_empty = isempty(effect_conc);
        Bolus_Saver = num2str(bolus_stepper);
        infusion_rate_saver = num2str(infusion_rate_stepper);
        plot(BSR(1:end),'b','LineWidth',1.5)
        title( [Bolus_Saver, ' mg bolus for 30 sec, ', infusion_rate_saver, ' mcg/kg/min Infusion @ ', num2str(Participant_Demographics(4)),' Kg'])
        xlim([0 30*60])
        xticks([0, 3*60, 6*60, 9*60, 12*60, 15*60, 18*60, 21*60, 24*60, 27*60, 30*60])
        xticklabels({'0','3', '6', '9', '12', '15', '18', '21', '24', '27', '30'})
        ylabel('Burst Suppresion Ratio (BSR)')
        xlabel('Time (Minutes)'); ylim([0 100])
        yline(70,'g','LineWidth',1.5)
        yline(90,'g','LineWidth',1.5)
        pause(0.05)
    end
    pause(1)
    x_empty = 0;
    y = 1; %reinitialize the looping variable
    
    while y == 1
        if x_empty == false
            infusion_rate_stepper = infusion_rate_stepper - 10;
            disp('infusion rate decreased')
        elseif x_empty == true
            y = 0;
        end
        Infusion_Saver(31:infusion_end) = infusion_rate_stepper*Participant_Demographics(4)/60;; 
        [Cp, Ce] = EleveldPKFun([Participant_Demographics(1) Participant_Demographics(2) Participant_Demographics(3) Participant_Demographics(4)], Infusion_Saver, PROP_Ke0(1));
        fun = @(z,xdata) (100*xdata.^z(1))  ./ (z(2).^z(1)+xdata.^z(1)); %this is the Hill equation
        BSR = fun([PROP_Hill(1) PROP_Ce50(1)],Ce); 
        effect_conc = find(BSR>80);
        x_empty = isempty(effect_conc);
        infusion_rate_saver = num2str(infusion_rate_stepper);
        plot(BSR(1:end),'b','LineWidth',1.5)
        title( [Bolus_Saver, ' mg bolus for 30 sec, ', infusion_rate_saver, ' mcg/kg/min Infusion @ ', num2str(Participant_Demographics(4)),' Kg'])
        xlim([0 30*60])
        xticks([0, 3*60, 6*60, 9*60, 12*60, 15*60, 18*60, 21*60, 24*60, 27*60, 30*60])
        xticklabels({'0','3', '6', '9', '12', '15', '18', '21', '24', '27', '30'})
        ylabel('Burst Suppresion Ratio (BSR)')
        xlabel('Time (Minutes)'); ylim([0 100])
        yline(70,'g','LineWidth',1.5)
        yline(90,'g','LineWidth',1.5)
        pause(0.01)
    end
    
    pause(1);
    infusion_rate_saver = num2str(infusion_rate_stepper);
    BSR_slope = (BSR(600) - BSR(350)) / (600-350);
    while BSR_slope > 0.1
        [Cp, Ce] = EleveldPKFun([Participant_Demographics(1) Participant_Demographics(2) Participant_Demographics(3) Participant_Demographics(4)], Infusion_Saver, PROP_Ke0(1));
        fun = @(z,xdata) (100*xdata.^z(1))  ./ (z(2).^z(1)+xdata.^z(1)); 
        BSR = fun([PROP_Hill(1) PROP_Ce50(1)],Ce); 
        BSR_slope = (BSR(600) - BSR(350)) / (600-350);
        Infusion_Saver(31:infusion_end) = infusion_rate_stepper*Participant_Demographics(4)/60;; 
        infusion_rate_stepper = infusion_rate_stepper - 5;
        disp('infusion rate decreased, slope')
        
        infusion_rate_saver = num2str(infusion_rate_stepper);
        plot(BSR(1:end),'b','LineWidth',1.5)
        title( [Bolus_Saver, ' mg bolus for 30 sec, ', infusion_rate_saver, ' mcg/kg/min Infusion @ ', num2str(Participant_Demographics(4)),' Kg'])
        xlim([0 30*60])
        xticks([0, 3*60, 6*60, 9*60, 12*60, 15*60, 18*60, 21*60, 24*60, 27*60, 30*60])
        xticklabels({'0','3', '6', '9', '12', '15', '18', '21', '24', '27', '30'})
        ylabel('Burst Suppresion Ratio (BSR)')
        xlabel('Time (Minutes)'); ylim([0 100])
        yline(70,'g','LineWidth',1.5)
        yline(90,'g','LineWidth',1.5)
        pause(0.05)
    end
    
    while BSR_slope < 0
        
        [Cp, Ce] = EleveldPKFun([Participant_Demographics(1) Participant_Demographics(2) Participant_Demographics(3) Participant_Demographics(4)], Infusion_Saver, PROP_Ke0(1));
        fun = @(z,xdata) (100*xdata.^z(1))  ./ (z(2).^z(1)+xdata.^z(1)); 
        BSR = fun([PROP_Hill(1) PROP_Ce50(1)],Ce); 
        BSR_slope = (BSR(600) - BSR(350)) / (600-350);
        Infusion_Saver(31:infusion_end) = infusion_rate_stepper*Participant_Demographics(4)/60;; 
        infusion_rate_stepper = infusion_rate_stepper + 5;
        disp('infusion rate increased')
        
        infusion_rate_saver = num2str(infusion_rate_stepper);
        plot(BSR(1:end),'b','LineWidth',1.5)
        title( [Bolus_Saver, ' mg bolus for 30 sec, ', infusion_rate_saver, ' mcg/kg/min Infusion @ ', num2str(Participant_Demographics(4)),' Kg'])
        xlim([0 30*60])
        xticks([0, 3*60, 6*60, 9*60, 12*60, 15*60, 18*60, 21*60, 24*60, 27*60, 30*60])
        xticklabels({'0','3', '6', '9', '12', '15', '18', '21', '24', '27', '30'})
        ylabel('Burst Suppresion Ratio (BSR)')
        xlabel('Time (Minutes)'); ylim([0 100])
        yline(70,'g','LineWidth',1.5)
        yline(90,'g','LineWidth',1.5)
        pause(0.05)
    end
    infusion_rate_saver = num2str(infusion_rate_stepper);
    pause(2)
    y = 1; %reinitialize the looping variable
    effect_conc = find(BSR>80);
    x_empty = isempty(effect_conc);
    while y == 1
        if x_empty == false
            infusion_rate_stepper = infusion_rate_stepper - 10;
            bolus_stepper = bolus_stepper - 5;
            disp('infusion rate decreased')
            disp('bolus decreased')
        elseif x_empty == true
            y = 0;
        end
        Infusion_Saver(1:Bolus_duration) = (bolus_stepper*1000/Bolus_duration);
        Infusion_Saver(31:infusion_end) = infusion_rate_stepper*Participant_Demographics(4)/60;; 
        [Cp, Ce] = EleveldPKFun([Participant_Demographics(1) Participant_Demographics(2) Participant_Demographics(3) Participant_Demographics(4)], Infusion_Saver, PROP_Ke0(1));
        fun = @(z,xdata) (100*xdata.^z(1))  ./ (z(2).^z(1)+xdata.^z(1)); 
        BSR = fun([PROP_Hill(1) PROP_Ce50(1)],Ce); 
        effect_conc = find(BSR>80);
        x_empty = isempty(effect_conc);
        
        infusion_rate_saver = num2str(infusion_rate_stepper);
        plot(BSR(1:end),'b','LineWidth',1.5)
        title( [Bolus_Saver, ' mg bolus for 30 sec, ', infusion_rate_saver, ' mcg/kg/min Infusion @ ', num2str(Participant_Demographics(4)),' Kg'])
        xlim([0 30*60])
        xticks([0, 3*60, 6*60, 9*60, 12*60, 15*60, 18*60, 21*60, 24*60, 27*60, 30*60])
        xticklabels({'0','3', '6', '9', '12', '15', '18', '21', '24', '27', '30'})
        ylabel('Burst Suppresion Ratio (BSR)')
        xlabel('Time (Minutes)'); ylim([0 100])
        yline(70,'g','LineWidth',1.5)
        yline(90,'g','LineWidth',1.5)
        pause(0.1)
    end
    
    [Cp, Ce] = EleveldPKFun([Participant_Demographics(1) Participant_Demographics(2) Participant_Demographics(3) Participant_Demographics(4)], Infusion_Saver, PROP_Ke0(1));
    fun = @(z,xdata) (100*xdata.^z(1))  ./ (z(2).^z(1)+xdata.^z(1)); 
    BSR = fun([PROP_Hill(1) PROP_Ce50(1)],Ce); 
    plot(BSR(1:end),'b','LineWidth',1.5)
    title( [Bolus_Saver, ' mg bolus for 30 sec, ', infusion_rate_saver, ' mcg/kg/min Infusion @ ', num2str(Participant_Demographics(4)),' Kg'])
    xlim([0 30*60])
    xticks([0, 3*60, 6*60, 9*60, 12*60, 15*60, 18*60, 21*60, 24*60, 27*60, 30*60])
    xticklabels({'0','3', '6', '9', '12', '15', '18', '21', '24', '27', '30'})
    ylabel('Burst Suppresion Ratio (BSR)')
    xlabel('Time (Minutes)'); ylim([0 100])
    yline(70,'g','LineWidth',1.5)
    yline(90,'g','LineWidth',1.5)
    BSR_good_looking = BSR;
    
    infusion_rate_stepper_mg = (infusion_rate_stepper/1000)*(infusion_end/60)*Participant_Demographics(4);
    infusion_rate_stepper = (infusion_rate_stepper_mg*1000/(infusion_end/60)/Participant_Demographics(4));
    infusion_rate_stepper = round( infusion_rate_stepper / 5 ) * 5;
    bolus_stepper = bolus_stepper;
    Infusion_Saver(1:Bolus_duration) = (bolus_stepper*1000/Bolus_duration);
    Infusion_Saver(31:infusion_end) = infusion_rate_stepper*Participant_Demographics(4)/60; 
    
    [Cp, Ce] = EleveldPKFun([Participant_Demographics(1) Participant_Demographics(2) Participant_Demographics(3) Participant_Demographics(4)], Infusion_Saver, PROP_Ke0(1));
    BSR = fun([PROP_Hill(1) PROP_Ce50(1)],Ce); 
    bolus_stepper = round(bolus_stepper/5)*5;
    Bolus_Saver = num2str(bolus_stepper);
    infusion_rate_stepper = round(infusion_rate_stepper/5)*5;
    infusion_rate_saver = num2str(infusion_rate_stepper);
    
    BSR_Good_Time = length(find(BSR>70));
    if BSR_Good_Time >=600
        BSR_Good_Time_Minutes = BSR_Good_Time/60;
        BSR_Good_Time_str = num2str(BSR_Good_Time_Minutes,'%4.2f');
        BSR_Good_Time_Minutes = BSR_Good_Time_str(1:2);
        BSR_Good_Seconds_decimal = str2num(BSR_Good_Time_str(3:5));
        BSR_Good_Seconds = round(BSR_Good_Seconds_decimal*60);
        BSR_Good_Seconds = num2str(BSR_Good_Seconds);
        disp(['BSR in target range for ', BSR_Good_Time_Minutes, ' minutes ', BSR_Good_Seconds, ' seconds'])
        
    else
        BSR_Good_Time_Minutes = BSR_Good_Time/60;
        BSR_Good_Time_str = num2str(BSR_Good_Time_Minutes,'%4.2f');
        BSR_Good_Time_Minutes = BSR_Good_Time_str(1);
        BSR_Good_Seconds_decimal = str2num(BSR_Good_Time_str(2:4));
        BSR_Good_Seconds = round(BSR_Good_Seconds_decimal*60);
        BSR_Good_Seconds = num2str(BSR_Good_Seconds);
        disp(['BSR in target range for ', BSR_Good_Time_Minutes, ' minutes ', BSR_Good_Seconds, ' seconds'])
    end
    
    plot(BSR_good_looking(1:end),'b','LineWidth',1.5)
    title( [Bolus_Saver, ' mg bolus for 30 sec, ', infusion_rate_saver, ' mcg/kg/min Infusion @ ', num2str(Participant_Demographics(4)),' Kg'])
    xlim([0 30*60])
    rectangle('position',[-20,70,2000,20],'EdgeColor', [0,1,0], 'LineWidth',1.5);
    
    xticks([0, 3*60, 6*60, 9*60, 12*60, 15*60, 18*60, 21*60, 24*60, 27*60, 30*60])
    xticklabels({'0','3', '6', '9', '12', '15', '18', '21', '24', '27', '30'})
    ylabel('Burst Suppresion Ratio (BSR)')
    xlabel('Time (Minutes)'); ylim([0 100])
    
    % Calculate the minutes and seconds associated with reaching 70% BSR
    vert = find(BSR_good_looking > 70, 1, 'first');
    vert_time = vert/60;
    vert_time_str = num2str(vert_time,'%4.2f');
    minutes = vert_time_str(1);
    seconds_decimal = str2num(vert_time_str(2:4));
    seconds = seconds_decimal*60;
    seconds = round(seconds, 0);
    seconds = num2str(seconds);
    xline(vert,'r','LineWidth', 1.5, 'Label', ['Projected 70% BSR @' minutes ' minutes ' seconds ' seconds'], 'LabelOrientation', 'horizontal')
end
