%% Script for the comparison of estimated second by second BSR with observed second by second BSR during high dose propofol treatments
clear
fitType = @(Hill, EC50, x) 100./(1+(EC50./x).^Hill); %This defines the equation used by the Fit function later, the Hill Equation used as our regression equation
fitOptions = fitoptions( 'Method', 'NonlinearLeastSquares', ...
    'Lower', [0 0], ... % minimum lower guess for hill and Ec50
    'Upper', [25 25], ... %maximum upper guess for hill and Ec50
    'Start', [7 6], ... %starting point of 7 and 6 for Hill and EC50
    'Robust', 'LAR');

participant_numbers = [21 22 24 25 28 29 31 34 36 37 38 39 41 44 45 46 47 48 50 54]; %Load all participants
total_treat_vector_expander = 1; %used to iteratively expand variables to save values across all treatments of all participants

Ke0_Saver_Full = []; % these save all retrospectively determined PK/PD parameters
Hill_Saver_Full = [];
EC50_Saver_Full = [];

for participant_selector = 1:1:length(participant_numbers) %initialize a loop that iterates for each participant
    Ke0_withinParticipant_saver = []; %initialize variables that save PK/PD parameters for each treatment of each participant, which reset for each iterated participant
    Hill_withinParticipant_saver = [];
    EC50_withinParticipant_saver = [];
    
    participant_number = participant_numbers(participant_selector); %select current participant
    participant_number_str = num2str(participant_number);
    
    try % set save/read directory
        save_directory = sprintf('C:\\Users\\u6027582\\OneDrive\\Masters Materials\\Matlab\\PROP%d_Matlab\\', participant_number);
        cd(save_directory);
        addpath( 'C:\Users\u6027582\OneDrive\Masters Materials\Matlab\');
        addpath(save_directory);
    catch
    end
    Participant_Demographics = csvread(['PROP',participant_number_str,'_Demographics.csv']); %read participant information
    number_of_treatments = [];
    try % try loading each possible treatment number of the participant under considersation
        csvread(['PROP',participant_number_str,'_T',num2str(1),'_Dosing.csv']); % read second by second dosing in mcg
        number_of_treatments = [1];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(2),'_Dosing.csv']);
        number_of_treatments = [number_of_treatments 2];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(3),'_Dosing.csv']); 
        number_of_treatments = [number_of_treatments 3];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(4),'_Dosing.csv']);
        number_of_treatments = [number_of_treatments 4];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(5),'_Dosing.csv']); 
        number_of_treatments = [number_of_treatments 5];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(6),'_Dosing.csv']); 
        number_of_treatments = [number_of_treatments 6];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(7),'_Dosing.csv']); 
        number_of_treatments = [number_of_treatments 7];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(8),'_Dosing.csv']); 
        number_of_treatments = [number_of_treatments 8];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(9),'_Dosing.csv']); 
        number_of_treatments = [number_of_treatments 9];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(10),'_Dosing.csv']); 
        number_of_treatments = [number_of_treatments 10];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(11),'_Dosing.csv']); 
        number_of_treatments = [number_of_treatments 11];
    catch
    end
    try
        csvread(['PROP',participant_number_str,'_T',num2str(12),'_Dosing.csv']); 
        number_of_treatments = [number_of_treatments 12];
    catch
    end
    within_treatment_vector_expander = 1;
    
    fitType = @(Hill, EC50, x) 100./(1+(EC50./x).^Hill); % This is the tissue response variation of the hill equation. Used to define our non linear least squares regression
    fitOptions = fitoptions( 'Method', 'NonlinearLeastSquares', ...
        'Lower', [0 0], ... 
        'Upper', [25 25], ... 
        'StartPoint', [7 6], ... 
        'Robust', 'LAR'); 
    first_treat = 1;
    
    for treatment_selector = 1:length(number_of_treatments) % n is big loop for each treatment
        BSR_100 = [];
        BSR_actual = zeros(1,4000); % initialize BSR saving vectors
        BSR_zeros = zeros(1,4000);
        
        Infusion_dose = csvread(['PROP',participant_number_str,'_T',num2str(number_of_treatments(treatment_selector)),'_Dosing.csv']); % second by second dosing in mcg
        BSR = csvread(['PROP',participant_number_str,'_T',num2str(number_of_treatments(treatment_selector)),'_BSR.csv']); % second by second actual BSR
        length_BSR = length(BSR); %save the length of the actual BSR
        BSR = BSR(1:1:end); % burst suppression ratio throughout the treatment.
        ke0_test_vector = [0.00:0.05:3]; 
        RMSEVector = []; % initialize error saving vectors
        y_error_saver = [];
        PDVector = [];
        
        for Ke0_test_interator = 1:length(ke0_test_vector) %guessing and checking different Ke0 values
            ke0 = ke0_test_vector(Ke0_test_interator); 
            [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, ke0);
            effect_conc = Ce_effect_concentration(1:1:length_BSR)'; % shorten Ce_effect concentration to the length of BSR and transpose
            [fitobject, gof, output] = fit(effect_conc, BSR, fitType, fitOptions); % regression of effect vs effect_conc curve
            PD = coeffvalues(fitobject); %grab pharmacokinetic parameters, Hill and EC50 from our regression of the effect vs effect conc curve
            PDVector = [PDVector; PD]; %store PD parameters
            y_error = mean((output.residuals).^2); % calculate the mean square error from this iteration of the regression with our selected Ke0 value
            y_error_saver = [y_error_saver y_error]; % save error from all iterations of Ke0 for this treatment
        end
        
        [min_RMSE, min_RMSE_index] = min(y_error_saver); % select least mean square error
        Oke0 = ke0_test_vector(min_RMSE_index); % grab Ke0 associated with least mean square error
        ke0_test_vector = [Oke0-0.05:0.025:Oke0+0.05]; % iterate Ke0 again, but with better resolution
        RMSEVector = [];
        PDVector = []; 
        y_error_saver = [];
        
        for Ke0_test_interator = 1:length(ke0_test_vector)
            ke0 = ke0_test_vector(Ke0_test_interator);
            [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, ke0);
            effect_conc = Ce_effect_concentration(1:1:length_BSR)';
            [fitobject, gof, output] = fit(effect_conc, BSR, fitType, fitOptions);
            PD = coeffvalues(fitobject);
            PDVector = [PDVector; PD];
            y_error = mean((output.residuals).^2);
            y_error_saver = [y_error_saver y_error];
        end
        
        [min_RMSE, min_RMSE_index] = min(y_error_saver);
        Oke0 = ke0_test_vector(min_RMSE_index);
        ke0_test_vector = [Oke0-0.005:0.001:Oke0+0.005];
        RMSEVector = [];
        PDVector = []; 
        y_error_saver = [];
        
        for Ke0_test_interator = 1:length(ke0_test_vector)
            ke0 = ke0_test_vector(Ke0_test_interator);
            [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, ke0);
            effect_conc = Ce_effect_concentration(1:1:length_BSR)';
            [fitobject, gof, output] = fit(effect_conc, BSR, fitType, fitOptions);
            PD = coeffvalues(fitobject);
            PDVector = [PDVector; PD];
            y_error = mean((output.residuals).^2);
            y_error_saver = [y_error_saver y_error];
        end
        [min_RMSE, min_RMSE_index] = min(y_error_saver); 
        Oke0 = ke0_test_vector(min_RMSE_index); % save lowest overall error Ke0
        
        [Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, Oke0); % compute the final effect site concentration using the most accurate Ke0 you have calculated
        effect_conc = Ce_effect_concentration(1:1:length_BSR)'; 
        [fitobject, gof, output] = fit(effect_conc, BSR, fitType, fitOptions);
        
        figure; %plot best fit Ke0 hysteresus (effect vs effect site concentration curve)
        plot(fitobject, effect_conc, BSR)
        legend('off')
        title('Final BSR')
        hold on;
        
        if first_treat == 1
            PROP_Ke0 = ((Participant_Demographics(4)/70)^(-0.25))*0.146; %calculate Ke0 using the equation found in the Eleveld paper.
            if participant_number == 21
                PROP_Hill = 6.314; % select the median Hill and EC50 values from the first treatment of all participants, excluding the participant under consideration
                PROP_EC50 =6.206;
            elseif participant_number == 22
                PROP_Hill = 6.601; 
                PROP_EC50 =6.167;
            elseif participant_number == 24
                PROP_Hill = 6.601; 
                PROP_EC50 =6.206;
            elseif participant_number == 25
                PROP_Hill = 6.601; 
                PROP_EC50 =6.206;
            elseif participant_number == 28
                PROP_Hill = 6.314; 
                PROP_EC50 =6.206;
            elseif participant_number ==29
                PROP_Hill = 6.314; 
                PROP_EC50 =6.167;
            elseif participant_number == 31
                PROP_Hill = 6.314;
                PROP_EC50 =6.167;
            elseif participant_number == 34
                PROP_Hill = 6.314; 
                PROP_EC50 =6.167;
            elseif participant_number == 36
                PROP_Hill = 6.314;
                PROP_EC50 =6.167;
            elseif participant_number == 37
                PROP_Hill = 6.601; 
                PROP_EC50 =6.206;
            elseif participant_number ==  38
                PROP_Hill = 6.314; 
                PROP_EC50 =6.206;
            elseif participant_number ==  39
                PROP_Hill = 6.314; 
                PROP_EC50 =6.167;
            elseif participant_number == 41
                PROP_Hill = 6.314;
                PROP_EC50 =6.167;
            elseif participant_number == 44
                PROP_Hill = 6.601; 
                PROP_EC50 =6.167;
            elseif participant_number == 45
                PROP_Hill = 6.601;
                PROP_EC50 =6.167;
            elseif participant_number == 46
                PROP_Hill = 6.314;
                PROP_EC50 =6.206;
            elseif participant_number == 47
                PROP_Hill = 6.601;
                PROP_EC50 =6.206;
            elseif participant_number == 48
                PROP_Hill = 6.601; 
                PROP_EC50 =6.167;
            elseif participant_number == 50
                PROP_Hill = 6.601;
                PROP_EC50 =6.206;
            elseif participant_number == 54
                PROP_Hill = 6.314;
                PROP_EC50 =6.167;
            end
        else
            PROP_Ke0 = Oke0; % for treatments 2-n, use the Ke0 Hill and CE50 values from the most recent treatment
            PROP_Hill = PDVector(min_RMSE_index,1);
            PROP_EC50 = PDVector(min_RMSE_index, 2);
        end
        
        Ke0_withinParticipant_saver = [Ke0_withinParticipant_saver, Oke0];% save PK/PD parameters from all treatments of this participant
        Hill_withinParticipant_saver = [Hill_withinParticipant_saver, PDVector(min_RMSE_index,1)];
        EC50_withinParticipant_saver = [EC50_withinParticipant_saver, PDVector(min_RMSE_index, 2)];
        
        if first_treat == 1
            Ke0_Saver_Full = [Ke0_Saver_Full,Oke0]; % save all PK/PD parameters
            Hill_Saver_Full = [Hill_Saver_Full,PDVector(min_RMSE_index,1)];
            EC50_Saver_Full = [EC50_Saver_Full,PDVector(min_RMSE_index, 2)];
            
            % for first treatment use median population Hill and EC50, along with Eleveld Ke0 
            [Cp, Ce] = EleveldPKFun([Participant_Demographics(1) Participant_Demographics(2) Participant_Demographics(3) Participant_Demographics(4)], Infusion_dose, PROP_Ke0);
            fun = @(z,xdata) (100*xdata.^z(1))  ./ (z(2).^z(1)+xdata.^z(1)); 
            BSR = fun([PROP_Hill PROP_EC50],Ce); 
            first_treat = 0;
        else          
            Ke0_Saver_Full = [Ke0_Saver_Full,PROP_Ke0]; %save all PK/PD parameters
            Hill_Saver_Full = [Hill_Saver_Full,PROP_Hill];
            EC50_Saver_Full = [EC50_Saver_Full,PROP_EC50];
            
            % for subsequent treatments, use most recent treatment PK/PD
            [Cp, Ce] = EleveldPKFun([Participant_Demographics(1) Participant_Demographics(2) Participant_Demographics(3) Participant_Demographics(4)], Infusion_dose, Ke0_withinParticipant_saver(end-1));
            fun = @(z,xdata) (100*xdata.^z(1))  ./ (z(2).^z(1)+xdata.^z(1)); %this specifies Hill curve
            BSR = fun([Hill_withinParticipant_saver(end-1) EC50_withinParticipant_saver(end-1)],Ce); %insert [Hill, EC50] from EstimateModelRSME
        end
        
        for i = 1:length(BSR)
            BSR_zeros(i) = BSR(i);
        end
        
        BSR_read = csvread(['PROP',participant_number_str,'_T',num2str(number_of_treatments(treatment_selector)),'_BSR.csv']); % read in actual BSR of the treatment under consideration
        BSR_read = BSR_read'; 
        BSR_read = movmean(BSR_read,60); %filter with a moving mean to smooth
        
        for i = 1:length(BSR_read)
            BSR_actual(i) = BSR_read(i); % make BSR estimated the same length as the read BSR
        end

        BSR_threshold_proj = 30; % set the threshold above which we will consider BSR
        BSR_threshold_actual = 30;
        curve_vector = find(BSR_actual>BSR_threshold_actual); % find all points that BSR is above threshold
        
        % begin generating square versions of our treatment curves. Compute
        % the width first
        if isempty(curve_vector)
            curve_vector = zeros(1,1000);
        end
        curve_start = curve_vector(1);
        curve_end = curve_vector(end);
        width_actual = curve_end - curve_start;
        curve_vector = find(BSR_zeros>BSR_threshold_proj);
        if isempty(curve_vector)
            curve_vector = zeros(1,1000);
        end
        % do the same for the estimated BSR curve
        curve_start = curve_vector(1);
        curve_end = curve_vector(end);
        width_proj = curve_end - curve_start;
        
%         figure; % plot actual vs estimated BSR curves
        title( ['PROP', participant_number_str, ' Treament ', num2str(treatment_selector)])
        ylim([0 100]);
        xlim([0 2000]);
        xlabel('Time From First Bolus (Seconds)')
        ylabel('Burst Supression Ratio')
        hold on;
        plot(BSR_actual, 'Linewidth',2.0)
        plot(BSR_zeros,'Linewidth',2.0)
        hold on;
        
        %compute height of estimated BSR and actual BSR curves
        error_abs = (abs(BSR_actual - BSR_zeros));
        Error_saver(within_treatment_vector_expander) = trapz(error_abs)/1000;
        a_under_curve_actual = trapz(BSR_actual);
        a_under_curve_proj =  trapz(BSR_zeros);
        norm_height_actual = round(a_under_curve_actual/width_actual,2);
        norm_height_proj = round(a_under_curve_proj/width_proj,2);
        within_treatment_vector_expander = within_treatment_vector_expander+1;

        infusion_duration = length(find(Infusion_dose>0));
        if (width_proj-infusion_duration)>600 % if the width of our curves (BSR>30) extends more than 10 minutes beyond the end of infusion, cut off the projected curve width
            width_proj = infusion_duration+600;
        end
        
        % create vectors that save the actual height and width and
        % estimated height and width of curves for every treatment of every participant
        participant_number = num2str(participant_number);
        height_actual_saver(total_treat_vector_expander) = norm_height_actual;
        height_proj_saver(total_treat_vector_expander) =  norm_height_proj;
        width_actual_saver(total_treat_vector_expander) = width_actual;
        width_proj_saver(total_treat_vector_expander) = width_proj;
        total_treat_vector_expander = total_treat_vector_expander+1;
        width_proj_square = zeros(1,4000);
        width_actual_square = zeros(1,4000);
        
        for w_index = 1:1:width_proj % move our square curves out a bit from the y axis to make plot look nicer
            width_proj_square(w_index+20) = norm_height_proj;
        end
        for w_index = 1:1:width_actual
            width_actual_square(w_index+20) = norm_height_actual;
        end
        
        figure;
        title( ['PROP', participant_number_str, ' Treament ', num2str(treatment_selector)])
        xlabel('Time From First Crossing of 30% Threshold (Seconds)')
        ylabel('Burst Supression Ratio')
        ylim([0 100]);
        xlim([0 2000]);
        hold on;
        plot(width_actual_square, 'Linewidth',2.0)
        plot(width_proj_square, 'Linewidth',2.0)
        disp(['Participant Number: ',participant_number])
    end
    
    width_percent_error = 0; % reset variables for next loop iteration
    height_percent_error = 0;
    percent_overlap = 0;
    Error_saver = 0;
    percent_height = 0;
end