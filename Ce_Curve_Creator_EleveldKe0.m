%% Use this scrip to generate effect site concentration curves for each treatment
clc;
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


number_of_treatments = input('Enter the treatment numbers of the past treatments that you want to include in this analysis as a vector\n in the following format [1 2 3 5] : ');
for treatment_selector = 1:length(number_of_treatments) % loop for each treatment

Infusion_dose = csvread(['PROP',participant_number_str,'_T',num2str(number_of_treatments(treatment_selector)),'_Dosing.csv']); % second by second dosing in mcg

Oke0 = ((Participant_Demographics(4)/70)^(-0.25))*0.146; %calculate Ke0 using the equation found in the Eleveld paper.
[Cp_plasma_concentration, Ce_effect_concentration] = EleveldPKFun(Participant_Demographics, Infusion_dose, Oke0);
effect_conc = zeros(1,4000);
for a = 1:length(Ce_effect_concentration)
effect_conc = Ce_effect_concentration(a)';
end
cd(save_directory)
csvwrite(['PROP',participant_number_str,'_T',num2str(number_of_treatments(treatment_selector)),'_Ce_EKe0.csv'], Ce_effect_concentration); % Save effect site concentration as a CSV file.
end

%Plot the BSR as a function of effect site concentration.
%plot(effect_conc, BSR);
% xlabel('Predicted Ce (mcg/mL)');
% ylabel('BSR');
% legend('Observed BSR', 'RMSE Fit', 'Location', "northwest")
% title('PROP PD Model (RMSE Optimization)')
%Pharm_parameters_saver = [PROP_Ke0; PROP_Hill; PROP_EC50; PROP_RMSE]'
%display pharm parameters saver is a more user friendly format
