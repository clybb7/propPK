%Jason Huang, University of Utah. 2021.

function [Cp, Ce, PK] = EleveldPKFun(Participant_Demographics, infusion_rate, ke0)
%% ELEVELDPKFUN Convert Infusion to Concentration
% Imported *Patient Demographics* = [sex (0 or 1 F/M), height (cm), age (yr), weight (kg)]
% *Infusion* is the second-by-second administered drug rate in mcg/s
% *ke0* must be defined by the user. I have a seperate function to optimize.
% *Cp* is the second-by-second, predicted plasma concentration (mcg/mL)
% *Ce* is the second-by-second, predicted effect-site concentration (mcg/mL),
% based on the ke0
% *PK* outputs parameters: central volume, fast-acting volume, slow-acting volume;
% central clearing rate, fast-acting clearing rate, and slow-acting clear rate.

% See: Eleveld DJ, Colin P, Absalom AR, Struys MMRF. Pharmacokinetic–pharmacodynamic
% model for propofol for broad application in anaesthesia and sedation. British Journal of Anaesthesia 2018;120:942–59.

% Unpack Patient demographics
Sex = Participant_Demographics(1); 
Height = Participant_Demographics(2); 
Age = Participant_Demographics(3); 
Weight = Participant_Demographics(4);

% PK Parameters obtained from Eleveld et al. 
O1 = 6.28; % V1ref (L)
O2 = 25.5; % V2ref (L)
O3 = 273;  % V3ref (L)
O4 = 1.79; % C:ref (male) L*min^-1
O5 = 1.75; % Q2ref L*min^-1
O6 = 1.11; % Q3ref L*min^-1
O7 = 0.191;% Typical Residual Error
O8 = 42.3; %CL maturation E50
O9 = 9.06; %CL maturation slope
O10 = -0.0156; %Smaller V2 with age
O11 = -0.00286; %Lower CL with age
O12 = 33.6; % Weight for 50% of maximal V1
O13 = -0.0138; %Smaller V3 with age
O14 = 68.3; % Maturation of Q3
O15 = 2.10; % CLref (female)
O16 = 1.30; % Higher Q2 for maturation of Q3
O17 = 1.42; % V1 venous samples (children)
O18 = 0.68; % Higher Q2 venous samples

% Reference Individual, a generic individual based on populations estimates
ref.sex = 1; % male
ref.height = 170; % height in cm
ref.age = 35; 
ref.weight = 70; % weight in Kg
ref.BMI = ref.weight/(ref.height/100)^2; %body mass index
ref.FFM = [0.88 + (1-0.88)/(1+(ref.age/13.4)^-12.7)]*[(9270*ref.weight)/(6680+216*ref.BMI)]; % fat free mass
ref.PMA = ref.age*52 + 40; % post-menstrual age. try 52.1429

% PK Equations
BMI = Weight/((Height/100)^2); %BMI calculation of this patient

if Sex == 1
    FFM = [0.88 + (1-0.88)/(1+(Age/13.4)^-12.7)]*[(9270*Weight)/(6680+(216*BMI))]; %fat free mass calculation from Al-Sallami et al, cited in Eleveld et al.
else
    FFM = [1.11 + (1-1.11)/(1+(Age/7.1)^-1.1)]*[9270*Weight/(8780+(244*BMI))];
end

%PK equations taken directly from the eleveld model. Many of them are not
%clearly defined in the paper. 
Fageing = @(x) exp(x*(Age-ref.age)); %takes into account the effect of aging on the pharmacokinetics, using the population model of eleveld
Fsigmoid = @(x, EC50, gam) (x^gam)/(x^gam + EC50^gam); % no clue what this does
Fcentral = @(x) Fsigmoid(x, O12, 1); 
PMA = Age*52 + 40; %post menstrual age
FCLmaturation = Fsigmoid(PMA, O8, O9); % No clue
ref.CLmat = Fsigmoid(ref.PMA, O8, O9); % Don't know
FQ3maturation = Fsigmoid(Age*52 + 40, O14, 1); % 
ref.Q3mat = Fsigmoid(ref.age*52 + 40, O14, 1); % 
V1arterial = O1*(Fcentral(Weight)/Fcentral(ref.weight)); % arterial volume, representing the central compartment of the three compartment model
V2_rapid = O2*(Weight/ref.weight)*Fageing(O10); % volume of the rapidly equilibrating compartment of the three compartment model
ref.V2 = O2; 
V3_slow = O3*(FFM/ref.FFM); % volume of the slowing equilibrating compartment of the three compartment model
ref.V3 = O3; 

if Sex == 1 % if male, _
    OCL = O4; % different constants used for male vs female
else
    OCL = O15;
end

CL = OCL*((Weight/ref.weight)^0.75)*(FCLmaturation/ref.CLmat); 
Q2a = O5*((V2_rapid/ref.V2)^0.75)*(1+O16*(1-FQ3maturation)); 
Q3 = O6*((V3_slow/ref.V3)^0.75)*(FQ3maturation/ref.Q3mat); 
PK = [V1arterial V2_rapid V3_slow CL Q2a Q3]; % V1=volume of compartment 1, latter 3 are general clearance rates. This code all comes from the paper.

% Execute Simulation
step = 1/60; % one second
mass_prop_central_1 = 0; % create variable to hold the mass of propofol in the central compartment; initialize with 0 g
mass_prop_rapid_2 = 0; % create variable to hold the mass of propofol in the rapid compartment; initialize with 0 g
mass_prop_slow_3 = 0; % create variable to hold the mass of propofol in the slow compartment; initialize with 0 g
Conc1_arterial = mass_prop_central_1/V1arterial; % concentration of propofol in the blood (central compartment)
Conc2_rapid = mass_prop_rapid_2/V2_rapid; % concentration of propofol in the rapidly equilibrating compartment
Conc3_slow = mass_prop_slow_3/V3_slow; % concentration of propofol in the slowly equilibrating compartment
Ce = 0; % Effect site concentration begins at 0
store1 = [Conc1_arterial]; % store and track the concentration of the different compartments over time .
store2 = [Conc2_rapid]; 
store3 = [Conc3_slow];
store4 = [Ce];
k10 = PK(4)/PK(1); % clearance rates between compartments. Compartment 1 to compartment 0.
k12 = PK(5)/PK(1);
k13 = PK(6)/PK(1); 
k21 = PK(5)/PK(2);
k31 = PK(6)/PK(3);

for n = 1:length(infusion_rate) %load in second-by-second amount of drug being administered to the central compartment. 
    
    %the rate of change of propofol in the various compartments
    d_mass_prop_cental_1 = [mass_prop_rapid_2*k21 + mass_prop_slow_3*k31 - mass_prop_central_1*(k10+k12+k13)]*step + infusion_rate(n); % the central compartment is where drug is infused into since it includes blood
    d_mass_prop_rapid_2 = [mass_prop_central_1*k12 - mass_prop_rapid_2*k21]*step; % change of prop concentration in rapid compartment
    d_mass_prop_slow_3 = [mass_prop_central_1*k13 - mass_prop_slow_3*k31]*step; % change of prop concentration in slow compartment
    
    mass_prop_central_1 = mass_prop_central_1 + d_mass_prop_cental_1; % tracking the mass of drug in each compartment.
    mass_prop_rapid_2 = mass_prop_rapid_2 + d_mass_prop_rapid_2;
    mass_prop_slow_3 = mass_prop_slow_3 + d_mass_prop_slow_3;
    
    Conc1_arterial = mass_prop_central_1/V1arterial; 
    Conc2_rapid = mass_prop_rapid_2/V2_rapid; 
    Conc3_slow = mass_prop_slow_3/V3_slow;
    dCe = [ke0*(Conc1_arterial - Ce)]*step;
    Ce = Ce + dCe;
    
    store1 = [store1 Conc1_arterial]; % variables to store the amount of drug in each of the three compartments
    store2 = [store2 Conc2_rapid];
    store3 = [store3 Conc3_slow];
    store4 = [store4 Ce]; 
end
Cp = store1/1000; 
Ce = store4/1000;
end
