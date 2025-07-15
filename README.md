# propPK

This code contains a software package largely written by Carter Lybbert and partially by Jason Huang, at University of Utah. It is an adaptation of the Eleveld PK/PD model and an evaluation of its accuracy when applied to high-dose propofol treatments for depression. The results of this analysis are published here: 

Lybbert C, Huang J, Jones KG, Mickey BJ, Tadler S, Odell D, Stanford J, Kuck K. Clinical validation of an adapted Eleveld Model for high-dose propofol treatments for depression. Journal of Clinical Monitoring and Computing. 2023 Oct;37(5):1369-1377. doi: 10.1007/s10877-023-00986-7. 

Descripton of the code contained herein:
EleveldPKFun: This code is a function implementation of the three compartment Eleveld PK/PD model, explained here: 

Eleveld, D. J., Colin, P., Absalom, A. R., & Struys, M. M. R. F. (2018). Pharmacokinetic–pharmacodynamic model for propofol for broad application in anaesthesia and sedation. British Journal of Anaesthesia, 120(5), 942–959. https://doi.org/10.1016/j.bja.2018.01.018

It requires a second-by-second propofol infusion rate, participant age, height, weight and sex, and a value for Ke0 to produce a second-by-second effect site concentration. 

Ce_Curve_Creator: This code is used to output an effect site concentration vector as a function of treatment time, making primary use of the EleveldPKFun function. 
Participant_Doser: This code is used to:
1.	Create a vector describing the second-by-second dosing of a recent propofol treatment
 and
2. Generate a dosing recommendation for all treatments of all participants.
   
Projection_Accuracy_Analysis: This code is used to analyze the effectiveness of the adapted Eleveld PK/PD model to predict second-by-second burst suppression ratio by optimizing Ke0, Hill and Ce50. 
