%% getSubjectStatistics

clear all; close all; clc;

%% Select subjects

subjects = 1;


%% Do stuff

age = NaN(1,length(subjects));
mass = NaN(1,length(subjects));
height = NaN(1,length(subjects));
for isubj = subjects
    getInfo_Template;
    
    age(isubj) = subjage;
    mass(isubj) = subjmass;
    height(isubj) = subjheight;
    
    clear subjage subjmass subjheight
end


age_m = nanmean(age);
age_std = nanstd(age);

mass_m = nanmean(mass);
mass_std = nanstd(mass);

height_m = nanmean(height);
height_std = nanstd(height);


disp(['Age: ' num2str(age_m) '+-' num2str(age_std)]);
disp(['Mass: ' num2str(mass_m) '+-' num2str(mass_std)]);
disp(['Height: ' num2str(height_m) '+-' num2str(height_std)]);