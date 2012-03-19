clear all; clc; close all;
%% Initializing model
% Car = 2003 Infinity G35 sedan
% max hp = 260 @ 6000 rpm
% max torque = 260 ftlbs @ 4800 rpm
% max engine speed 6600 rpm
% curb weight = 3336 lbs
% wheelbase= 112.2 in
% width (max) = 69 in
% height = 57.7 in
% tires = P205/65VR16
% 1/4 mi time =14.68 @97.08 mph
% top gear rpm=1900
% weight distro 52/48
%guessed
r_eff=.35;
c_rr=.03;
jc=.4;
cg_height=.5*57.7*2.54/100;
rho=1.1;

%given from spec sheet
gear_ratios=[3.54,2.264,1.471,1.000,0.834];
final_drive=3.357;
max_es=6600/9.549;
mass=(3336/2.2);
wheel_base_width=69*2.54/100;
wheel_base_length=112.2*2.54/100;
frontal_area=.75*wheel_base_width*2*cg_height;
cd=.26*.5*rho;


%simulation constants/ICs
road_grade=0;
xdot_naught=0;
ydot_naught=0;
del=0;
yaw_rate=0;

%initializing generic car model, using Car.m object.  Object is created
%with the constants relevant to the car itself.  In seperate functions, the
%simulation parameters and IC's are set, and then the landscape is set to
%be a perfectly flat plane.  The drags and transmission values are also set 
%here.  They are modified througout the assignment. setEta sets the overall 
%engine efficiency term.  This model is used for the remainder of the homework.  
model=Car(mass, wheel_base_width,wheel_base_length,r_eff,cg_height,del);
model=model.SetupSim(xdot_naught, ydot_naught, yaw_rate);
model.setPlanarLandscape(road_grade);
model.setDrags(cd,c_rr,frontal_area);
model.setTransValues(gear_ratios(1), final_drive, max_es,'EffectiveDamping', .001, 'FrictionLoss', .01);
model.setEta(1);

shiftpoints=[6000,6000,6000,6000,inf];
downshiftpoints=[-inf, 3000,3000,3000,3000];
n=0;
e=0;
d=0;
gear=1;
F_brake=0;
model.setDrags(.3,c_rr,frontal_area);
for j=1:1000;
    if model.engine_speed(end)*9.549>shiftpoints(gear)
        gear=gear+1;
        if gear<6
            model.setTransValues(gear_ratios(gear), final_drive, max_es);
            model.stepSim(.01, 0, 0,'PropulsionInputType', 'EngineTorque');
            [~,~,n(end+1),e(end+1),d(end+1)]=model.getGlobalState;
        else
            gear=5;
        end
    end
    model.stepSim(.01, 200,F_brake,'BrakeInputType','PureForce','PropulsionInputType', 'EngineTorque');
    [~,~,n(end+1),e(end+1),d(end+1)]=model.getGlobalState;
end
for j=1:1000;
    if model.engine_speed(end)*9.549>shiftpoints(gear)
        gear=gear+1;
        if gear<6
            model.setTransValues(gear_ratios(gear), final_drive, max_es);
            steerin(j)=.1*pi/180; %sin((2*pi/1000)*j)*1*pi/180;
            model.stepSim(.01, 0, 0,'PropulsionInputType', 'EngineTorque');
            [~,~,n(end+1),e(end+1),d(end+1)]=model.getGlobalState;
        else
            gear=5;
        end
    end
    model.stepSim(.01, 100,F_brake,'BrakeInputType','PureForce','PropulsionInputType', 'EngineTorque');
    steerin(j)=.1*pi/180; %sin((2*pi/1000)*j)*2*pi/180;
    model.steerAdj(steerin(j), 'SteerAngle');
    [~,~,n(end+1),e(end+1),d(end+1)]=model.getGlobalState;
end

figure;
plot3(n, e, d);
xlabel('NORTH');
ylabel('EAST');
zlabel('DOWN');

figure;
subplot(3,1,1);
plot(0:.01:length(model.xdot)*.01-.01, model.xdot.*2.23693629);
title('Speed over time');
ylabel('MPH');
xlabel('Time (seconds)');
subplot(3,1,2);
plot(0:.01:length(model.xdot)*.01-.01, model.engine_speed.*9.549);
title('Engine Speed over time');
ylabel('RPM');
xlabel('Time (seconds)');
subplot(3,1,3);
plot(0:.01:length(model.xdot)*.01-.01, model.xdoubledot./9.81);
title('Acceleration over time');
ylabel('G''s');
xlabel('Time (seconds)');

