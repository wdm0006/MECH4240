clear all; clc;
%% Simulating a local gradient accelleration, helix
% road_grade=-(20/360)*2*pi;
% mass=1500; 
% r_eff=.35; 
% terminal_velocity=140*.44704;
% xdot_naught=20;
% ydot_naught=1;
% del=(5/360)*2*pi;
% yaw_rate=(50/360)*2*pi;
% wheel_base_width=1.5;
% wheel_base_length=2.5;
% c_rr=0.01;
% cg_height=.5;
% frontal_area=2;
% 
% %initializing generic car model
% model=Car(road_grade,mass,terminal_velocity,xdot_naught,ydot_naught,...
%             del,yaw_rate,wheel_base_width,wheel_base_length,r_eff,c_rr,cg_height,frontal_area);
% 
% x=0;
% xdot=xdot_naught;
% yaw=0;
% y=0;
% ydot=ydot_naught;
% n=0;
% e=0;
% d=0;
% 
% F_brake=-100;
% 
% model.setCD(.3);
% for j=1:1000;
%     model.stepBrakeSim(.01, 2, F_brake);
%     [yaw(end+1), x(end+1),xdot(end+1),~,y(end+1),ydot(end+1),~]=model.getState;
%     [~,~,n(end+1),e(end+1),d(end+1)]=model.getGlobalState;
% end
% 
% figure; 
% plot3(n, e, d);
% xlabel('NORTH');
% ylabel('EAST');
% zlabel('DOWN');
%% Simulating driving/turning on random surface
tsize=250;
terrain=genSampleTerrain(tsize+1, tsize+1);
road_grade=0;
mass=1500; 
r_eff=.35; 
terminal_velocity=140*.44704;
xdot_naught=20;
ydot_naught=.1;
del=(5/360)*2*pi;
yaw_rate=(50/360)*2*pi;
wheel_base_width=1.5;
wheel_base_length=2.5;
c_rr=0.01;
cg_height=.5;
frontal_area=2;

%initializing generic car model
model=Car(road_grade,mass,terminal_velocity,xdot_naught,ydot_naught,...
            del,yaw_rate,wheel_base_width,wheel_base_length,r_eff,c_rr,cg_height,frontal_area);

x=0;
xdot=xdot_naught;
yaw=0;
y=0;
ydot=ydot_naught;
n=0;
e=0;
d=0;

F_brake=-100;

model.setCD(.3);
model.setPlanarLandscape(terrain);
for j=1:1000;
    model.stepTerrainBrakeSim(.1, 2, F_brake);
    [yaw(end+1), x(end+1),xdot(end+1),~,y(end+1),ydot(end+1),~]=model.getState;
    [~,~,n(end+1),e(end+1),d(end+1)]=model.getGlobalState;
end
n=n+tsize/2;
e=e+tsize/2;

figure; 
surfl(terrain);
hold on;
plot3(n(2:end), e(2:end), d(2:end),'r');
xlabel('NORTH');
ylabel('EAST');
zlabel('DOWN');
xlim([0,tsize+1]);
ylim([0,tsize+1]);