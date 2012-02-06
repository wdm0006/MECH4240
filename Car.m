%% Car class
% Car class represents a vehicle rolling on an arbitrary surface.  Initial
% conditions and inputs allow various simulations, and get functions allow
% the monitoring and analysis of the system's response.

%Version 3
% -Order parameter in simulation deprecated, neglect air and rolling
% resistance by setting corresponding coefficients to zero, using
% setDrags() function
% -Added in input force designation in simulation, not yet functioning,
% uses the inputs as pure forces no matter what, frame work in place for it
% though.

%Version 2 notes
% -Added in a braking force input to simulation, negative is propulsive
% -Deprecated coasting sim, just use plain sim with no input
% -Added in terrain map ability, road grade integrated into that,
% deprecated as a seperate parameter
% -Seperated Constructor and simulation setup
% -Terminal Velocity property deprecated, user inputs Cd and Af


classdef Car < handle
    
    %% Class Properties
    properties
        %System Properties: the parameters which define the car itself.
        %Mass: kg
        mass
        %Gravity is assumed to be 9.81 m/s
        g=9.81;
        %Drag Coefficient
        Cd
        %Width of the wheels of the car
        width
        %Wheel to wheel length of car
        length
        %Effective Radius of drive wheel
        r_eff
        %Coefficient of rolling resistance
        C_rr
        %Height in meters of the car's center of gravity
        cg_height
        %Frontal Area
        fa
        %trans and diff ratios
        N1=1;
        N2=1;
        %trans loss coeffs
        transLoss=0;
        transB=0;
        
        
        %Parameters for the simulation of the system
        x
        xdot
        xdotnaught
        xdoubledot
        y
        ydot
        ydoubledot
        ydotnaught
        yaw
        yawdot
        yawdotnaught
        del
        N
        E
        D
        terrain
    end
    %% Methods
    methods
        %% Constructor
        %Sets parameters of the car itself
        function obj=Car(m,wbw, wbl,rf,cgh,d)
            obj.mass=m;
            obj.del=d;
            obj.length=wbl;
            obj.width=wbw;
            obj.r_eff=rf;
            obj.C_rr=0;
            obj.cg_height=cgh;
            obj.fa=0;
            obj.Cd=0;
            
            %simulation parameters set to 0 just in case
            obj.x=0;
            obj.xdot=0;
            obj.xdotnaught=0;
            obj.xdoubledot=0;
            obj.y=0;
            obj.ydot=0;
            obj.ydoubledot=0;
            obj.ydotnaught=0;
            obj.yaw=0;
            obj.yawdot=0;
            obj.yawdotnaught=0;
            obj.N=0;
            obj.D=0;
            obj.E=0;
        end
        %% Setup
        %sets the initial conditions of the simulation
        function obj=SetupSim(obj,xdot0,ydot0, yawrate )
            obj.x=0;
            obj.xdot=xdot0;
            obj.xdotnaught=xdot0;
            obj.xdoubledot=0;
            obj.y=0;
            obj.ydot=ydot0;
            obj.ydoubledot=0;
            obj.ydotnaught=ydot0;
            obj.yaw=0;
            obj.yawdot=yawrate;
            obj.yawdotnaught=yawrate;
            obj.N=0;
            obj.D=0;
            obj.E=0;
        end
        %% Step Forward Terrain-Braking Simulation
        %Does support rolling resistance, does not support yaw double dot.
        %turning values are still pretty primative, just uses in inputed
        %ydot, and moves that direction.  A lateral force summation should
        %be used in the nearish future, that just gets ignored if there are
        %not enough parameters or there is an input ydot or ydoubledot.
        %Handles multiple forms of braking and propulsion input, but for
        %now does the same thing no matter what.
        function out=stepSim(obj,time_step,F_prop,F_brake,varargin)
            [vn1, ve1]=obj.getGlobalVelocity;
            gradient=obj.getLocalGradient;
            
            %INPUT TYPE PARSING********************************************
            if mod(size(varargin,2),2)==1 && size(varargin,2)>0
                disp(varargin);
                disp(size(varargin,2));
                error('Must have an even number of force designation parameters');
            end
            %If one input designation was given
            if size(varargin,2)==2
                %if it is designating the brake input type, prop assumed to
                %be a pure force;
                if strcmp(varargin{1,1},'BrakeInputType')
                    if strcmp(varargin{1,2},'PureForce')
                        obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                    elseif strcmp(varargin{1,2}, 'PedalPressure')
                        F_brake=F_brake; %braking force should equal the pure force resultant from a pedal pressure
                        obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                    else
                        obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                    end
                    %if it is designating the propulsion input type, braking is
                    %assumed to be a pure force
                elseif strcmp(varargin{1,1}, 'PropulsionInputType')
                    if strcmp(varargin{1,2},'PureForce')
                        obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                    elseif strcmp(varargin{1,2}, 'EngineTorque')
                        if obj.N1==1 && obj.N2==1
                            warning('Looks like you never set your gear ratios, buddy');
                        end
                        B_eff=((obj.N1^2*obj.N2^2*obj.transB)/(obj.r_eff^2))+((obj.transB*obj.N2^2)/(obj.r_eff^2))+obj.transB/(obj.r_eff);
                        Floss=((obj.N1*obj.N2*obj.transLoss)/(obj.r_eff^2))+((obj.N2*obj.transLoss)/obj.r_eff^2)+(obj.transLoss/obj.r_eff);
                        obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop*(obj.N1*obj.N2)/(obj.r_eff))-(B_eff*obj.xdot(end))-(Floss);
                    else
                        
                    end
                end
                % if two force designations are given
            elseif size(varargin,2)==4
                
                %if the first set of designations are for brakes
                if strcmp(varargin{1,1},'BrakeInputType')
                    if strcmp(varargin{1,2},'PureForce')
                        if strcmp(varargin{1,4},'PureForce')
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        elseif strcmp(varargin{1,4}, 'EngineTorque')
                            if obj.N1==1 && obj.N2==1
                                warning(1,'Looks like you never set your gear ratios, buddy');
                            end
                            B_eff=((obj.N1^2*obj.N2^2*obj.transB)/(obj.r_eff^2))+((obj.transB*obj.N2^2)/(obj.r_eff^2))+obj.transB/(obj.r_eff);
                            Floss=((obj.N1*obj.N2*obj.transLoss)/(obj.r_eff^2))+((obj.N2*obj.transLoss)/obj.r_eff^2)+(obj.transLoss/obj.r_eff);
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop*(obj.N1*obj.N2)/(obj.r_eff))-(B_eff*obj.xdot(end))-(Floss);
                        else
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        end
                    elseif strcmp(varargin{1,2}, 'PedalPressure')
                        if strcmp(varargin{1,4},'PureForce')
                            F_brake=F_brake; %braking force should equal the pure force resultant from a pedal pressure
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        elseif strcmp(varargin{1,4}, 'EngineTorque')
                            if obj.N1==1 && obj.N2==1
                                warning(1,'Looks like you never set your gear ratios, buddy');
                            end
                            F_brake=F_brake; %braking force should equal the pure force resultant from a pedal pressure
                            B_eff=((obj.N1^2*obj.N2^2*obj.transB)/(obj.r_eff^2))+((obj.transB*obj.N2^2)/(obj.r_eff^2))+obj.transB/(obj.r_eff);
                            Floss=((obj.N1*obj.N2*obj.transLoss)/(obj.r_eff^2))+((obj.N2*obj.transLoss)/obj.r_eff^2)+(obj.transLoss/obj.r_eff);
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop*(obj.N1*obj.N2)/(obj.r_eff))-(B_eff*obj.xdot(end))-(Floss);
                        else
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        end
                    else
                        if strcmp(varargin{1,4},'PureForce')
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        elseif strcmp(varargin{1,4}, 'EngineTorque')
                            if obj.N1==1 && obj.N2==1
                                warning(1,'Looks like you never set your gear ratios, buddy');
                            end
                            B_eff=((obj.N1^2*obj.N2^2*obj.transB)/(obj.r_eff^2))+((obj.transB*obj.N2^2)/(obj.r_eff^2))+obj.transB/(obj.r_eff);
                            Floss=((obj.N1*obj.N2*obj.transLoss)/(obj.r_eff^2))+((obj.N2*obj.transLoss)/obj.r_eff^2)+(obj.transLoss/obj.r_eff);
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop*(obj.N1*obj.N2)/(obj.r_eff))-(B_eff*obj.xdot(end))-(Floss);
                        else
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        end
                    end
                    
                    %if the first set of designations are for propulsion
                elseif strcmp(varargin{1,1}, 'PropulsionInputType')
                    if strcmp(varargin{1,2},'PureForce')
                        if strcmp(varargin{1,4},'PureForce')
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        elseif strcmp(varargin{1,4}, 'PedalPressure')
                            F_brake=F_brake; %braking force should equal the pure force resultant from a pedal pressure
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        else
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        end
                    elseif strcmp(varargin{1,2}, 'EngineTorque')
                        if strcmp(varargin{1,4},'PureForce')
                            if obj.N1==1 && obj.N2==1
                                warning(1,'Looks like you never set your gear ratios, buddy');
                            end
                            B_eff=((obj.N1^2*obj.N2^2*obj.transB)/(obj.r_eff^2))+((obj.transB*obj.N2^2)/(obj.r_eff^2))+obj.transB/(obj.r_eff);
                            Floss=((obj.N1*obj.N2*obj.transLoss)/(obj.r_eff^2))+((obj.N2*obj.transLoss)/obj.r_eff^2)+(obj.transLoss/obj.r_eff);
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop*(obj.N1*obj.N2)/(obj.r_eff))-(B_eff*obj.xdot(end))-(Floss);
                        elseif strcmp(varargin{1,4}, 'PedalPressure')
                            if obj.N1==1 && obj.N2==1
                                warning(1,'Looks like you never set your gear ratios, buddy');
                            end
                            F_brake=F_brake;
                            B_eff=((obj.N1^2*obj.N2^2*obj.transB)/(obj.r_eff^2))+((obj.transB*obj.N2^2)/(obj.r_eff^2))+obj.transB/(obj.r_eff);
                            Floss=((obj.N1*obj.N2*obj.transLoss)/(obj.r_eff^2))+((obj.N2*obj.transLoss)/obj.r_eff^2)+(obj.transLoss/obj.r_eff);
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop*(obj.N1*obj.N2)/(obj.r_eff))-(B_eff*obj.xdot(end))-(Floss);
                        else
                            B_eff=((obj.N1^2*obj.N2^2*obj.transB)/(obj.r_eff^2))+((obj.transB*obj.N2^2)/(obj.r_eff^2))+obj.transB/(obj.r_eff);
                            Floss=((obj.N1*obj.N2*obj.transLoss)/(obj.r_eff^2))+((obj.N2*obj.transLoss)/obj.r_eff^2)+(obj.transLoss/obj.r_eff);
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop*(obj.N1*obj.N2)/(obj.r_eff))-(B_eff*obj.xdot(end))-(Floss);
                        end
                    else
                        if strcmp(varargin{1,4},'PureForce')
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        elseif strcmp(varargin{1,4}, 'PedalPressure')
                            F_brake=F_brake; %braking force should equal the pure force resultant from a pedal pressure
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        else
                            obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
                        end
                    end
                end
                %if nothing, assume 2 pure forces
            else
                obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^2))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr)+(F_prop/obj.mass);
            end
            
            
            %integrating to get the other x values
            obj.xdot(end+1)=obj.xdot(end)+(.5*time_step*(obj.xdoubledot(end)+obj.xdoubledot(end-1)));
            obj.x(end+1)=obj.x(end)+(.5*time_step*(obj.xdot(end)+obj.xdot(end-1)));
            
            %y values, found from given constant y dot, needs to be derived
            %properly from ydoubledot
            obj.ydoubledot(end+1)=0;
            obj.ydot(end+1)=obj.ydotnaught;
            obj.y(end+1)=obj.y(end)+(.5*time_step*(obj.ydot(end)+obj.ydot(end-1)));
            
            %integrated from given yaw rate, not really calculated in any
            %real way.
            obj.yaw(end+1)=obj.yaw(end)+(.5*time_step*(obj.yawdot+obj.yawdot));
            
            %logging global coordinates
            [vn2, ve2]=obj.getGlobalVelocity;
            obj.N(end+1)=obj.N(end)+(time_step*.5*(vn1+vn2));
            obj.E(end+1)=obj.E(end)+(time_step*.5*(ve1+ve2));
            if size(obj.terrain,1)==1 && size(obj.terrain,2)==1
                obj.D(end+1)=obj.D(end)+(sin(obj.terrain)*(obj.x(end)-obj.x(end-1)));
            else
                obj.D(end+1)=obj.terrain(ceil(obj.N(end)+(size(obj.terrain,1)/2)+.5), ceil(obj.E(end)+(size(obj.terrain,2)/2)+.5));
            end
            
            %return
            out=1;
        end
        %% Calculate Tire Velocities
        %returns car velocity in car tire frames, assuming low speed
        %turning criteria
        function [rlx,rrx,flx,frx,rly,rry,fly,fry]=calcTireVelocities(obj)
            cvx=obj.xdot(end);
            cvy=obj.ydot(end);
            r=sqrt((obj.width/2)^2+(obj.length/2)^2);
            
            rlx=cvx+obj.yawdot*r*sin(atan(obj.length/obj.width));
            rrx=cvx+obj.yawdot*r*sin(atan(obj.length/obj.width));
            flx=cvx*cos(obj.del)+cvy*sin(obj.del)+obj.yawdot*r*sin((pi/2)-atan(obj.length/obj.width)-obj.del);
            frx=cvx*sin(obj.del)+cvy*cos(obj.del)+obj.yawdot*r*sin((pi/2)-atan(obj.length/obj.width)-obj.del);
            
            rly=cvy+obj.yawdot*r*cos(atan(obj.length/obj.width));
            rry=cvy+obj.yawdot*r*cos(atan(obj.length/obj.width));
            fly=cvx*cos(obj.del)+cvy*sin(obj.del)+obj.yawdot*r*cos((pi/2)-atan(obj.length/obj.width)-obj.del);
            fry=cvx*sin(obj.del)+cvy*cos(obj.del)+obj.yawdot*r*cos((pi/2)-atan(obj.length/obj.width)-obj.del);
        end
        %% Get Global Velocity
        %returns the velocity of COM in global coordinates.  Assumes
        %car starts at [0,0,0] does not yet support down dimension.  Does
        %not yet support Down dimension
        function [vn, ve]=getGlobalVelocity(obj)
            vn=obj.xdot(end)*cos(obj.yaw(end))+obj.ydot(end)*cos(pi/4-obj.yaw(end));
            ve=obj.xdot(end)*sin(obj.yaw(end))+obj.ydot(end)*sin(pi/4-obj.yaw(end));
        end
        %% Simulate Coasting Models
        %simulates and returns the output from a coasting model for a
        %quarter of a Km.  Uses designated aerodynamic drag model
        %order.
        function out=simulateCoastingModels(obj, time_step,order)
            timevec=0;
            obj.resetValues;
            while obj.x(end)<250
                obj.stepSim(time_step,order,0);
                timevec(end+1)=timevec(end)+time_step;
            end
            out=horzcat(timevec',obj.x', obj.xdot');
        end
        %% Set C_d
        %Helper function to set the drag coefficient.
        function out=setDrags(obj,varargin)
            if isempty(varargin)
                obj.Cd=0;
                obj.C_rr=0;
                obj.fa=0;
            else
                obj.Cd=varargin{1,1};
                obj.C_rr=varargin{1,2};
                obj.fa=varargin{1,3};
            end
            out=[obj.Cd,obj.C_rr,obj.fa];
        end
        %% Reset Values
        %resets the simulation vectors to initial state
        function out=resetValues(obj)
            obj.x=0;
            obj.xdot=obj.xdotnaught;
            obj.xdoubledot=0;
            obj.y=0;
            obj.ydot=obj.ydotnaught;
            obj.ydoubledot=0;
            obj.yaw=0;
            obj.yawdot=obj.yawdotnaught;
            obj.N=0;
            obj.E=0;
            obj.D=0;
            out=1;
        end
        %% Get State
        %returns relevant current state parameters
        function [h, xo,xdo,xddo,yo,ydo,yddo]=getState(obj)
            h=obj.yaw(end);
            xo=obj.x(end);
            xdo=obj.xdot(end);
            xddo=obj.xdoubledot(end);
            yo=obj.y(end);
            ydo=obj.ydot(end);
            yddo=obj.ydoubledot(end);
        end
        %% Get Global State
        %returns state in global frame
        function [vn,ve,n,e,d]=getGlobalState(obj)
            [vn, ve]=obj.getGlobalVelocity;
            n=obj.N(end);
            e=obj.E(end);
            d=obj.D(end);
        end
        %% setPlanarLandscape
        %imports a matrix, where rows correspond to East, in meters
        %columns correspond to North in meters, and the values of the
        %indecies are Down, in meters.  Num rows and num cols must be odd,
        %as the center point is the dimension (0,0,0) in NED global
        %coorinates.  Function normalizes the Down values to the center
        %being 0.
        
        %if the size of the matrix is one, the value is assumed to be a
        %localized gradient in radians
        function out=setPlanarLandscape(obj, inMat)
            if(mod(size(inMat, 1),2)~=1)
                out=0;
                error('matrix must have an odd number of rows');
            elseif (mod(size(inMat, 2),2)~=1)
                out=0;
                error('matrix must have an odd number of columns');
            else
                obj.terrain=inMat;
                out=1;
            end
            if size(obj.terrain,1)==1 && size(obj.terrain,2)==1
                obj.D=0;
            else
                obj.D=inMat(floor(size(inMat,1)/2), floor(size(inMat,2)/2));
            end
        end
        %% Get Local Grandient
        %uses the global position and velocity to determine current
        %location on the terrain map, as well as the likely next position,
        %and exrapolates a local gradient from that.
        function out=getLocalGradient(obj)
            if size(obj.terrain,1)==1 && size(obj.terrain,2)==1
                out=obj.terrain;
            else
                [vn,ve,n,e,~]=obj.getGlobalState;
                n=n+((size(obj.terrain,1)/2)+.5);
                e=e+((size(obj.terrain,2)/2)+.5);
                
                curr_elevation=obj.terrain(ceil(n), ceil(e));
                theta=atan(vn/ve);
                if (theta>=0 && theta<=pi/8) || (theta<=0 && theta>=(15*pi)/8)
                    next_elevation=obj.terrain(ceil(n), ceil(e)+1);
                elseif (theta>pi/8) && (theta<=(3*pi)/8)
                    next_elevation=obj.terrain(ceil(n)+1, ceil(e)+1);
                elseif (theta>(3*pi)/8) && (theta<=(5*pi)/8)
                    next_elevation=obj.terrain(ceil(n)+1, ceil(e));
                elseif (theta>(5*pi)/8) && (theta<=(7*pi)/8)
                    next_elevation=obj.terrain(ceil(n)+1, ceil(e)-1);
                elseif (theta>(7*pi)/8) && (theta<=(9*pi)/8)
                    next_elevation=obj.terrain(ceil(n), ceil(e)-1);
                elseif (theta>(9*pi)/8) && (theta<=(11*pi)/8)
                    next_elevation=obj.terrain(ceil(n)-1, ceil(e)-1);
                elseif (theta>(11*pi)/8) && (theta<=(13*pi)/8)
                    next_elevation=obj.terrain(ceil(n)-1, ceil(e));
                else
                    next_elevation=obj.terrain(ceil(n)-1, ceil(e)+1);
                end
                %if each entry is one meter;
                dist=1;
                out=atan((next_elevation-curr_elevation)/(dist));
            end
        end
        %% Set Trans Values
        %allows the user to set various gear ratios in the transmission and
        %differential.  Also set the friction and visous losses, which are
        %assumed to be zero if not included.
        %When the gear ratios are updated, transmission losses and friction
        %are not changed unless noted.
        function out=setTransValues(obj, n1, n2, varargin)
            obj.N1=n1;
            obj.N2=n2;
            if size(varargin,2)==2
                if strcmp(varargin{1,1},'EffectiveDamping')
                    obj.transB=varargin{1,2};
                elseif strcmp(varargin{1,1},'FrictionLoss')
                    obj.transLoss=varargin{1,2};
                end
            elseif size(varargin,2)==4
                if strcmp(varargin{1,1},'EffectiveDamping')
                    obj.transB=varargin{1,2};
                    obj.transLoss=varargin{1,4};
                elseif strcmp(varargin{1,1},'FrictionLoss')
                    obj.transLoss=varargin{1,2};
                    obj.transB=varargin{1,4};
                end
            end
            out=1;
        end
    end
    
end