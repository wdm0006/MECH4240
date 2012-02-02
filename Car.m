%% Car class 
% Car class represents a vehicle rolling on an arbitrary surface.  Initial 
% conditions and inputs allow various simulations, and get functions allow 
% the monitoring and analysis of the system's response. 
% Version 2
% -Added in a braking force input to simulation, negative is propulsive
% -Deprecated coasting sim, just use plain sim with no input
% -Added in terrain map ability, road grade integrated into that,
% deprecated as a seperate parameter
% -Seperated Constructor and simulation setup
% -Terminal Velocity property deprecated, user inputs Cd and Af

classdef Car < handle
    
%% Class Properties 
% System Properties: the parameters which define the car itself.
% Parameters for the simulation of the system
    properties 
        mass
        g=9.81;
        Cd
        width
        length 
        r_eff
        C_rr
        cg_height
        fa
        
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
%         Sets parameters of the car itself
%         simulation parameters set to 0 just in case
        function obj=Car(m,wbw, wbl,rf,crr,cgh, Af,d)
            obj.mass=m;
            obj.del=d;
            obj.length=wbl;
            obj.width=wbw;
            obj.r_eff=rf;
            obj.C_rr=crr;
            obj.cg_height=cgh;
            obj.fa=Af;
            obj.Cd=0;
            
            
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
%         sets the initial conditions of the simulation
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
%         Functionally the same as the coasting simulation, but with a
%         parameter for force applied by braking.  Because of the simplicity
%         of the model at this point, a negative braking force input would
%         simulate propulsion.  Does support rolling resistance, does not
%         support yaw double dot.
%         negative road grade is downhill.
        function out=stepSim(obj,time_step,order,F_brake)
                [vn1, ve1]=obj.getGlobalVelocity;
                gradient=obj.getLocalGradient;
                if order==-1
                    obj.xdoubledot(end+1)=-obj.g*sin(gradient)-(F_brake/obj.mass);
                    obj.xdot(end+1)=obj.xdot(end)+(.5*time_step*(obj.xdoubledot(end)+obj.xdoubledot(end-1)));
                    obj.x(end+1)=obj.x(end)+(.5*time_step*(obj.xdot(end)+obj.xdot(end-1)));
                    obj.ydoubledot(end+1)=0;
                    obj.ydot(end+1)=obj.ydot(end);
                    obj.y(end+1)=obj.y(end)+(.5*time_step*(obj.ydot(end)+obj.ydot(end)));
                    obj.yaw(end+1)=obj.yaw(end)+(.5*time_step*(obj.yawdot+obj.yawdot));
                else
                    obj.xdoubledot(end+1)=-obj.g*sin(gradient)-((obj.Cd*obj.fa*(obj.xdot(end)^order))/obj.mass)-(F_brake/obj.mass)-(obj.g*obj.C_rr);
                    obj.xdot(end+1)=obj.xdot(end)+(.5*time_step*(obj.xdoubledot(end)+obj.xdoubledot(end-1)));
                    obj.x(end+1)=obj.x(end)+(.5*time_step*(obj.xdot(end)+obj.xdot(end-1)));
                    obj.ydoubledot(end+1)=0;
                    obj.ydot(end+1)=obj.ydotnaught;
                    obj.y(end+1)=obj.y(end)+(.5*time_step*(obj.ydot(end)+obj.ydot(end-1)));
                    obj.yaw(end+1)=obj.yaw(end)+(.5*time_step*(obj.yawdot+obj.yawdot));
                end
                [vn2, ve2]=obj.getGlobalVelocity;
                obj.N(end+1)=obj.N(end)+(time_step*.5*(vn1+vn2));
                obj.E(end+1)=obj.E(end)+(time_step*.5*(ve1+ve2));
                if size(obj.terrain,1)==1 && size(obj.terrain,2)==1
                    obj.D(end+1)=obj.D(end)-(sin(obj.terrain)*(obj.x(end)-obj.x(end-1)));
                else
                    %obj.D(end+1)=obj.terrain(ceil(obj.N(end)+(size(obj.terrain,1)/2)+.5), ceil(obj.E(end)+(size(obj.terrain,2)/2)+.5));
                    obj.D(end+1)=obj.D(end)-(sin(gradient)*(obj.x(end)-obj.x(end-1)));
                end
                out=1;
        end
        %% Calculate Tire Velocities
%         returns car velocity in car tire frames, assuming low speed
%         turning criteria
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
%         returns the velocity of COM in global coordinates.  Assumes
%         car starts at [0,0,0] does not yet support down dimension.  Does
%         not yet support Down dimension
        function [vn, ve]=getGlobalVelocity(obj)
            vn=obj.xdot(end)*cos(obj.yaw(end))+obj.ydot(end)*cos(pi/4-obj.yaw(end));
            ve=obj.xdot(end)*sin(obj.yaw(end))+obj.ydot(end)*sin(pi/4-obj.yaw(end));
        end
        %% Simulate Coasting Models
%         simulates and returns the output from a coasting model for a
%         quarter of a Km.  Uses designated aerodynamic drag model
%         order.
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
%         Helper function to set the drag coefficient.
        function out=setCD(obj,varargin)
            if isempty(varargin)
                obj.Cd=0;
            else
                obj.Cd=varargin{1,1};
            end
            out=obj.Cd;
        end
        %% Reset Values
%         resets the simulation vectors to initial state
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
%         returns relevant current state parameters
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
%         returns state in global frame
        function [vn,ve,n,e,d]=getGlobalState(obj)
           [vn, ve]=obj.getGlobalVelocity;
           n=obj.N(end);
           e=obj.E(end);
           d=obj.D(end);
        end
        %% setPlanarLandscape
%         imports a matrix, where rows correspond to East, in meters
%         columns correspond to North in meters, and the values of the
%         indecies are Down, in meters.  Num rows and num cols must be odd,
%         as the center point is the dimension (0,0,0) in NED global
%         coorinates.  Function normalizes the Down values to the center
%         being 0.
%         if the size of the matrix is one, the value is assumed to be a
%         localized gradient in radians
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
%         uses the global position and velocity to determine current
%         location on the terrain map, as well as the likely next position,
%         and exrapolates a local gradient from that.
        function out=getLocalGradient(obj)
            if size(obj.terrain,1)==1 && size(obj.terrain,2)==1
                out=1*obj.terrain; %invert gradient, in accordance with Down coordiante axis
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
    end
    
end