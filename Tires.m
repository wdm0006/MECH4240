classdef Tires < handle
    properties
        a0
        a1
        a2
        a3
        a4
        a5
        a6
        a7
        a8
        k_phi_f
        k_phi_r
        w
        a
        b
        t_f
        t_r
        h1
        h_cg
    end
    
    methods
        function obj=Tires(k_phi_fi, k_phi_ri, wi, ai, bi, t_fi, t_ri, h1i,h_cgi)
            obj.a0=0;
            obj.a1=-22.1;
            obj.a2=1011;
            obj.a3=1078;
            obj.a4=1.82;
            obj.a5=.208;
            obj.a6=0;
            obj.a7=-.354;
            obj.a8=.707;
            if nargin~=0
                obj.k_phi_f=k_phi_fi;
                obj.k_phi_r=k_phi_ri;
                obj.w=wi;
                obj.a=ai;
                obj.b=bi;
                obj.t_f=t_fi;
                obj.t_r=t_ri;
                obj.h1=h1i;
                obj.h_cg=h_cgi;
            else %use default values
                nmdeg=57.3;
                obj.k_phi_f=480*nmdeg;
                obj.k_phi_r=250*nmdeg;
                obj.w=1475*9.81;
                obj.a=1.5;
                obj.b=1.5;
                obj.t_f=1.44;
                obj.t_r=1.44;
                obj.h1=.1;
                obj.h_cg=.6;
            end
        end
        function alpha=calc_alpha(obj,fzin, fzon,fydesn)
            [alpha,precision]=fzero(@(alpha) calc_fy_des(alpha, fzin, fzon,fydesn), 0, [0,1,1,1]);  %#ok<NASGU>
        end
        function fy=calc_fy_des(obj,alpha, fzi, fzo,fydes)
            fzi=fzi/1000;
            fzo=fzo/1000;
            C=1.3;
            Do=obj.a0.*fzo.*fzo.*fzo+obj.a1.*fzo.*fzo+obj.a2.*fzo;
            Bo=(obj.a3.*sin(obj.a4.*atan(obj.a5.*fzo)))./(C.*Do);
            Eo=obj.a6.*fzo.*fzo+obj.a7.*fzo+obj.a8;
            
            Di=obj.a0.*fzi.*fzi.*fzi+obj.a1.*fzi.*fzi+obj.a2.*fzi;
            Bi=(obj.a3.*sin(obj.a4.*atan(obj.a5.*fzi)))./(C.*Di);
            Ei=obj.a6.*fzi.*fzi+obj.a7.*fzi+obj.a8;
            
            phio=(1-Eo).*alpha+(Eo./Bo).*atan(Bo.*alpha);
            phii=(1-Ei).*alpha+(Ei./Bi).*atan(Bi.*alpha);
            
            fyo=Do.*sin(C.*atan(Bo.*phio));
            fyi=Di.*sin(C.*atan(Bi.*phii));
            if fzi==0
                fyi=0;
            end
            if fzo==0
                fyo=0;
            end
            fy=(fyo+fyi)-fydes;
        end
        function ay=calc_ay(obj,delin,aynaught, vx, vy,yawdot,mass) 
            L=obj.a+obj.b;
            w_sf=((obj.w*obj.b)/L); %NEWTONS
            w_sr=((obj.w*obj.a)/L);  %NEWTONS
            phi=((obj.w*obj.h1)/(obj.k_phi_f+obj.k_phi_r-obj.w*obj.h1))*aynaught; %RADIANS            
            del_fzf=(1/obj.t_f)*(obj.k_phi_f*(phi)+w_sf*obj.h_cg*aynaught); %NEWTONS
            del_fzr=(1/obj.t_r)*(obj.k_phi_r*(phi)+w_sr*obj.h_cg*aynaught); %NEWTONS
            fzof=max(((w_sf/2)+del_fzf),0);%NEWTONS
            fzif=max(((w_sf/2)-del_fzf),0);%NEWTONS
            fzor=max(((w_sr/2)+del_fzr),0);%NEWTONS
            fzir=max(((w_sr/2)-del_fzr),0);%NEWTONS
            [af,ar]=obj.calc_alphas(vy, yawdot, vx,delin);
            fyf=-1*obj.calc_fy(af, fzif, fzof);
            fyr=-1*obj.calc_fy(ar, fzir, fzor);
            ay=-vx*yawdot+fyr/mass+fyf/mass;
        end
        function del=calc_del_des(obj,ay,R,delin)
            %ay in g's
            
            L=obj.a+obj.b;
            w_sf=((obj.w*obj.b)/L); %NEWTONS
            w_sr=((obj.w*obj.a)/L);  %NEWTONS
            
            fyf=(((obj.w)*obj.b)/L)*ay; %NEWTONS
            fyr=(((obj.w)*obj.a)/L)*ay; %NEWTONS
            
            phi=((obj.w*obj.h1)/(obj.k_phi_f+obj.k_phi_r-obj.w*obj.h1))*ay; %RADIANS
            
            del_fzf=(1/obj.t_f)*(obj.k_phi_f*(phi)+w_sf*obj.h_cg*ay); %NEWTONS
            del_fzr=(1/obj.t_r)*(obj.k_phi_r*(phi)+w_sr*obj.h_cg*ay); %NEWTONS
            
            
            fzof=max(((w_sf/2)+del_fzf),0);%NEWTONS
            fzif=max(((w_sf/2)-del_fzf),0);%NEWTONS
            fzor=max(((w_sr/2)+del_fzr),0);%NEWTONS
            fzir=max(((w_sr/2)-del_fzr),0);%NEWTONS
            
            [~,fymaxf]=calc_max_fy(fzif, fzof);%NEWTONS
            [~,fymaxr]=calc_max_fy(fzir, fzor);%NEWTONS
            
            if(abs(fyf)>abs(fymaxf))
                error('Front end saturated: %9.2fm/s^2',ay*9.81);
            elseif(abs(fyr)>abs(fymaxr))
                error('rear end saturated: %9.2fm/s^2',ay*9.81);
            end
            
            alphaf=calc_alpha(fzif, fzof, fyf);%DEGREES
            alphar=calc_alpha(fzir, fzor, fyr);%DEGREES
            
            del=-asind(L/R)+(alphar-alphaf)-delin;
        end
        function del=calc_del(obj,ay,R)
            %ay in g's
            
            L=obj.a+obj.b;
            w_sf=((obj.w*obj.b)/L); %NEWTONS
            w_sr=((obj.w*obj.a)/L);  %NEWTONS
            
            fyf=(((obj.w)*obj.b)/L)*ay; %NEWTONS
            fyr=(((obj.w)*obj.a)/L)*ay; %NEWTONS
            
            phi=((obj.w*obj.h1)/(obj.k_phi_f+obj.k_phi_r-obj.w*obj.h1))*ay; %RADIANS
            
            del_fzf=(1/obj.t_f)*(obj.k_phi_f*(phi)+w_sf*obj.h_cg*ay); %NEWTONS
            del_fzr=(1/obj.t_r)*(obj.k_phi_r*(phi)+w_sr*obj.h_cg*ay); %NEWTONS
            
            
            fzof=max(((w_sf/2)+del_fzf),0);%NEWTONS
            fzif=max(((w_sf/2)-del_fzf),0);%NEWTONS
            fzor=max(((w_sr/2)+del_fzr),0);%NEWTONS
            fzir=max(((w_sr/2)-del_fzr),0);%NEWTONS
            
            [~,fymaxf]=calc_max_fy(fzif, fzof);%NEWTONS
            [~,fymaxr]=calc_max_fy(fzir, fzor);%NEWTONS
            
            if(abs(fyf)>abs(fymaxf))
                error('Front end saturated: %9.2fm/s^2',ay*9.81);
            elseif(abs(fyr)>abs(fymaxr))
                error('rear end saturated: %9.2fm/s^2',ay*9.81);
            end
            
            alphaf=calc_alpha(fzif, fzof, fyf);%DEGREES
            alphar=calc_alpha(fzir, fzor, fyr);%DEGREES
            
            del=-asind(L/R)+(alphar-alphaf);
        end
        %function to find max fy value on the tire curve for a certain tire model
        function [alpha,fyout]=calc_max_fy(obj,fzin, fzon)
            [alpha,fyout]=fminsearch(@(alpha) calc_fy(alpha, fzin, fzon), 0, [0,1,1]);
            fyout=abs(fyout);
            
        end
        %subfunction
        function fy=calc_fy(obj,alpha, fzi, fzo)
            fzi=fzi/1000;
            fzo=fzo/1000;
            C=1.3;
            Do=obj.a0.*fzo.*fzo.*fzo+obj.a1.*fzo.*fzo+obj.a2.*fzo;
            Bo=(obj.a3.*sin(obj.a4.*atan(obj.a5.*fzo)))./(C.*Do);
            Eo=obj.a6.*fzo.*fzo+obj.a7.*fzo+obj.a8;
            
            Di=obj.a0.*fzi.*fzi.*fzi+obj.a1.*fzi.*fzi+obj.a2.*fzi;
            Bi=(obj.a3.*sin(obj.a4.*atan(obj.a5.*fzi)))./(C.*Di);
            Ei=obj.a6.*fzi.*fzi+obj.a7.*fzi+obj.a8;
            
            phio=(1-Eo).*alpha+(Eo./Bo).*atan(Bo.*alpha);
            phii=(1-Ei).*alpha+(Ei./Bi).*atan(Bi.*alpha);
            
            fyo=Do.*sin(C.*atan(Bo.*phio));
            fyi=Di.*sin(C.*atan(Bi.*phii));            
            if fzi==0
                fyi=0;
            end
            if fzo==0
                fyo=0;
            end
            fy=-(fyo+fyi);
        end
        function [calphaf, calphar]=getCorneringStiffness(obj,ay)
            L=obj.a+obj.b;
            w_sf=((obj.w*obj.b)/L); %NEWTONS
            w_sr=((obj.w*obj.a)/L);  %NEWTONS
            
            phi=((obj.w*obj.h1)/(obj.k_phi_f+obj.k_phi_r-obj.w*obj.h1))*ay; %RADIANS
            
            del_fzf=(1/obj.t_f)*(obj.k_phi_f*(phi)+w_sf*obj.h_cg*ay); %NEWTONS
            del_fzr=(1/obj.t_r)*(obj.k_phi_r*(phi)+w_sr*obj.h_cg*ay); %NEWTONS
            
            
            fzof=max(((w_sf/2)+del_fzf),0);%NEWTONS
            fzif=max(((w_sf/2)-del_fzf),0);%NEWTONS
            fzor=max(((w_sr/2)+del_fzr),0);%NEWTONS
            fzir=max(((w_sr/2)-del_fzr),0);%NEWTONS
            
            fzif=fzif/1000;
            fzof=fzof/1000;
            fzir=fzir/1000;
            fzor=fzor/1000;
            
            Cf=1.3;
            Dof=obj.a0.*fzof.*fzof.*fzof+obj.a1.*fzof.*fzof+obj.a2.*fzof;
            Bof=(obj.a3.*sin(obj.a4.*atan(obj.a5.*fzof)))./(Cf.*Dof);
            
            Dif=obj.a0.*fzif.*fzif.*fzif+obj.a1.*fzif.*fzif+obj.a2.*fzif;
            Bif=(obj.a3.*sin(obj.a4.*atan(obj.a5.*fzif)))./(Cf.*Dif);
            
            Cr=1.3;
            Dor=obj.a0.*fzor.*fzor.*fzor+obj.a1.*fzor.*fzor+obj.a2.*fzor;
            Bor=(obj.a3.*sin(obj.a4.*atan(obj.a5.*fzor)))./(Cr.*Dor);
            
            Dir=obj.a0.*fzir.*fzir.*fzir+obj.a1.*fzir.*fzir+obj.a2.*fzir;
            Bir=(obj.a3.*sin(obj.a4.*atan(obj.a5.*fzir)))./(Cr.*Dir);
            
            calphaf=(Bof*Cf*Dof)+(Bif*Cf*Dif);
            calphar=(Bor*Cr*Dor)+(Bir*Cr*Dir);
        end
        function [yawdoubledot,alphaf, alphar]=calcYawRateDot(obj, delin, vx, vy, yawdot, mass, aynaught)
            Iz=(obj.t_f*(obj.a+obj.b)^3)/12;
            L=obj.a+obj.b;
            w_sf=((obj.w*obj.b)/L); %NEWTONS
            w_sr=((obj.w*obj.a)/L);  %NEWTONS
            phi=((obj.w*obj.h1)/(obj.k_phi_f+obj.k_phi_r-obj.w*obj.h1))*aynaught; %RADIANS            
            del_fzf=(1/obj.t_f)*(obj.k_phi_f*(phi)+w_sf*obj.h_cg*aynaught); %NEWTONS
            del_fzr=(1/obj.t_r)*(obj.k_phi_r*(phi)+w_sr*obj.h_cg*aynaught); %NEWTONS
            fzof=max(((w_sf/2)+del_fzf),0);%NEWTONS
            fzif=max(((w_sf/2)-del_fzf),0);%NEWTONS
            fzor=max(((w_sr/2)+del_fzr),0);%NEWTONS
            fzir=max(((w_sr/2)-del_fzr),0);%NEWTONS
            [af,ar]=obj.calc_alphas(vy, yawdot, vx,delin);
            fyf=-1*obj.calc_fy(af, fzif, fzof);
            fyr=-1*obj.calc_fy(ar, fzir, fzor);
            yawdoubledot=obj.a*fyf/Iz-obj.b*fyr/Iz;
            alphaf=(vy+obj.a*yawdot)/vx-delin*(180/pi);
            alphar=(vy-obj.b*yawdot)/vx;
            
        end
        function [alphaf, alphar]=calc_alphas(obj, vy, yawdot, vx,delin)
            alphaf=(vy+obj.a*yawdot)/vx-delin*(180/pi);
            alphar=(vy+obj.b*yawdot)/vx;
        end
        
    end    
    
end