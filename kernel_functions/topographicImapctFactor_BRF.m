function [ue,use,ge,T] = topographicImapctFactor_BRF(slope,aspect,solarZenith,solarAzimuth,viewZenith,viewAzimuth,shadowBinaryFctorofSolar,shadowBinaryFctorofView)
% generate topographic impact factor under clear sky
% the meaning of each input variable is obvious!
% output: ue is cos(VZA); use is cos(SZA); ge is cos(RAA); T is TopographicImapctFactor
% written by dalei
% 2017-10-10

%% initialize all variables
%s sunlit; v visible
count_i = 33;
count_j = 33;
Asv = 0;
Av_uj=0;
Asv_uj=0;
Asv_usj_us = 0;
Asv_vsj_vs=0;
Asv_gj_g=0;
for i = 1:count_i
    for j = 1:count_j
        
        relativeSolarZenith =  cos(slope(i,j))*cos(solarZenith)+sin(slope(i,j))*sin(solarZenith)*cos(aspect(i,j) -solarAzimuth );
        relativeViewZenith =  cos(slope(i,j))*cos(viewZenith)+sin(slope(i,j))*sin(viewZenith)*cos(aspect(i,j) -viewAzimuth );
        relativeSolarAzimuth = relativeAzimuthCal(solarZenith,solarAzimuth,slope(i,j),aspect(i,j));
        relativeViewAzimuth = relativeAzimuthCal(viewZenith,viewAzimuth,slope(i,j),aspect(i,j));
        
        % limit
        if(relativeSolarZenith <0)
            shadowBinaryFctorofSolar(i,j) = 0;
            relativeSolarZenith = 0;
        end
        if(relativeViewZenith <0)
            shadowBinaryFctorofView(i,j) = 0;
            relativeViewZenith = 0;
        end
        if(relativeSolarZenith>1)
            relativeSolarZenith = 1;
        end
        if(relativeViewZenith>1)
            relativeViewZenith = 1;
        end
        
        % cumulate
        Asv = Asv+shadowBinaryFctorofSolar(i,j)*shadowBinaryFctorofView(i,j)/ cos(slope(i,j));
        %Asv_uj_u =Asv_uj_u+(relativeViewZenith-cos(viewZenith))*shadowBinaryFctorofView(i,j)*shadowBinaryFctorofSolar(i,j)/ cos(slope(i,j));
        Av_uj=Av_uj+shadowBinaryFctorofView(i,j)*relativeViewZenith/ cos(slope(i,j));
        Asv_uj=Asv_uj+shadowBinaryFctorofSolar(i,j)*shadowBinaryFctorofView(i,j)*relativeViewZenith/ cos(slope(i,j));
        Asv_usj_us = Asv_usj_us+(relativeSolarZenith-cos(solarZenith))*shadowBinaryFctorofSolar(i,j)*shadowBinaryFctorofView(i,j)/ cos(slope(i,j));
        Asv_vsj_vs=Asv_vsj_vs+(cos(relativeSolarAzimuth)-cos(solarAzimuth))*shadowBinaryFctorofSolar(i,j)*shadowBinaryFctorofView(i,j)/ cos(slope(i,j));
        Asv_gj_g=Asv_gj_g+(cos(relativeViewAzimuth-relativeSolarAzimuth)-cos(viewAzimuth-solarAzimuth))*shadowBinaryFctorofSolar(i,j)*shadowBinaryFctorofView(i,j)/ cos(slope(i,j));
    end
end

ue = Asv_uj/Asv;
use = Asv_usj_us/Asv+cos(solarZenith);

% range limit
if(use>1)
    use=1;
end
if(use<0)
    use = 0;
end
if(ue>1)
    ue=1;
end
if(ue<0)
    ue = 0;
end

ge = (cos(solarZenith)*cos(viewZenith)+ sin(solarZenith)*sin(viewZenith)*cos(solarAzimuth-viewAzimuth)-ue*use)/(sin(acos(ue))*sin(acos(use)));

if(ge>1)
    ge=1;
end
if(ge<-1)
    ge = -1;
end

T = use*ue/cos(solarZenith)* Asv/Av_uj;
end