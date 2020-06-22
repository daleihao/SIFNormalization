function [ue,T2] = TopographicImapctFactor_HDRF(slope,aspect,viewZenith,viewAzimuth,shadowBinaryFctorofView,skyViewFactor)
% generate topographic impact factor under overcast sky
% the meaning of each input variable is obvious!
% output: ue is cos(VZA); T2 is TopographicImapctFactor
% written by dalei
% 2017-10-10

%% initialize all variables
%v visible
count_i = 33;
count_j = 33;
Av_V=0;
Av_uj_u_V =0;
Av_uj=0;
Av_uj_V = 0;

%% calculate TopographicImapctFactor
for i = 1:count_i
    for j = 1:count_j
        relativeViewZenith =  cos(slope(i,j))*cos(viewZenith)+sin(slope(i,j))*sin(viewZenith)*cos(aspect(i,j) -viewAzimuth );
        
        % limit
        if(relativeViewZenith <0)
            shadowBinaryFctorofView(i,j) = 0;
            relativeViewZenith = 1;
        end
        if(relativeViewZenith>1)
            relativeViewZenith = 1;
        end
        
        % cumulate
        Av_V= Av_V+shadowBinaryFctorofView(i,j)*skyViewFactor(i,j)/ cos(slope(i,j));
        Av_uj_u_V=Av_uj_u_V+(relativeViewZenith-cos(viewZenith))*shadowBinaryFctorofView(i,j)*skyViewFactor(i,j)/ cos(slope(i,j));
        Av_uj=Av_uj+shadowBinaryFctorofView(i,j)*relativeViewZenith/ cos(slope(i,j));
        Av_uj_V = Av_uj_V+relativeViewZenith*shadowBinaryFctorofView(i,j)*skyViewFactor(i,j)/ cos(slope(i,j));
    end
end

ue = Av_uj_V/Av_V;

% range limit
if(ue>1)
    ue=1;
end
if(ue<0)
    ue = 0;
end

T2 = ue*Av_V/Av_uj;

end