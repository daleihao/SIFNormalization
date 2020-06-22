function relativeAzimuth = relativeAzimuthCal(solarZenith,solarAzimuth,slope,aspect)
x=sin(solarZenith).*cos(slope).*cos(solarAzimuth-aspect)-cos(solarZenith).*sin(slope);
y=sin(solarAzimuth-aspect).*sin(solarZenith);
relativeAzimuth=atan2(y,x);
if relativeAzimuth<0
    relativeAzimuth = relativeAzimuth + 2*pi;
end
if relativeAzimuth>2*pi
    relativeAzimuth = relativeAzimuth - 2*pi;
end

