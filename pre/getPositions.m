function pos = getPositions(Resolution,radiusMultiplier, angles)
    Center = [Resolution(2) Resolution(1)]/2;
    radius = min(Center);
    radius = (radius - 0.15*radius)*radiusMultiplier;
    if isempty(angles)
        angles = [50 170 290];
    end
        
    pos = rand(3,2);
    for i = 1:length(angles)
        pos(i,:) = [Center(1)+radius*cosd(angles(i)) Center(2)+radius*sind(angles(i))];
    end
end