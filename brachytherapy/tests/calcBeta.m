function beta = calcBeta(r, theta,L)
    % calculates beta from r[cm], theta [deg] and L[cm] 
    % array inputs are allowed for theta
   
    assert(0 <= min(theta,[],"all") & 180 >= max(theta,[],"all"), 'theta out of bounds')
    
    
    r1 = sqrt(r.^2 + (L/2)^2 - r.*L.*cosd(180 - theta)); % cosine theorem
    r2 = sqrt(r.^2 + (L/2)^2 - r.*L.*cosd(theta)); % cosine theorem
    
    beta1 = asind(sind(180-theta).*L/2./r1); % sine theorem
    beta2 = asind(sind(theta).*L/2./r2); % sine theorem
    
    beta = beta1 + beta2;
end