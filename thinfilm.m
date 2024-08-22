% Define the wavelengths of light
wavelengths = linspace(400, 700, 300);

% Define the incident angles (in degrees)
angles = linspace(0, 90, 91);

% Define the refractive indices of the ambient medium and the substrate
n_amb = 1; % air
n_sub = 3.5; % substrate

% Define the refractive indices and thicknesses of the thin film layers
n_layers = [1.62, 1.8]; % refractive indices
d_layers = [0,140]; % thicknesses in nm

% Initialize arrays to store the reflectance data
reflectance_data = zeros(length(wavelengths), length(angles));

% Loop over the incident angles
for i = 1:length(angles)
    theta = angles(i); % current incident angle
    for j = 1:length(wavelengths)
        lambda = wavelengths(j); % current wavelength
        
        % Calculate the reflection and transmission coefficients for each layer
        r_coeff = (n_amb.*cosd(theta) - n_layers.*sqrt(1 - (n_amb./n_layers).^2.*sind(theta).^2)) ./ ...
            (n_amb.*cosd(theta) + n_layers.*sqrt(1 - (n_amb./n_layers).^2.*sind(theta).^2));
        t_coeff = 2*n_amb.*cosd(theta) ./ (n_amb.*cosd(theta) + n_layers.*sqrt(1 - (n_amb./n_layers).^2.*sind(theta).^2));
        
        % Calculate the phase shift for each layer
        delta = 2*pi*n_layers.*d_layers.*cosd(theta) ./ lambda;
        
        % Calculate the reflectance for the entire structure
        reflectance_data(j, i) = abs(sum(r_coeff.*exp(1i*delta).*t_coeff.^2)).^2;
    end
end

% Plot the reflectance graph
figure;
imagesc(angles, wavelengths, reflectance_data);
colorbar;
xlabel('Incident Angle (degrees)');
ylabel('Wavelength (nm)');
title('Angular Spectroscopic Reflectometry');
axis xy;

% Find the thicknesses of the thin film layers from the reflectance graph
thicknesses = zeros(1, length(n_layers));
for i = 1:length(n_layers)
    [~, index] = max(reflectance_data(:, i)); % find the index of the maximum reflectance
    thicknesses(i) = d_layers(i) * index / length(wavelengths); % calculate the thickness
end

% Display the thicknesses of the thin film layers
disp('Thicknesses of Thin Film Layers (nm):');
disp(thicknesses);
