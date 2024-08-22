% % Define the wavelength range and sampling interval
% lambda_min = 1500;  % nm
% lambda_max = 1600;  % nm
% d_lambda = 0.1;     % nm
% 
% % Generate the wavelength vector
% lambda = lambda_min:d_lambda:lambda_max;
% 
% % Define the refractive index of the substrate and the thin films
% n_substrate = 1.5;  % silicon
% n_film1 = 1.8;      % thin film 1
% n_film2 = 1.6;      % thin film 2
% 
% % Define the thicknesses of the substrate and the thin films
% t_substrate = 1000; % nm
% t_film1 = 200;      % nm
% t_film2 = 300;      % nm
% 
% % Calculate the optical path length difference between the reference arm and the sample arm
% delta_opd = 2 * (n_substrate * t_substrate + n_film1 * t_film1 + n_film2 * t_film2);
% 
% % Generate the interference spectrum
% interference_spectrum = cos(2*pi*delta_opd./(lambda./n_substrate));
% 
% % Plot the interference fringes
% figure;
% plot(lambda, interference_spectrum);
% xlabel('Wavelength (nm)');
% ylabel('Interference Intensity');
% title('Interference Fringes');
% Define parameters
n = 3; % Number of layers
lambda0 = 10.6e-6; % Central wavelength of optical comb (infrared region)
deltaLambda = 1e-10; % Spectral resolution of optical comb
L = 1; % Length of interferometer arm
theta = 45; % Angle of incidence
m = 1000; % Number of points in wavelength axis
t = zeros(1, n); % Thickness of each layer
for i = 1:n
    t(i) = (i-1) * lambda0 / (4 * sqrt(2) * cosd(theta)); % Calculate thickness of each layer
end

% Generate optical comb
frep = 80e6; % Repetition rate of comb
k = 2 * pi / lambda0; % Wave vector
omega0 = k * 3e8; % Central frequency of comb
deltaOmega = 2 * pi * frep; % Frequency spacing of comb
omega = linspace(omega0 - m/2*deltaOmega, omega0 + m/2*deltaOmega, m); % Frequency axis
lambda = 2 * pi * 3e8 ./ omega; % Wavelength axis
comb = exp(1i * deltaOmega * L / 3e8 * (cosd(theta) * sqrt(k^2 - (omega - omega0).^2/c^2) - k)); % Optical comb

% Calculate reflectance and phase shift for each layer
R = zeros(1, n);
phi = zeros(1, n);
figure;
for i = 1:n
    delta = 2 * pi * t(i) * sqrt(k^2 - (omega - omega0).^2/c^2) / lambda0; % Phase shift
    phi(i) = sum(delta); % Total phase shift
    R(i) = abs(sum(comb .* exp(1i * delta))); % Reflectance
    subplot(n, 1, i);
    plot(lambda*1e6, abs(R(i)*comb).^2, 'b-', 'LineWidth', 1); % Display spectrum
    xlabel('Wavelength (\mum)');
    ylabel('Intensity');
    title(sprintf('Layer %d Reflectance Spectrum', i));
end

% Calculate total physical thickness
deltaPhi = phi(n) - phi(1); % Total phase shift across all layers
d = deltaPhi * lambda0 / (4 * pi * sqrt(2) * cosd(theta)); % Total physical thickness

% Display results
fprintf('Total physical thickness: %.3f microns\n', d*1e6);
fprintf('Thickness of each layer:\n');
for i = 1:n
    fprintf('Layer %d: %.3f microns\n', i, t(i)*1e6);
end
