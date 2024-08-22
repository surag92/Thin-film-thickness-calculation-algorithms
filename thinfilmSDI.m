% Define the wavelength range and sampling interval
lambda_min = 1500;  % nm
lambda_max = 1600;  % nm
d_lambda = 0.1;     % nm

% Generate the wavelength vector
lambda = lambda_min:d_lambda:lambda_max;

% Define the refractive index of the substrate and the thin films
n_substrate = 1.5;  % silicon
n_film1 = 1.8;      % thin film 1
n_film2 = 1.6;      % thin film 2

% Define the thicknesses of the substrate and the thin films
t_substrate = 1000; % nm
t_film1 = 200;      % nm
t_film2 = 300;      % nm

% Calculate the optical path length difference between the reference arm and the sample arm
delta_opd = 2 * (n_substrate * t_substrate + n_film1 * t_film1 + n_film2 * t_film2);

% Generate the interference spectrum
interference_spectrum = cos(2*pi*delta_opd./(lambda./n_substrate));

% Plot the interference spectrum
figure;
plot(lambda, interference_spectrum);
xlabel('Wavelength (nm)');
ylabel('Interference Intensity');
title('Interference Spectrum');

% Calculate the thicknesses of the thin films using the interference fringe spacing
fringe_spacing = (lambda./2)./(n_substrate * (t_film1 + t_film2));
t_film1_calc = (fringe_spacing/2) - t_film2;
t_film2_calc = (fringe_spacing/2) - t_film1;

% Calculate the substrate thickness using the phase shift method with NIR illumination
lambda_nir = 1550;        % nm
n_substrate_nir = 3.4;    % silicon (at NIR wavelength)
delta_opd_nir = (2*pi/lambda_nir) * (n_substrate_nir * t_substrate);
phase_shift = angle(cos(2*pi*delta_opd_nir/(lambda_nir/n_substrate_nir)) + 1i*sin(2*pi*delta_opd_nir/(lambda_nir/n_substrate_nir)));
t_substrate_calc = (phase_shift*lambda_nir)/(4*pi*n_substrate_nir);

% Display the calculated thicknesses
disp(['Thin Film 1 Thickness: ', num2str(t_film1_calc), ' nm']);
disp(['Thin Film 2 Thickness: ', num2str(t_film2_calc), ' nm']);
disp(['Substrate Thickness: ', num2str(t_substrate_calc), ' nm']);
