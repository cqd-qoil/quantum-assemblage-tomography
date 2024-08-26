function [MTestResult, Assem, T, epsx, Lag, Lag0] = tomography_m2(bobData)
    % m2 - Perform assemblage tomography on quantum measurement data (bobData) using
    % conical optimization techniques to estimate: quantum assemblage and reduced
    % quantum state of the trusted pary, detection efficiencies of the
    % untrusted party, and their corresponding likelihoods.
    % Assumes the intermediate loss model described in the manuscript "Quantum
    % Assemblage Tomography," where detection efficiencies are only influenced
    % by measurement choices.
    % It utilizes the YALMIP toolbox and the MOSEK solver to perform a series of
    % calculations, including an implementation of a multinomial test on
    % the data. The output includes the results of the multinomial test, the
    % optimized quantum assemblage, and other relevant parameters.

    % MTestResult: Outcome of the multinomial test, indicating whether
    %               the observed frequencies in bobData significantly
    %               deviate from the expected frequencies of a multinomial
    %               distribution.
    %               See https://en.wikipedia.org/wiki/Multinomial_test
    %
    % Assem: Quantum assemblage that represents the collection of quantum
    %        states corresponding to the measurement outcomes in bobData.
    %
    % T: The reduced state of Bob (rho_B).
    %
    % epsx: Array of Alice's estimated detection efficiencies.
    %
    % Lag: Log likelihood of the optimization.
    %
    % Lag0: Log likelihood of the baseline model, serving as a reference
    %       point for comparison.

    % Copyright (c) 2024 Yuanlong Wang (and Luis Villegas-Aguilar)

    % Permission is hereby granted, free of charge, to any person obtaining a copy
    % of this software and associated documentation files (the "Software"), to deal
    % in the Software without restriction, including without limitation the rights
    % to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    % copies of the Software, and to permit persons to whom the Software is
    % furnished to do so, subject to the following conditions:
    %
    % The above copyright notice and this permission notice shall be included in
    % all copies or substantial portions of the Software.
    %
    % THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    % IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    % FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    % AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    % LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    % FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
    % IN THE SOFTWARE.

    % ---------------------------------------------------------------------------%

    % Ensure MOSEK is installed and set up.

    % Set YALMIP solver settings to use MOSEK
    yalmip('clear');

    % Set up solver options
    ops = sdpsettings('solver', 'mosek', 'verbose', 0);

    % Manual configuration
    LastTwoSum = 1; % Indicates if the last two columns of the data are total measurements (1) or null results (0)

    % Load data from bobData into a structured array
    data = zeros(9, 6); % 9: Alice's POVMs * Bob's POVMs, 6: outcomes
                        % hh hv vh vv nullH nullV
    data(:,:) = bobData;

    if LastTwoSum ~= 0 % 1 means last two columns are total measurement shot numbers
        data(:, 5:6) = data(:, 5:6) - data(:, 1:2) - data(:, 3:4);
    end

    % Initialize single-qubit Pauli matrices
    sigma = zeros(2, 2, 4);
    sigma(:,:,1) = [1, 0; 0, -1];
    sigma(:,:,2) = [0, 1; 1, 0];
    sigma(:,:,3) = [0, -1i; 1i, 0];
    sigma(:,:,4) = eye(2);

    % Calculate total number of copies for each POVM group
    NEG = sum(data, 2);

    % Define the number of outcomes and POVMs
    aMax = 3; % Number of different outcomes for Alice's POVMs
    xMax = size(data, 1) / 3; % Number of different POVMs on Alice's side
    bMax = 2; % Number of different outcomes for Bob's POVMs
    yMax = 3; % Number of different POVMs on Bob's side

    % Initialize POVM matrices for Bob
    E = zeros(2, 2, bMax, yMax);
    for j = 1:yMax
        E(:,:,1,j) = (eye(2) + sigma(:,:,j)) / 2;
        E(:,:,2,j) = (eye(2) - sigma(:,:,j)) / 2;
    end

    % Check parameter consistency
    if yMax * xMax ~= size(data, 1) || bMax * aMax ~= size(data, 2)
        warndlg('Specified parameters are not consistent!', 'Warning');
    end

    % Estimate the general detection efficiency
    epsilon = sum(data(:, 1:4), 'all') / sum(data, 'all');
    % epsilon = sum(sum(sum(data(:,1:4)))) / sum(sum(sum(data(:,:))));

    % Initialize assemblage and detection efficiency
    Assem = zeros(2, 2, aMax, xMax);
    epsx = ones(size(data, 1) / 3, 1) * sum(sum(sum(data(:,1:4)))) / sum(sum(sum(data)));

    % Initial guess for assemblage
    T = zeros(2, 2);
    rho0 = eye(4) / 4;
    for a = 1:aMax
        for x = 1:xMax
            T(:,:) = rho0(1:2, 1:2) + rho0(3:4, 3:4);
            Assem(:,:,a,x) = T(:,:) * ((a < aMax) * epsilon / (aMax - 1) + (a == aMax) * (1 - epsilon));
        end
    end

    % Calculate initial likelihood
    LagStart = 0;
    for a = 1:aMax
        for b = 1:bMax
            for x = 1:xMax
                for y = 1:yMax
                    LagStart = LagStart + data(yMax*(x-1)+y, bMax*(a-1)+b) * log(trace(Assem(:,:,a,x) * E(:,:,b,y)));
                end
            end
        end
    end

    % Optimization setup
    Sig = cell(aMax, xMax);
    for a = 1:aMax
        for x = 1:xMax
            Sig{a,x} = sdpvar(2, 2, 'full', 'complex');
            Sig{a,x} = (Sig{a,x} + Sig{a,x}') / 2;
        end
    end

    T = sdpvar(2, 2, 'full', 'complex');
    T = (T + T') / 2;
    Fcons = [];

    for x = 1:xMax
        Fcons = [Fcons, Sig{1,x} + Sig{2,x} == T * epsx(x), ...
                 Sig{1,x} + Sig{2,x} + Sig{3,x} == T, ...
                 Sig{1,x} >= 0, Sig{2,x} >= 0, Sig{3,x} >= 0];
    end

    Fcons = [Fcons, T >= 0, trace(T) == 1];

    % Objective function
    obj = 0;
    for x = 1:xMax
        for b = 1:bMax
            for y = 1:yMax
                obj = obj + data(yMax*(x-1)+y, bMax*(1-1)+b) * log(real(trace(Sig{1,x} * E(:,:,b,y)))) + ...
                             data(yMax*(x-1)+y, bMax*(2-1)+b) * log(real(trace(Sig{2,x} * E(:,:,b,y)))) + ...
                             data(yMax*(x-1)+y, bMax*(3-1)+b) * log(real(trace(Sig{3,x} * E(:,:,b,y))));
            end
        end
    end

    optimize(Fcons, -obj / 10^5, ops);

    T = value(T);

    for x = 1:xMax
        for a = 1:aMax
            Assem(:,:,a,x) = value(Sig{a,x});
        end
    end

    % Calculate the corresponding likelihood
    Lag = 0;
    Lag0 = 0;
    for b = 1:bMax
        for x = 1:xMax
            for y = 1:yMax
                Lag = Lag + data(yMax*(x-1)+y, bMax*(1-1)+b) * log(real(trace(Assem(:,:,1,x) * E(:,:,b,y)))) + ...
                             data(yMax*(x-1)+y, bMax*(2-1)+b) * log(real(trace(Assem(:,:,2,x) * E(:,:,b,y)))) + ...
                             data(yMax*(x-1)+y, bMax*(3-1)+b) * log(real(trace(Assem(:,:,3,x) * E(:,:,b,y))));

                ndata = sum(data(yMax*(x-1)+y, :));
                Lag0 = Lag0 + data(yMax*(x-1)+y, bMax*(1-1)+b) * log(data(yMax*(x-1)+y, bMax*(1-1)+b) / ndata) + ...
                               data(yMax*(x-1)+y, bMax*(2-1)+b) * log(data(yMax*(x-1)+y, bMax*(2-1)+b) / ndata) + ...
                               data(yMax*(x-1)+y, bMax*(3-1)+b) * log(data(yMax*(x-1)+y, bMax*(3-1)+b) / ndata);
            end
        end
    end

    % Check non-signaling condition
    teste1 = zeros(xMax);
    errMax = 0;
    nonsignal = zeros(2, 2, xMax);

    for x = 1:xMax
        for a = 1:aMax
            nonsignal(:,:,x) = nonsignal(:,:,x) + Assem(:,:,a,x);
            [vs, ds] = eig(Assem(:,:,a,x));
            if sum(diag(ds) < 0) > 0
                errMax = 1; % Check if assemblage is positive semidefinite
            end
        end
        errMax = max(errMax, norm(nonsignal(:,:,1) - nonsignal(:,:,x), 'fro'));
        teste1(x) = norm(nonsignal(:,:,1) * epsx(x) - Assem(:,:,1,x) - Assem(:,:,2,x), 'fro');
    end
    errMax = max(errMax, abs(1 - trace(nonsignal(:,:,1))));

    % Calculate ideal probabilities that would be generated by our estimation
    % in the absence of additional noise
    idealp = zeros(size(data));
    for y = 1:yMax
        for b = 1:bMax
            for x = 1:xMax
                for a = 1:aMax
                    idealp(yMax*(x-1)+y, bMax*(a-1)+b) = real(trace(Assem(:,:,a,x) * E(:,:,b,y)));
                end
            end
        end
    end

    % Calculate experimental frequency and multinomial test
    expf = zeros(size(data));
    MultinomialTest = zeros(size(data, 1), 1);
    for povmn = 1:size(data, 1)
        pmle0 = idealp(povmn, :);
        ndata1 = data(povmn, :);
        ndata = sum(ndata1);
        tdata = ndata1 / ndata;
        expf(povmn, :) = tdata;
        MultinomialTest(povmn) = -2 * sum(ndata1 .* log(pmle0 ./ (tdata + (tdata == 0)))) / ...
            (1 + (sum(1 ./ pmle0) - 1) / (6 * ndata * (size(data, 2) - 1)));
    end

    MTestResult = sum(MultinomialTest); % Chi-squared distribution
    temv1 = (data > -1);
    temv2 = (idealp - expf) .* temv1;
    temv3 = expf .* temv1;
    RelativeFreqErr = (norm(temv2(:), 'fro')) / (norm(temv3(:), 'fro'));

    % Uncomment to display results
    % disp('[errMax, RelativeFreqErr, MTest]')
    % disp([errMax, RelativeFreqErr, MTestResult])
end
