%% To construct a sphere metric:
%(Domain and metric have been coded in basically the same way)

    R = 0.7;

    ms0 = sphereMetric; % initializes a default sphere metric
    ms0.radius = R; % we can modify parameters by name after creation

    ms1 = sphereMetric(radius = R); %equivically, this initializes the same metric as above

%you can also construct the same metric directly from a function handle expression of log(g), 
%but this is less efficient and can be prone to error
    ms3 = Metric( @(x,y) log(4*R^4) - 2*log(x.*x + y.*y + R*R) );



%% To construct a cosine domain:
% The order and number of parameters referenced doesn't matter (this is also true for metrics)
    dcos0 = cosineDomain; % a domain constructed with the default parameters
    dcos1 = cosineDomain(amplitude = 0.5, radius = 7); % doesnt reference 'cycles' property
    dcos2 = cosineDomain(cycles = 4, radius = 7 , amplitude = 0.5); % does reference 'cycles' property

% After initialization, domains can be transformed using the transform
% method or by changing the variables individually
    dcos0 = dcos0.transform(1,0,pi*0.25);

    %equivically,
    dcos1.originX = 1;
    dcos1.theta = pi*0.25;



%% To construct a riemann surface:

    rs0 = RiemannSurface; % a surface constructed with the default parameters
    rs1 = RiemannSurface(dcos2, ms1); % a surface constructed with the above domain and metric (domain must always preceed the metric)

    rs2 = RiemannSurface(dcos2, stepType = 'EE'); % we can also modify the stepping method properties during or after initialization 

    rs1.stepType = 'EE';


