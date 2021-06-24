function MetricTestB3()
    met = Metric.build('sphere', 'radius', rand(1));
    
    metOLD.R = met.radius;
    
    x = rand(1000);
    y = rand(size(x));
    
    pause(1);
    tic
        %metricValuesCurvature3(x, y, metOLD);
        met.metricValsCurv(x, y);
    toc
    pause(1);
    tic
        metricValuesCurvature3(x, y, metOLD);
        %met.metricValsCurv(x, y);
    toc
        
end