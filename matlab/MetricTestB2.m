function MetricTestB3()
    met = Metric.build('polynomial', 'coeffs', rand(1,6));
    
    metOLD.cfs = met.coeffs;
    
    x = rand(200);
    y = rand(size(x));
    
    pause(1);
    tic
        met.metricValsCurv(x, y);
        %metricValuesCurvature3(x, y, metOLD);
        
    toc
    pause(1);
    tic
        metricValuesCurvature2(x, y, metOLD);       
        %met.metricValsCurv(x, y);
    toc
        
end