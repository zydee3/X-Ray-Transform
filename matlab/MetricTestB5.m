function MetricTestB5()
    met = Metric.build('euclid');
    
    metOLD.cfs = 12312;
    
    x = rand(200);
    y = rand(size(x));
    
    pause(1);
    tic
        met.metricValsCurv(x, y);
        %metricValuesCurvature3(x, y, metOLD);
        
    toc
    pause(1);
    tic
        metricValuesCurvature5(x, y, metOLD);       
        %met.metricValsCurv(x, y);
    toc
        
end