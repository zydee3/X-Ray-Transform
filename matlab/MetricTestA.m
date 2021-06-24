function MetricTestA(met) % Call from command line
   
    [lg,dxlg,dylg,curv] = met.getHandles();

    metNumeric = Metric(lg);
    
    met.plotALL;
    metNumeric.plotALL;
end