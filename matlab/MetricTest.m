function MetricTest(met)
    
    [lg,dxlg,dylg,curv] = met.getHandles()

    metNumeric = Metric(lg)
    
    met.plotALL
    metNumeric.plotALL
end