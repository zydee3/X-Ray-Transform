        function [fDataO, betaO,alphaO] = geoA_interp(obj, Beta,Alpha,FData, method,sign)
            [BetaScatt,AlphaScatt] = obj.scatteringRelation(Beta,Alpha);
            [fDataO, betaO,alphaO] = obj.geoA_precomp_interp(Beta,Alpha,FData, BetaScatt,AlphaScatt, method,sign);
        end
        
        function [fDataO, betaO,alphaO] = geoA_precomp_interp(~, Beta,Alpha,FData, BetaScatt,AlphaScatt, method, sign)
            
            % BetaScatt,AlphaScatt must have the same size as Beta,Alpha,FData
            % Values are expected to be in NDGRID format for gridded
            % interpolation (importantly, same beta occupies same row)
            
            % outputs are shaped in preperation for hilbert transform
            
            % methods: 
            %   gridded_nearest     - handles Beta,Alpha,FData with griddedinterpolant nearest method
            %   gridded_linear      - handles Beta,Alpha,FData with griddedinterpolant linear method
            %   scattered_nearest   - handles Beta,Alpha,FData with scatteredinterpolant nearest method
            %   scattered_linear    - handles Beta,Alpha,FData with scatteredinterpolant linear method
            
            switch method
                case {'gridded_nearest','gridded_linear'}
                    f = griddedInterpolant(Beta,Alpha,FData); % defaults to linear
                    if (strcmp(method,'gridded_nearest'))
                        f.Method = 'nearest'; 
                    end    
                    
                case {'scattered_nearest','scattered_linear'}    
                    f = scatteredInterpolant(Beta,Alpha,FData); % defaults to linear
                    if (strcmp(method,'scattered_nearest'))
                        f.Method = 'nearest'; 
                    end  
                    
                otherwise    
                    error('wrong method')
            end
            
            
            po2 = pi/2;
            fDataO = [FData, sign*f(mod(BetaScatt,2*pi), mod(AlphaScatt+po2,pi)-po2)];
            betaO = [Beta, Beta];
            alphaO = [Alpha, AlphaScatt]; % alpha on [0,2pi]
        end      
        
        function [fDataO, betaO,alphaO] = geoA_simple(~, Beta,Alpha,FData, sign)
            
            % For the case that FData(Beta,Alpha) = -FData(BetaScatt,AlphaScatt)
            % Values are expected so that same betas share same row
            
            fDataO = [FData, sign*FData];
            betaO = [Beta, Beta];
            alphaO = [Alpha, Alpha+pi]; % alpha on a weird interval
        end      
        
        function [fO] = geoA_simplefunc(~, F, sign)
            
            % For the case that -FData(Beta,Alpha) = FData(BetaScatt,AlphaScatt)
            % F is handled as an anonymous function
                        
            fO = @Af;
            
            function out = Af(Beta,Alpha)
                
                Beta = mod(Beta,2*pi);
                Alpha = mod(Alpha+pi,2*pi)-pi;
                
                bool = abs(Alpha)>pi/2;
                indS = bool;
                ind = ~bool;
                
                out = zeros(size(Alpha));
                
                out(indS) = sign*F(Beta,Alpha);
                out(ind) = F(Beta,Alpha);
               
            end
        end  
        
 
        function fO = geoHilbert(obj, FData, m, bw) 
            % it is assumed that FData is arranged so that same betas share same rows and that alphas are equispaced
            % additionally, alpha is expected to be cyclically ordered wrt 
            
            sze = size(FData);
            nal = sze(2)/2;
            nbt = sze(1);
            
            H = 2*nal*[zeros(nbt,1), -1j*ones(nbt,nal), 1j*ones(nbt,nal-1)];
            
            % filtering row by row
            fO = ifft(H .* fft(FData,2), 2) /(2*nal);
        end   
        
        function fO = geoHilbert_func(obj, F)
            % F is handled as an anonymous function
            % descretizes the function along alpha in order to perform FFT
            % output is a function handle with respect to beta that returns
            % a list of alphas
            
            
            fO = @Hilbertf;
            
            function out = Hilbertf(Beta,nal)
                nal = floor(nal/2);
                sze = size(Beta);
                Beta = Beta(:);
                [bt,al] = meshgrid(unique(Beta), linspace(-pi/2,3*pi/2,nal));
                
                H = 2*nal*[zeros(nbt,1), -1j*ones(nbt,nal), 1j*ones(nbt,nal-1)];

                % descretize the function along alpha, perform fast
                % hilbert transform
                out = ifft(H .* fft( F(bt,al), nal )) /(2*nal);
                
            end
        end   
        
        function fO = geoHilbert_funcinterp(obj, F, nal)
            % F is handled as an anonymous function
            % descretizes the function along alpha in order to perform FFT
            % output is a function handle
            
            al = linspace(-pi/2,3*pi/2,nal);
            
            fO = @Hilbertf;
            
            function out = Hilbertf(Beta,Alpha)
                sze = size(Beta);
                Beta = Beta(:);
                [ubt] = unique(Beta);
                
                H = 2*nal*[zeros(nbt,1), -1j*ones(nbt,nal), 1j*ones(nbt,nal-1)];

                out = NaN(sze);
                for bt = ubt
                    alind = (bt == Beta); 
                    % descretize the function along alpha, perform fast
                    % hilbert transform and interpolate results
                    
                    % could speed up with interp1q or by removing for loop
                    % with vectorization
                    out(alind) = interp1(al, ifft(H .* fft( F(ones(1,nal)*bt,al), nal )) /(2*nal),...
                                         Alpha(alind), 'cubic');
                end
                
            end
        end   
        

        function [fDataO, betaO,alphaO] = geoAstar_interp(obj, Beta,Alpha,FData, method,sign)
            [BetaScatt,AlphaScatt] = obj.scatteringRelation(Beta,Alpha);
            [fDataO, betaO,alphaO] = obj.geoAstar_precomp_interp(Beta,Alpha,FData, BetaScatt,AlphaScatt, method,sign);
        end
        
        function [fDataO, betaO,alphaO] = geoAstar_precomp_interp(~, Beta,Alpha,FData, BetaScatt,AlphaScatt, method, sign)
            
            % BetaScatt,AlphaScatt must have the same size as Beta,Alpha,FData
            % Values are expected to be outputs of geoHilbert
                        
            % methods: 
            %   gridded_nearest     - handles Beta,Alpha,FData with griddedinterpolant nearest method
            %   gridded_linear      - handles Beta,Alpha,FData with griddedinterpolant linear method
            %   scattered_nearest   - handles Beta,Alpha,FData with scatteredinterpolant nearest method
            %   scattered_linear    - handles Beta,Alpha,FData with scatteredinterpolant linear method
            
            nal = length(F)/2;
            
            switch method
                case {'gridded_nearest','gridded_linear'}
                    f = griddedInterpolant(Beta,Alpha,FData(:,nal+1:end)); % defaults to linear
                    if (strcmp(method,'gridded_nearest'))
                        f.Method = 'nearest'; 
                    end    
                    
                case {'scattered_nearest','scattered_linear'}    
                    f = scatteredInterpolant(Beta,Alpha,FData(:,nal+1:end)); % defaults to linear
                    if (strcmp(method,'scattered_nearest'))
                        f.Method = 'nearest'; 
                    end  
                    
                otherwise    
                    error('wrong method')
            end
            
            fDataO = [FData(:,nal+1:end), + sign*f(BetaScatt(:,nal+1:end),AlphaScatt(:,nal+1:end))];
            betaO = Beta(:,1:nal);
            alphaO = Alpha(:,1:nal);
        end      
        
        function [fDataO, betaO,alphaO] = geoAstar_simple(~, Beta,Alpha,FData, sign)
            
            % For the case that FData(Beta,Alpha) = -FData(BetaScatt,AlphaScatt)
            % Values are expected so that same betas share same row
            
            sze = size(FData);
            nal = sze(2)/2;
            nbt = sze(1);
            
            fDataO = [FData(:,1:nal, sign*FData(:,nal+1:end];
            betaO = Beta(:,1:nal);
            alphaO = Alpha(:,1:nal); % alpha on a weird interval
        end 
        
        function [fO] = geoAstar_simplefunc(~, F, sign)
            
            % For the case that FData(Beta,Alpha) = -FData(BetaScatt,AlphaScatt)
            % F is handled as an anonymous function
                        
            fO = @Astarf;
            
            function out = Astarf(Beta,Alpha)
                Beta = mod(Beta,2*pi);
                
                out = F(Beta,mod(Alpha+pi,2*pi)-pi) + sign*F(Beta,mod(Alpha,2*pi)-pi);
                
            end
            
        end          
        