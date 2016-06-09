classdef MyvarHolder
    % Class for arbitrary time series variables. 
    
    properties (SetAccess=private, Hidden) 
        p_var       % Backup of the original variable, from the time of making the object
        p_ts        % Backup of the original timestamps, from the time o fmaking the object
        checkRoot   % Boolean, checks whether or not we have aligned MyVarHolder.b_ts to root.b_ts
    end
    
    properties (Hidden) %All of ts & var, interpolated to another variable.
        b_ts    % Unfiltered vector of timestamps after interpolating to root.ts
        b_var   % Unfiltered vector of variable after interpolating
        inds % Filters down to an epoch (cell array)
    end
    
    properties (Dependent=false, SetAccess=private,Hidden=false)
        name
        isCirc
    end
    
    properties (Dependent=true, SetAccess=private) %When epoch is set
        var
        ts
    end
    
    methods
       
        %% Constructor
        function self = MyvarHolder(name,p_var,p_ts,tsVar,isCirc)

            self.p_var = p_var(:);
            self.p_ts = p_ts(:);
            
            if length(self.p_var) ~=length(self.p_ts)
               error('var and ts must be the same length! \n') 
            end
            
            if ~exist('tsVar','var')
                path=0;
            else 
                if isempty(tsVar)
                    path = 0;
                else
                    path = 1;
                end
            end
            
            %%%%
            if exist('isCirc','var')
                if ~isempty(isCirc)
                    if ~((isCirc==1) || (isCirc==0))
                        error('isCirc must be 0 or 1')
                    end
                else
                    isCirc=0;
                end
            else
                isCirc = 0;
            end
            
            self.isCirc = isCirc;
            
            %%%%
            
            if (self.isCirc) && (((max(p_var) > 2*pi) && (min(p_var) < 0)) || ((max(p_var) > pi) && (min(p_var) < -pi)))
                warning('Must use variables in radians ranging from -pi to pi. You appear to be using degrees. Converting Automatically');
                p_var = deg2rad(p_var);
            end
            
            if (self.isCirc)
                p_var = wrapToPi(p_var);
            end
            
            %%%%
            
            if (path==0)
                fprintf('no tsVar input. Assuming same as signal ts \n');
                tsVar = p_ts(:);
                self.checkRoot = 0; 
                self.b_var = p_var;
                self.b_ts = p_ts;
            elseif (path==1)
                self.checkRoot = 1;
                extra = (self.p_ts > max(tsVar)) | (self.p_ts < min(tsVar)); % these are outside of root.b_ts
                self.b_ts = self.p_ts(~extra);
                self.b_var = self.p_var(~extra);
                self.b_var = interp1(self.b_ts,self.b_var,tsVar);
                self.b_ts  = tsVar;
            end
            
            self.inds = (1:length(tsVar))';
            self.name = name;
            
        end
        
        %% Sets
        function self = setEpoch(self,epochInds,root_b_ts)
            if self.checkRoot == 0
                % if we didn't get the root.b_ts during initialization of
                % the variable then need to check now that Myvar.b_ts is
                % same as root.b_ts.
                self.checkRoot = 1;
                tsVar = root_b_ts;
                extra = (self.p_ts > max(tsVar)) | (self.p_ts < min(tsVar)); % these are outside of root.b_ts
                self.b_ts = self.p_ts(~extra);
                self.b_var = self.p_var(~extra);
                self.b_var = interp1(self.b_ts,self.b_var,tsVar);
                self.b_ts  = tsVar;
            end
            
            self.inds = epochInds;
            
        end
        
        %% Gets
        function var = get.var(self)
            if iscell(self.inds)
                var = cellfun(@(x) self.b_var(x), self.inds,'UniformOutput',0); 
            else  
                var = self.b_var(self.inds);
            end
        end
        
        function ts = get.ts(self)
            if iscell(self.inds)
                ts = cellfun(@(x) self.b_ts(x), self.inds,'UniformOutput',0); 
            else
                ts = self.b_ts(self.inds); 
            end
        end
        
        function p_var = get.p_var(self)
           p_var = self.p_var; 
        end
        
        function p_ts = get.p_ts(self)
           p_ts = self.p_ts; 
        end
        
        function namev = get.name(self)
            namev=self.name;
        end
        
        function isCirc = get.isCirc(self)
            isCirc = self.isCirc;
        end
        
        %% Mean functions
        function mv = mean(self,inds,epNumber)
            % For a single myVar, finds the the mean at defined points. 
            % Outputs: 
            %       mv: The mean of the variable (vec)
            % Inputs:
            %       inds: Indinces for each epoch (vec/cell)
            %       epNumber: If using multiple epochs then must pass in a
            %                 vector of which epoch each index cell is for.
            
            if ~exist('inds','var')
               inds = [];
            end 
            
            if ~exist('epNumber','var')
                epNumber = [];
            end
            
            if ~iscell(self.var)
               if isempty(inds)
                inds = self.inds;
               end
               
               if self.isCirc
                   var = self.var(inds);
                   var(isnan(var)) = [];
                   mv = CMBHOME.Utils.circ_mean(var);
               else
                   mv = nanmean(self.var(inds));
               end
               
            else
                
                if isempty(inds) && isempty(epNumber)
                    for i = 1:length(self.var)
                        if self.isCirc
                            var = self.var{i};
                            bads = isnan(var);
                            var(bads) = [];
                            mv(i) = CMBHOME.Utils.circ_mean(var);
                        else
                            mv(i) = nanmean(self.var{i});
                        end
                    end
                elseif isempty(inds) && ~isempty(epNumber)
                    for i = 1:length(epNumber)
                        if self.isCirc
                            var = self.var{epNumber(i)};
                            bads = isnan(var);
                            var(bads) = [];
                            mv(i) = CMBHOME.Utils.circ_mean(var);
                        else
                            mv(i) = nanmean(self.var{epNumber(i)});
                        end
                    end
                elseif ~isempty(inds) && isempty(epNumber)
                    error('If using myvar2.mean when multiple epochs, must pass in epNumber as well. (doc CMBHOME.Utils.MyvarHolder.mean)')
                else % ~isempty(inds) & exist('epNumber','var')
                    for i = 1:length(epNumber)
                       if self.isCirc
                            var = self.var{epNumber(i)}(inds{i});
                            bads = isnan(var);
                            var(bads) = [];
                            mv(i) = CMBHOME.Utils.circ_mean(var);
                       else
                           mv(i) = nanmean(self.var{epNumber(i)}(inds{i}));
                       end
                    end
                end    
                
                mv = mv(:);
            end
            
            
        end
        
    end
    
end