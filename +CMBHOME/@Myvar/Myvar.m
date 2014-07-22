classdef Myvar
% Standard format for time-varying measurements not defined elsewhere in
% the CMBHOME package. Allows the storage and retrieval of time-varying
% variables measured at arbitary timescales. 
%
% Variables will be interpolated to the same time scale as root.fs_video
% and will undergo epoch filtering just the same as other variables in the
% Session object. The raw timestamps and variable measurements are always
% stored and can be accessed by root.Myvar2.p_get('variableName') and
% root.Myvar2.ts_get('variableName'). See documentation for more
% information ('doc CMBHOME.Myvar')
%
% See Documentation for more info: doc CMBHOME.Myvar

% wchapman 22/07/2013

    properties (SetAccess=private, Hidden) 
        data
    end
    
    properties (Hidden) 
        
    end
    
    properties (Dependent=false, SetAccess=private,Hidden=false)
        map % Cell array of variable names stored in object
    end
    
    properties (Dependent=true, SetAccess=private) %When epoch is set
      
    end
    
    methods
        %% MyvarHolder()
        function self = Myvar(name,p_var,p_ts,tsVar,isCirc)
             import CMBHOME.Utils.* %need CMBHOME.Utils.MyvarHolder
             
             if nargin == 0     
                 % This is default contructor, just returns blank
                 self.map = [];
             elseif (isa(name,'CMBHOME.Utils.MyvarHolder'))
                 self = self.add(name);
             else       % Calls CMBHOME.Utils.MyvarHolder first and returns non-blank
                 if ~exist('tsVar','var'),tsVar = [];end
                 if ~exist('isCirc','var'),isCirc = 0;end
                 self = self.add(name,p_var,p_ts,tsVar,isCirc);
             end
        end
        
        %% Myvar.get()
        function var = get(self,fv)
            % x = root.myvar2.get('x'); 
            % Finds & returns the variable 'x' from myvar2. Filtered by
            % epochs (like root.x)
           
            ind = find(strcmp(fv,self.map));
            if isempty(ind) || ind==0
                error('That variable is not loaded. Check root.myvar2.map');
            end
            
            var = self.data(ind).var;
        end
        
        %% Myvar.b_get()
        function b_var = b_get(self,fv,inds)
            % interpolated_x = root.myvar2.b_get('x'); 
            % Finds & returns the variable 'x' from myvar2. "b_" refers to
            % the up/down sampled data that has not yet been filtered to
            % root.epoch.

            ind = strcmp(fv,self.map);
            if isempty(ind)
                error('That variable is not loaded. Check MyVar.map');
            end

            b_var = self.data(ind).b_var;    
            
            if exist('inds','var') && ~isempty('inds','var')
                b_var = b_var(inds);
            end
        end
        
        %% Myvar.p_get()
        function p_var = p_get(self,fv)
            % raw_x = root.myvar2.p_get('x'); 
            % Finds & returns the variable 'x' from myvar2. "p_" refers to
            % the RAW variable, before interpolation or epochs

            ind = strcmp(fv,self.map);
            if isempty(ind)
                error('That variable is not loaded. Check MyVar.map');
            end

            p_var = self.data(ind).p_var;    
        end
        
        %% Myvar.ts_get()
        function raw_ts = ts_get(self,fv)
            % x = root.myvar2.get('x'); 
            % Finds & returns the RAW time stamps from
            % root.myvar2.data(ind). For the time stamps of the spikes or
            % up/downsamples just look at root.ts/root.b_ts.

            ind = strcmp(fv,self.map);
            if isempty(ind)
                error('That variable is not loaded. Check MyVar.map');
            end

            raw_ts = self.data(ind).p_ts;    
        end
        
        %% isCirc_get()
        function isCirc = isCirc_get(self,fv)
            ind = strcmp(fv,self.map);
            isCirc = self.data(ind).isCirc;
        end
        %% Myvar.remove()
        function [self success] = remove(self,fv)
            % root.myvar2 = root.myvar2.remove('x');
            %
            % Removes the input field and all related fields (eg
            % b_fieldname) from myvar2.
            
            ind = strcmp(fv,self.map);
            
            if isempty(ind)
                warning('That variable is not loaded. Check MyVar.map');
                success = 0;
            else
                self.data(ind) = [];
                self.map(ind) = [];  
                success = 1;
            end
        end
        
        %% Myvar.add()
        function self = add(self,myvar,p_var,p_ts,tsVar,isCirc)
            % root.myvar2 = root.myvar2.add(myvar_holder)
            
            %correctClass = isa(p_var,'single') | isa(p_var,'double') | isa(p_var,'float') | isa(p_var,'double');
            
            %if ~correctClass
            %    error('CMBHOME.Myvar.add(): p_var must be of class double,single or float')
            %end
            
            if nargin >2 % If we have more than self & myvar then construct first
               name = myvar;
               if ~exist('tsVar','var'),tsVar = [];end
               if ~exist('isCirc','var'),isCirc = 0;end
               myvar = CMBHOME.Utils.MyvarHolder(name,p_var,p_ts,tsVar,isCirc);
            end
            
            
            if isempty(self.data)
                self.data = myvar;
                self.map = {};
                self.map{1} = myvar.name;
            else
                self.data(end+1) = myvar;
                self.map(end+1) = {myvar.name};
                self.map = self.map(:);
            end
        end
        
        %% Myvar.setEpoch()
        function self = setEpoch(self,inds,root_b_ts)
            % root.myvar2 = root.myvar2.setEpoch(inds);
            %
            % Sets the epoch using indices from root.ind. This function is
            % automatically called whenever the epoch in root is changed.
            
            for i = 1:length(self.map)
                self.data(i) = self.data(i).setEpoch(inds,root_b_ts);
            end
            
        end
        
        %% Myvar.mean()
        function mv = mean(self,fv,inds,epNumber)
            ind = strcmp(fv,self.map);
            
            if ~exist('inds','var')
                inds = [];
            end
            
            if ~exist('epNumber','var')
                epNumber = [];
            end
            
            mv = self.data(ind).mean(inds,epNumber);
        end
    end
    
end