function [ results ] = pass_index( self, varargin )
if isempty(self.cel)
    self.cel = varargin{1};
    varargin = varargin(2:end);    
end
if isempty(self.active_lfp)
    self.active_lfp = self.cel(1);
        varargin = varargin(2:end);    
end
results = CMBHOME.PASS_INDEX.pass_index(self.ts,[self.sx self.sy],self.spk.ts,self.lfp.ts,self.lfp.signal,varargin{:});

end

