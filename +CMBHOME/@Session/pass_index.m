function [ results ] = pass_index( self, varargin )
% How far through the field? 

if isempty(self.cel)
    self.cel = varargin{1};
    varargin = varargin(2:end);    
    
end
if isempty(self.active_lfp)
    self.active_lfp = self.cel(1);
        varargin = varargin(2:end);    
end
results = CMBHOME.PASS_INDEX.pass_index('pos_ts',self.ts,'pos',[self.sx self.sy],'spk_ts',self.cel_ts{1},'lfp_ts',self.lfp.ts,'lfp_sig',self.lfp.signal,varargin{:});

end

