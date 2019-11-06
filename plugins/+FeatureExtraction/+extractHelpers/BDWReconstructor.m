% This file is part of DAVE, a MATLAB toolbox for data evaluation.
% Copyright (C) 2018-2019 Saarland University, Author: Manuel Bastuck
% Website/Contact: www.lmt.uni-saarland.de, info@lmt.uni-saarland.de
% 
% The author thanks Tobias Baur, Tizian Schneider, and Jannis Morsch
% for their contributions.
% 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU Affero General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU Affero General Public License for more details.
% 
% You should have received a copy of the GNU Affero General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>. 

classdef BDWReconstructor < FeatureExtraction.extractHelpers.ReconstructorInterface
    %BDWRECONSTRUCTOR Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        ind = [];
        sensor = [];
    end
    
    methods
        function this = BDWReconstructor(ind, sensor, target)
            if nargin > 0
                this.ind = ~ind;
                this.sensor = sensor;
                this.sensor.addNextElement(this);
            end
            if nargin > 2
                this.target = target;
            end
        end
        
        function res = getResiduals(this, cycNum)
            %get original cycle Data
            metaInfo = containers.Map('KeyType','char','ValueType','any');
            metaInfo('ID') = '';
            this.sensor.step(cycNum, metaInfo);
            data = this.originalCycle;
            
            %perform wavelet transformation
            [data, l] = wavedec(data, wmaxlev(length(data),'db2'), 'db2');
            
            %zero unselected coefficients
            data(this.ind) = 0;
            
            %reconstruct signal
            recData = waverec(data, l, 'db2');
            
            res = recData - this.originalCycle;
        end
    end
    
end

