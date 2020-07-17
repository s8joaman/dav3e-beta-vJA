project = dave.project;

%% Ordner in dem die Dateien gespeichert werden sollen
folder = 'C:\Users\Johannes\Desktop\test';
mkdir(folder)

if ~exist('models', 'var')
    models = struct();
end

%% Alle Sensoren die berechnet werden sollen
sensors = {
%   'sgp30_sensor_0',...
%    'sgp30_sensor_1',...
%    'sgp30_sensor_2',...
%    'sgp30_sensor_3',...
 %  'sgp30_sensor_all',...
   % 'bme680',...
   %'zmod4410',...
   % 'as_mlv_p2_vergiftet',...
  %  'as_mlv_p2',...
   % 'as_mlv',...
   % 'ust_proto',...
    %'BME680_569736729_channel_0' ...
   % 'BME680_569710616_channel_1' ...
   % 'ZMOD4410_77241862418580_channel_1' ...
   % 'ZMOD4410_77241862398351_channel_0' ...
    'SGP30_all_channel_0' ...
    'SGP30_all_channel_1' ...
    };

%% Alle Gruppen die berechnet werden sollen
groups = {...
    ' acetone', ...
    ' toluene', ...
    ' carbon_monoxide', ...
    ' formaldehyde',...
    ' hydrogen', ...
    ' ethanol', ...
    ' TVOC_ppb', ...
    ' TVOC_ugm3' ...
    };

%% Alle Rergessions-Methoden die berechnet werden sollen
method = {
    % 'RFESVR', ...
     'RFEplsr', ...
     'RFEleastsquares', ...
    };

evalF={"sensor","gas","method",...
        "numFeatMinOneStd","nCompMinOneStd",...
        "TrainingError", 'TestingError',... 
       "rank","Cluster","featCap",...
       "fullModelData", "POT"};

for s=1:length(sensors)
    
    disp('**************************************************************')
    disp('**************************************************************')
	disp(['compute ', num2str(s), '/',num2str(length(sensors))])
	disp('**************************************************************')
    
    % deactivate all sensors
    project.getSensors().setActive(false)
    
    sensor = sensors{s};
    data = struct();
    
    %% setzt die unterschiedlichen Sensoren aktiv, die ich haben will
    switch sensor
        case 'sgp30_sensor_0'
            system = 'asic1';
            sel_sensors = {'sgp30_sensor_0'};
            setSensors(project, system, sel_sensors)
        case 'sgp30_sensor_1'
            system = 'asic1';
            sel_sensors = {'sgp30_sensor_1'};
            setSensors(project, system, sel_sensors)
        case 'sgp30_sensor_2'
            system = 'asic1';
            sel_sensors = {'sgp30_sensor_2'};
            setSensors(project, system, sel_sensors)
        case 'sgp30_sensor_3'
            system = 'asic1';
            sel_sensors = {'sgp30_sensor_3'};
            setSensors(project, system, sel_sensors)
        case 'sgp30_sensor_all'
            system = 'asic1';
            sel_sensors = {'sgp30_sensor_0', 'sgp30_sensor_1', 'sgp30_sensor_2', 'sgp30_sensor_3'};
            setSensors(project, system, sel_sensors)
        case 'bme680'
            system = 'asic2';
            sel_sensors = {'bme680'};
            setSensors(project, system, sel_sensors)
        case 'zmod4410'
            system = 'asic2';
            sel_sensors = {'zmod4410'};
            setSensors(project, system, sel_sensors)
       case 'as_mlv_p2'
            system = 'hs1';
            sel_sensors = {'as_mlv_p2'};
            setSensors(project, system, sel_sensors)
        case 'as_mlv_p2_vergiftet'
            system = 'hs1';
            sel_sensors = {'as_mlv_p2_vergiftet'};
            setSensors(project, system, sel_sensors)
        case 'ust_proto'
            system = 'hs2';
            sel_sensors = {'ust_proto'};
            setSensors(project, system, sel_sensors)
        case 'as_mlv'
            system = 'hs2';
            sel_sensors = {'as_mlv'};
            setSensors(project, system, sel_sensors)
        case 'BME680_569736729_channel_0'
            system = 'calibration0';
            sel_sensors = {'BME680_569736729_channel_0'};
            setSensors(project, system, sel_sensors)
        case 'BME680_569710616_channel_1'
            system = 'calibration0';
            sel_sensors = {'BME680_569710616_channel_1'};
            setSensors(project, system, sel_sensors)
        case 'ZMOD4410_77241862418580_channel_1'
            system = 'calibration0';
            sel_sensors = {'ZMOD4410_77241862418580_channel_1'};
            setSensors(project, system, sel_sensors)
        case 'ZMOD4410_77241862398351_channel_0'
            system = 'calibration0';
            sel_sensors = {'ZMOD4410_77241862398351_channel_0'};
            setSensors(project, system, sel_sensors)
        case 'SGP30_all_channel_0'
            system = 'calibration0';
            sel_sensors = {'SGP30_8692385_sensor0_channel_0', 'SGP30_8692385_sensor1_channel_0', 'SGP30_8692385_sensor2_channel_0', 'SGP30_8692385_sensor3_channel_0'};
            setSensors(project, system, sel_sensors)
        case 'SGP30_all_channel_1'
            system = 'calibration0';
            sel_sensors = {'SGP30_8692778_sensor0_channel_1', 'SGP30_8692778_sensor1_channel_1', 'SGP30_8692778_sensor2_channel_1', 'SGP30_8692778_sensor3_channel_1'};
            setSensors(project, system, sel_sensors)
    end
    
    % compute Features
    project.computeFeatures();
    project.mergeFeatures();


    for g=1:length(groups)
        disp('      **************************************************************')
        disp('      **************************************************************')
        disp(['      compute ', num2str(g), '/',num2str(length(groups))])
        disp('      **************************************************************')
        
        for m=1:length(method)
            disp('            **************************************************************')
            disp('            **************************************************************')
            disp(['            compute ', num2str(m), '/',num2str(length(method))])
            disp('            **************************************************************')
            % model = computeModel(project, groups{g}, method{m});
            [model, eval] = computeModel(project, groups{g}, method{m}, sensors{s}, evalR);
            nbtk = size(evalF,1)+1;
            evalF{nbtk,1} = sensors{s};
            evalF{nbtk,2} = groups{g};
            evalF{nbtk,3} = method{m};
            evalF{nbtk,4} = eval.nFeatMinOneStd;
            evalF{nbtk,5} = eval.nCompMinOneStd;
            evalF{nbtk,6} = eval.err.training;
            evalF{nbtk,7} = eval.err.testing;
            evalF{nbtk,8} = eval.rank;
            evalF{nbtk,9} = eval.currentCluster;
            evalF{nbtk,10} = eval.featCap;
            evalF{nbtk,11} = eval.fullModelData;
            evalF{nbtk,12} = eval.pred;

            % file = fullfile(folder, [sensor,'-', groups{g},'-', method{m}, '.mat']);
            % save(file, 'model');
        end
    end
    file = fullfile(folder, [sensor,'-', 'eval' '.mat']);
    save(file, 'evalF');    
end
clearvars ans data eval file folder g groups m method model models nbtk project s sel_sensors sensor sensors system

%% funktion die die passenden Sensoren aktiv setzt
function setSensors(project, system, sensors)
    for i=1:length(sensors)
        sensor = sensors{i};
        switch system
            case 'asic1'
                % setSensorActive(project,'20190218_180207_asic1', sensor) % SEQ1
                setSensorActive(project,'mea_asic1', sensor) % 1
%                 setSensorActive(project,'sgp30_sensor_1', sensor) % 2
%                 setSensorActive(project,'sgp30_sensor_2', sensor) % 3
%                 setSensorActive(project,'sgp30_sensor_3', sensor) % 3
 %               setSensorActive(project,'20190228_090643_asic1', sensor) % 4
%                 setSensorActive(project,'20190302_113000_asic1', sensor) % 5
%                 setSensorActive(project,'20190304_090609_asic1', sensor) % 6
%                 setSensorActive(project,'20190306_130509_asic1', sensor) % 7
%                  setSensorActive(project,'20190308_144021_asic1', sensor) % SEQ2
%                 setSensorActive(project,'20190314_164205_asic1', sensor) % 10
%                 setSensorActive(project,'20190316_090407_asic1', sensor) % 11
%                 setSensorActive(project,'20190318_130225_asic1', sensor) % 12
%                 setSensorActive(project,'20190320_131548_asic1', sensor) % 13
%                 setSensorActive(project,'20190322_134640_asic1', sensor) % 14
            case 'asic2'
                % setSensorActive(project,'20190218_180207_asic2', sensor) % SEQ1
                setSensorActive(project,'20190221_180223_asic2', sensor) % 1
                setSensorActive(project,'20190224_162756_asic2', sensor) % 2
                setSensorActive(project,'20190226_091454_asic2', sensor) % 3
                setSensorActive(project,'20190228_090643_asic2', sensor) % 4
                setSensorActive(project,'20190302_113000_asic2', sensor) % 5
                setSensorActive(project,'20190304_090609_asic2', sensor) % 6
                setSensorActive(project,'20190306_130509_asic2', sensor) % 7
                % setSensorActive(project,'20190308_144021_asic1', sensor) % SEQ2
                setSensorActive(project,'20190314_164205_asic2', sensor) % 10
                setSensorActive(project,'20190316_090407_asic2', sensor) % 11
                setSensorActive(project,'20190318_130225_asic2', sensor) % 12
                setSensorActive(project,'20190320_131548_asic2', sensor) % 13
                setSensorActive(project,'20190322_134640_asic2', sensor) % 14 
            case 'hs1'
                % setSensorActive(project,'20190213_174429_hs1', sensor) % SEQ0
                setSensorActive(project,'mea_hs1', sensor) % SEQ1
                % setSensorActive(project,'20190221_180318_hs1', sensor) % 1
                % setSensorActive(project,'20190224_162850_hs1', sensor) % 2
                % setSensorActive(project,'20190226_091518_hs1', sensor) % 3
                % setSensorActive(project,'20190228_090744_hs1', sensor) % 4
                % setSensorActive(project,'20190302_113106_hs1', sensor) % 5
                % setSensorActive(project,'20190304_090628_hs1', sensor) % 6
                % setSensorActive(project,'20190306_130526_hs1', sensor) % 7
                %  setSensorActive(project,'20190308_144033_hs1', sensor) % SEQ2
                % setSensorActive(project,'20190314_164219_hs1', sensor) % 10
                % setSensorActive(project,'20190316_090416_hs1', sensor) % 11
                % setSensorActive(project,'20190318_130321_hs1', sensor) % 12
                % setSensorActive(project,'20190320_131605_hs1', sensor) % 13
                % setSensorActive(project,'20190322_134659_hs1', sensor) % 14 
            case 'hs2'
                % setSensorActive(project,'20190213_174429_hs2', sensor) % SEQ0
                 setSensorActive(project,'20190218_180309_hs2', sensor) % SEQ1
                setSensorActive(project,'20190221_180318_hs2', sensor) % 1
                setSensorActive(project,'20190224_162850_hs2', sensor) % 2
                setSensorActive(project,'20190226_091518_hs2', sensor) % 3
                setSensorActive(project,'20190228_090744_hs2', sensor) % 4
                setSensorActive(project,'20190302_113106_hs2', sensor) % 5
                setSensorActive(project,'20190304_090628_hs2', sensor) % 6
                setSensorActive(project,'20190306_130526_hs2', sensor) % 7
                 setSensorActive(project,'20190308_144033_hs2', sensor) % SEQ2
                setSensorActive(project,'20190314_164219_hs2', sensor) % 10
                setSensorActive(project,'20190316_090416_hs2', sensor) % 11
                setSensorActive(project,'20190318_130321_hs2', sensor) % 12
                setSensorActive(project,'20190320_131605_hs2', sensor) % 13
                setSensorActive(project,'20190322_134659_hs2', sensor) % 14
            case 'calibration0'
                 setSensorActive(project,'calibration0', sensor) % SEQ0
                  setSensorActive(project,'calibration1', sensor) % SEQ0
%                 setSensorActive(project,'field0', sensor) % SEQ0
%                 setSensorActive(project,'field1', sensor) % SEQ0
%                 setSensorActive(project,'field2', sensor) % SEQ0
%                 setSensorActive(project,'field3', sensor) % SEQ0
%              setSensorActive(project,'calibrationAfterField', sensor) % SEQ0
%              setSensorActive(project,'calibrationLowConcentrations', sensor) % SEQ0
%                setSensorActive(project,'calibrationHighConcentrations', sensor) % SEQ0
                setSensorActive(project,'field4', sensor) % SEQ0
                 setSensorActive(project,'field5', sensor) % SEQ0
                 setSensorActive(project,'field6', sensor) % SEQ0
        end
    end
end

%%
% Modell berechnung, hier können  die Parameter gesetzt werden. Wenn man
% was neues hat einfach sich hParam (Zeile 204) anschauen, wie es aufgebaut
% ist und dann kann man es passend ändern. Das Modell kann in Zeile 201
% ausgewählt werden

function [modelStruct, eval] = computeModel(project, group, method, sensor, evalR)
    % get model
    data = project.mergedFeatureData;
    % data.reducedGrouping = data.groupings;
    model = getModel(project, 'model');

    % get Parameter
    hParam = model.processingChain.getChainHyperParameters();

    % change grouping
    hParam.Annotation_annotate.grouping = group;
    grouping = project.getGroupingByCaption(group);
    hParam.Annotation_annotate.groups = cellstr(grouping.getCategories());
    %
    sgmInd = find(contains(string(evalR(:,1)),sensor) & contains(string(evalR(:,2)),group) & contains(string(evalR(:,3)),method));
    featCap = evalR{sgmInd, 14};
    rank = evalR{sgmInd, 15};
    newFeat=featCap(rank(1:evalR{sgmInd,4}));
    hParam.Annotation_annotate.features = cellstr(newFeat);
    %
    hParam.Annotation_annotate.target = 'same as grouping';
    
    % change DimensionalityReduction method method
%     if strcmp(method, 'RFESVR')
%         Hmethod = 1;
%     elseif strcmp(method,'RFEplsr')
%         Hmethod = 5;
%     elseif strcmp(method,'RFEleastsquares')
%         Hmethod = 6;
%     else
%         warning('Something went wrong')
%     end
    
%    hParam.DimensionalityReduction_automatedMethods.methods = Hmethod;
%    sgmInd = find(contains(string(evalR(:,1)),sensor) & contains(string(evalR(:,2)),group) & contains(string(evalR(:,3)),method));
%    hParam.DimensionalityReduction_automatedMethods.numFeat = evalR{sgmInd,4};
%     hParam.Regression_automatedSelection.groupbased = true;
%     hParam.Regression_automatedSelection.grouping = 'count';
%     hParam.Regression_automatedSelection.evaluator = 'PLSR';
%     hParam.Regression_automatedSelection.nCompPLSR = '19';
%       hParam.Regression_automatedSelection.criterion = 'MinOneStd+OptNComp';

    % change plsr components
    hParam.Regression_plsr.nComp = int32(evalR{sgmInd,5});

    % change k-Fold Validation
    % hParam.Validation_kfold.groupbased = true;
    % hParam.Validation_kfold.grouping = 'count';
    % hParam.Validation_kfold.folds = 10;
    % hParam.Validation_kfold.trainAlways = {};
    % hParam.Validation_kfold.iterations = 1;

    % change k-Fold Testing
%     hParam.Testing_kfold.groupbased = true;
%     hParam.Testing_kfold.grouping = 'count';
%     hParam.Testing_kfold.folds = 10;
%     hParam.Testing_kfold.trainAlways = {};
%     hParam.Testing_kfold.iterations = 1;

    % change groups Testing
    % grouping = project.getGroupingByCaption('grouping');
    %Param.Testing_groups.grouping = 'grouping';
    %hParam.Testing_groups.groups = {'2','3','5','6','7','10','11','12','13','14'};    
    
    % change groups Testing
    % grouping = project.getGroupingByCaption('grouping')
    % hParam.Testing_groups.grouping = 'rand_test';
    % hParam.Testing_groups.groups = {'1'};
    
    % hParam.DataReduction_useOnlyDataInGroup.grouping = 'rand_train_rand_test';
    

    % set hyper parameter
    project.currentModel.processingChain.setChainHyperParameters(hParam);
  
    model.train(data);
    modelStruct = struct(model);
    % modelStruct = rmfield(modelStruct,'processingChain');
    % modelStruct = rmfield(modelStruct,'fullModelData');
    modelStruct.datas = struct(model.datas(1));
    
    % save for evaluation 1: MinOneStd
  
    eval.nFeatMinOneStd = evalR{sgmInd,4};
    eval.nCompMinOneStd = evalR{sgmInd,5};
    
    eval.err.training = sqrt(mean((modelStruct.datas.trainingPrediction(modelStruct.datas.trainingSelection)-modelStruct.datas.target(modelStruct.datas.trainingSelection)).^2));
    eval.err.testing = sqrt(mean((modelStruct.datas.testingPrediction(modelStruct.datas.testingSelection)-modelStruct.datas.target(modelStruct.datas.testingSelection)).^2));
    eval.fullModelData= modelStruct.datas;
    
    eval.pred.trainData = modelStruct.datas.target(modelStruct.datas.trainingSelection);
    eval.pred.trainPred = modelStruct.datas.trainingPrediction(modelStruct.datas.trainingSelection);
    eval.pred.testPred = modelStruct.datas.testingPrediction(modelStruct.datas.testingSelection);
    eval.pred.testData = modelStruct.datas.target(modelStruct.datas.testingSelection);
    
    eval.pred.trainOffset = modelStruct.datas.offsets(modelStruct.datas.trainingSelection);
    eval.pred.testOffset = modelStruct.datas.offsets(modelStruct.datas.testingSelection);
%      % save for evaluation 2: Elbow
%     valErr = eval.err.validation(end,:);
%     x = 1:1:(length(valErr));
%     y = valErr;
%     p1 = [x(1),y(1)];
%     p2 = [x(end),y(end)];
%     dpx = p2(1) - p1(1);
%     dpy = p2(2) - p1(2);
%     dp = sqrt(sum((p2-p1).^2));
%     dists = abs(dpy*x - dpx*y + p2(1)*p1(2) - p2(2)*p1(1)) / dp;
%     [~,idx] = max(dists);
% 
%     eval.errorTestElbow1 = eval.err.testing(end,idx);
%     eval.nFeatElbow1 = idx;
%     eval.nCompElbow1 = size(eval.err.validation,1);
%     
%      % save for evaluation 3: Double Elbow
%     valErr = eval.err.validation(:,idx)';
%     x = 1:1:(length(valErr));
%     y = valErr;
%     p1 = [x(1),y(1)];
%     p2 = [x(end),y(end)];
%     dpx = p2(1) - p1(1);
%     dpy = p2(2) - p1(2);
%     dp = sqrt(sum((p2-p1).^2));
%     dists = abs(dpy*x - dpx*y + p2(1)*p1(2) - p2(2)*p1(1)) / dp;
%     [~,idx2] = max(dists);
%     
%     eval.errorTestElbow2 = eval.err.testing(idx2,idx);
%     eval.nFeatElbow2 = idx;
%     eval.nCompElbow2 = idx2;
    
     % save for evaluation 4: plot selected features
     % select the correct cluster for sensor 
     for i=1:length(project.clusters)
        clust=project.clusters(1,i);
        if ~isempty(clust.featureData)
            cluster = clust;
        end
     end
     
     % select sensor and save the important datas
    cluster = struct(cluster);
    clusterStruct.nCyclePoints = cluster.nCyclePoints;
    clusterStruct.samplingPeriod = cluster.samplingPeriod;
    
    for i=1:length(cluster.sensors)
        getGroup=single(cluster.featureData.groupings);
        setGroup=getGroup(:,1);
        setGroup(~isnan(setGroup))=1;
        setGroup(isnan(setGroup))=0;
        setGroup=logical(setGroup);
    
        allData=cluster.sensors(1, i).data;
        redDat=allData(setGroup,:);
        clusterStruct.sensors{i}.data=redDat(1,:); 
        clusterStruct.sensors{i}.caption=cluster.sensors(1, i).caption;
        
        dpbmean=cluster.sensors(1, i).featureDefinitionSet.featureDefinitions.getByCaption('mean').dataProcessingBlock;
        iPosmean=dpbmean.parameters.getByCaption('iPos').value;
        clusterStruct.sensors{i}.iPosmean=iPosmean;
        
        dpbpoly=cluster.sensors(1, i).featureDefinitionSet.featureDefinitions.getByCaption('polyfit').dataProcessingBlock;
        iPospoly=dpbpoly.parameters.getByCaption('iPos').value;
        clusterStruct.sensors{i}.iPospoly=iPospoly;
        
    end
    
    eval.currentCluster = clusterStruct ;
    eval.featCap = project.mergedFeatureData.featureCaptions';
    eval.rank = evalR{sgmInd, 15};
    
    for i=2:length(model.datas)
        modelStruct.datas(i) = struct(model.datas(i));
        modelStruct.datas(i).data = [];
    end
end

function setSensorActive(project, cluster, sensor)
    clusters = project.clusters;
    captions = clusters.getCaption();
    c = clusters(ismember(captions,cluster));
    sensors = c.sensors();
    captions = sensors.getCaption();
    s = sensors(ismember(captions,sensor));
    s.setActive(true);
end

function m = getModel(project, model)
    models = project.models;
    captions = models.getCaption();
    m = models(ismember(captions,model));
end

% function p = getByCaption(objArray,caption)
%     p = objArray(objArray.getCaption() == caption);
% end