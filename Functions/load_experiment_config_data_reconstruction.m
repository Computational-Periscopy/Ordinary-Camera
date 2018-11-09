%% This script hold the measurement configurations of the experimental tests
% performed. Given the experiment test letter, the correct dataset will be
% loaded for MATLAB processing.

% John Murray-Bruce at Boston University


switch TestLetter
    
    case 'D11'
        filename = 'TestPosD11';
        subblocksperaxis = [1 1];
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36];
        D = 1.03;
        
        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            %% Estimated using occluder localization scripts
            % Projection for mushroom
            Occ_LLcorner = [0.4583 0.5408 0.2026];
            
            %Projection for Tommy
            %Occ_LLcorner = [0.4569 0.5744 0.2080];
            
            % Projection for BU red (bur)
            %Occ_LLcorner = [0.4733   0.5661    0.2072];
          
            % Projection for RGB bars (colbar)
            %Occ_LLcorner = [0.4693 0.5629 0.2080];
            
            %%
        else
            % True
            Occ_LLcorner = [0.475 D-0.460 0.214];
        end
        
        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];
        
        %CAMERA CONFIG
        FOV_size = [0.4372 0.4372];
        FOV_LLCorner = [0.521 0.048];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];
        
        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];
        
        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024];
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);
        
        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));
        
        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];
       
        
        
    case 'D11Video'
        filename = 'TestPosD11Video';
        subblocksperaxis = [1 1];
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36];
        D = 1.03;
        
        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            Occ_LLcorner = [0.4750    0.5583    0.2000];
        else
            Occ_LLcorner = [0.477 D-0.466 0.21];
        end
        
        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];
        
        %CAMERA CONFIG
        FOV_size = [0.4672 0.4672];
        FOV_LLCorner = [0.496 0.014];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];
        
        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];
        
        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024]; %[hor, ver] pixels
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);
        
        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));
        
        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];
        
        
case {'Chair'}
        filename = 'TestPosChair';
        subblocksperaxis = [1 1];
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36];
        D = 1.03;

        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            Occ_LLcorner = [0.503 D-0.436 0.244];
        else
            Occ_LLcorner = [0.506 D-0.436 0.242];
        end

        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];

        %CAMERA CONFIG
        FOV_size = [0.401 0.401];
        FOV_LLCorner = [0.543 0.064];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];

        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];

        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024]; %[hor, ver] pixels
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);

        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));

        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];

    case {'ChairRealScene'}
        filename = 'TestPosChair';
        subblocksperaxis = [1 ,1];
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36].*subblocksperaxis;
        D = 1.03;

        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            Occ_LLcorner = [0.498 D-0.415 0.244];
        else

            Occ_LLcorner = [0.506 D-0.436 0.242];
        end

        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];

        %CAMERA CONFIG
        FOV_size = [0.401 0.401];
        FOV_LLCorner = [0.543 0.064];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];

        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];

        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024]; %[hor, ver] pixels
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);

        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));

        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];
 
    
    case {'ChairAmb'}
        filename = 'TestPosChair';
        subblocksperaxis = [1 ,1]; %Change to [2,2] to double resolution
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36].*subblocksperaxis;
        D = 1.03-0.006;

        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            Occ_LLcorner = [0.498 D-0.415 0.244];
        else

            Occ_LLcorner = [0.506 D-0.436 0.242];
        end

        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];

        %CAMERA CONFIG
        FOV_size = [0.401 0.401];
        FOV_LLCorner = [0.543 0.064];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];

        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];

        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024]; %[hor, ver] pixels
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);

        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));

        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];
    
        
 case {'Chair3DScene'}
        filename = 'TestPosChair';
        subblocksperaxis = [1 , 1]; %Change to [2,2] to double resolution
        NumBlocks_sim = [29 36].*subblocksperaxis;
        NumBlocks_cal = [29 36].*subblocksperaxis;
        D = 1.03;

        % OCCLUDER DEFINITION
        Occ_size = [0.075 0 0.075];
        if useEstimatedOccPos
            Occ_LLcorner = [0.498 D-0.415 0.244];
        else

            Occ_LLcorner = [0.506 D-0.436 0.242];
        end

        Occluder = [Occ_LLcorner; Occ_LLcorner + Occ_size];

        %CAMERA CONFIG
        FOV_size = [0.40 0.40];
        FOV_LLCorner = [0.557 0.062];
        FOV_cord = [FOV_LLCorner; FOV_LLCorner + FOV_size];

        % MONITOR CONFIG
        Mon_Offset = [0.0210 0.136];

        NumBlocks_col = NumBlocks_cal(2);
        NumBlocks_row = NumBlocks_cal(1);
        ScreenSize = [0.408 0.3085];
        ScreenResolution = [1280 1024]; %[hor, ver] pixels
        NumPixelsPerMonitorBlock = 35;
        PixelSize_m = (ScreenSize./ScreenResolution);

        Mon_Offset(2) = Mon_Offset(2) + PixelSize_m(2)*mod(ScreenResolution(2),NumBlocks_cal(1));
        Mon_Offset(1) = Mon_Offset(1) + 1*PixelSize_m(1)*mod(ScreenResolution(1),NumBlocks_cal(2));

        IlluminationBlock_Size = PixelSize_m.*NumPixelsPerMonitorBlock./[subblocksperaxis(2) subblocksperaxis(1)];
 
        
    otherwise
        disp('Not available');
end


% Set simulation of forward model parameters!

simuParams.NumBlocks = NumBlocks_sim;
simuParams.Ndiscr_mon = Ndiscr_mon;
simuParams.numPixels = numPixels;
simuParams.D = D;
simuParams.FOV_cord = FOV_cord;
simuParams.Occluder = Occluder;
simuParams.viewAngleCorrection = viewAngleCorrection;

simuParams.IlluminationBlock_Size = IlluminationBlock_Size;
simuParams.Mon_Offset = Mon_Offset;
calibParams.scaling = 1;

if ismac
    calibParams.filepath = ['', filename,'/'];
elseif ispc
    calibParams.filepath = ['', filename,'\'];
end

calibParams.x_max = NumBlocks_cal(2);
calibParams.y_max = NumBlocks_cal(1);

