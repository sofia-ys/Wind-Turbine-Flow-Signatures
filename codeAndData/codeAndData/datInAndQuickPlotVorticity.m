close all;
clear all;

timeStepSize = 0.013778915 * 15;

tic

saveFigure = 1;

closeAllIndex = 1;

deltaGrid = 1/80;

%% Which cases to compare, basic settings

D = 126;
V0 = 11.4;


figureFileNameMaster = "./NREL_FXXXXX_5D_000_00025_copy/exportedCSV/";

fileNumber = 1;

figNumber = fileNumber;


%% set up coordinate

meshGrid_xMin = -2 + deltaGrid/2;
meshGrid_xMax = 8 + deltaGrid/2;

meshGrid_zMin = -1 + deltaGrid/2;
meshGrid_zMax = 1 + deltaGrid - deltaGrid/2;

[xq,zq] = meshgrid( ...
    meshGrid_xMin:deltaGrid:meshGrid_xMax, ...
    meshGrid_zMin:deltaGrid:meshGrid_zMax);


coordinate_yPlane.xq = xq;
coordinate_yPlane.zq = zq;



save("./NREL_FXXXXX_5D_000_00025_copy/exportedMatMean/coordinate_yPlane.mat", ...
        "coordinate_yPlane");


%% interp data to uniform grid

% indexOfCSVfileArray = 0:1:2400;
indexOfCSVfileArray = 0:1:2400;

% switch for if always plot and save
alwayPlotting = 1;

for indexCSVHolder = indexOfCSVfileArray

    indexOfCSVfile = indexCSVHolder;
    indexOfCSVfileString = sprintf( '%04d', indexOfCSVfile );

    csvTail = {};
    csvTail{1} = "./NREL_FXXXXX_5D_000_00025_copy/exportedCSV/exportedCSV_y00D/yPlane_" + indexOfCSVfileString + ".csv";
    csvTail{2} = "./NREL_FXXXXX_5D_000_00025_copy/exportedCSV/exportedCSV_x02D/xPlane_2D_" + indexOfCSVfileString + ".csv";
    csvTail{3}=  "./NREL_FXXXXX_5D_000_00025_copy/exportedCSV/exportedCSV_x05D/xPlane_5D_" + indexOfCSVfileString + ".csv";
    csvTail{4}=  "./NREL_FXXXXX_5D_000_00025_copy/exportedCSV/exportedCSV_x08D/xPlane_8D_" + indexOfCSVfileString + ".csv";
    numberOfPlanes = length(csvTail);

    fileName_All = cell(fileNumber,numberOfPlanes);

    for pp = 1:fileNumber
        for gg = 1:numberOfPlanes
            fileName_All{pp, gg} = csvTail{gg};
        end
    end


    %% Read in all the needed data

    dataIn = cell(fileNumber,numberOfPlanes);


    for tt = 1:fileNumber

        for gg = 1:numberOfPlanes
            dataIn{tt, gg} = readtable(fileName_All{tt, gg}, 'VariableNamingRule', 'preserve');
        end

    end


    openFaomU = cell(fileNumber,numberOfPlanes);
    openFaomV = cell(fileNumber,numberOfPlanes);
    openFaomW = cell(fileNumber,numberOfPlanes);


    openFaomX = cell(fileNumber,numberOfPlanes);
    openFaomY = cell(fileNumber,numberOfPlanes);
    openFaomZ = cell(fileNumber,numberOfPlanes);



    for pp = 1:fileNumber

        for gg = 1:numberOfPlanes

            openFaomU{pp, gg} = dataIn{pp, gg}.("vorticityField:0");
            openFaomV{pp, gg} = dataIn{pp, gg}.("vorticityField:1");
            openFaomW{pp, gg} = dataIn{pp, gg}.("vorticityField:2");

            openFaomX{pp, gg} = dataIn{pp, gg}.("Points:0");
            openFaomY{pp, gg} = dataIn{pp, gg}.("Points:1");
            openFaomZ{pp, gg} = dataIn{pp, gg}.("Points:2");

        end


    end


    clear dataIn;

    %% yPlane

   [xq,zq] = meshgrid( ...
        meshGrid_xMin:deltaGrid:meshGrid_xMax, ...
        meshGrid_zMin:deltaGrid:meshGrid_zMax);

    indexOfPlane = 1;

    wq = griddata(openFaomX{pp, indexOfPlane}/D, openFaomZ{pp, indexOfPlane}/D, ...
        openFaomW{pp, indexOfPlane},xq,zq);

    vq = griddata(openFaomX{pp, indexOfPlane}/D, openFaomZ{pp, indexOfPlane}/D, ...
        openFaomV{pp, indexOfPlane},xq,zq);

    uq = griddata(openFaomX{pp, indexOfPlane}/D, openFaomZ{pp, indexOfPlane}/D, ...
        openFaomU{pp, indexOfPlane},xq,zq);

    interporlatedData.uq = uq;
    interporlatedData.vq = vq;
    interporlatedData.wq = wq;

    save("./NREL_FXXXXX_5D_000_00025_copy/exportedMatMean/yPlane_" + indexOfCSVfileString + ".mat", ...
        "interporlatedData");

    toc


    %% Plotting

    % Plotting switch
    % always plotting or plot when the last Index
    if (alwayPlotting == 1 || indexCSVHolder == indexOfCSVfileArray(end))
        %% Plot settings

        %close all;

        lable_font_size = 32;
        title_font_size = 20;
        legend_font_size = 16;
        gca_font_size = 20;

        colorbar_font_size = 38;
        colorbarTicks_fontSize = 20;

        line_Width = 1.8;
        Marker_Size = 20.0;

        annotation_font_size = 32.0;

        colobar_string = "$\omega$";


        annotationPosition =  [34.6, 4.01];


        maxVel = 0.5;
        middleVel = 1;
        minVel = 0.0;
        NUMBER_LEVELS = 16;


        middleIndex = NUMBER_LEVELS * (middleVel - minVel) / (maxVel - minVel);
        levelListIcrementLower = (middleVel - minVel) / (middleIndex);
        levelListIcrementHigher = (maxVel - middleVel) / (NUMBER_LEVELS-middleIndex);

        levelListSetLower = (minVel + levelListIcrementLower):levelListIcrementLower:middleVel;
        levelListSetHigher = middleVel:levelListIcrementHigher:(maxVel - levelListIcrementHigher);

        levelListSet = [levelListSetLower(1:end-1), levelListSetHigher(1:end)];
        levelListSet = [-100.0, levelListSet, 100.0];


        [COLORMAP_u] = slanCM('jet', NUMBER_LEVELS);


        figureSizeSet = [0 -300 1550 750];
        FigMaster = figure('Visible', 'off');
        FigMaster.Position = figureSizeSet;


        arrowScale = 0.3;
        arrowTipAngle = 20;
        arrowBackScale = 0.35;
        arrowbackLength = arrowScale * arrowBackScale;

        xplaneIntersted = [2 5 8];

        %% yPlane

        indexOfPlane = 1;

        meshGrid_xMin = -2 + deltaGrid/2;
        meshGrid_xMax = 8 + deltaGrid/2;

        meshGrid_zMin = -1 + deltaGrid/2;
        meshGrid_zMax = 1 + deltaGrid - deltaGrid/2;

        [xq,zq] = meshgrid( ...
            meshGrid_xMin:deltaGrid:meshGrid_xMax, ...
            meshGrid_zMin:deltaGrid:meshGrid_zMax);


        wq = griddata(openFaomX{pp, indexOfPlane}/D, openFaomZ{pp, indexOfPlane}/D, ...
            openFaomW{pp, indexOfPlane}/V0,xq,zq);

        uq = griddata(openFaomX{pp, indexOfPlane}/D, openFaomZ{pp, indexOfPlane}/D, ...
            openFaomU{pp, indexOfPlane}/V0,xq,zq);

        C = surf(squeeze(xq), squeeze(zq) .* 0, ...
            squeeze(zq), squeeze(uq));
        C.EdgeColor = 'none';
        hold on;

        x_End = 8;

        y_Left = -1;
        y_Right = 1;

        quiverDelta = 1/4;
        quiverClear1 = 1*quiverDelta;
        quiverClear2 = 1*quiverDelta;

        quiverFactor = 12;

        [quiverXq,quiverZq] = meshgrid( ...
            (-2 + 1/4):(quiverDelta*2):(8 - 1/4), ...
            (-1 + quiverClear2):quiverDelta:(1 - quiverClear2) ...
            );

        quiverU = griddata(openFaomX{pp, indexOfPlane}/D, openFaomZ{pp, indexOfPlane}/D, ...
            openFaomU{pp, indexOfPlane}/V0 * quiverFactor,quiverXq,quiverZq);
        quiverW = griddata(openFaomX{pp, indexOfPlane}/D, openFaomZ{pp, indexOfPlane}/D, ...
            openFaomW{pp, indexOfPlane}/V0 * quiverFactor,quiverXq,quiverZq);

        quiverU = quiverU ./ sqrt( norm(quiverU + quiverW) );
        quiverW = quiverW ./ sqrt( norm(quiverU + quiverW) );


        quiverYq = quiverXq * 0.0;
        quiverV = quiverU * 0.0;
        %     view(0, 0);

        for ii = 1:size(quiverXq,1)
            for jj = 1:size(quiverXq,2)


                nowX = quiverXq(ii, jj);
                nowY = quiverYq(ii, jj) - 1/120;
                nowZ = quiverZq(ii, jj);

                arrowStart = [nowX nowY nowZ];
                arrowDir = [quiverU(ii, jj) quiverV(ii, jj) quiverW(ii, jj)];
                arrowEnd = arrowStart + arrowScale * arrowDir;

                %             arrow('Start', arrowStart, 'Stop', arrowEnd, 'NormalDir', [0 -1 0], ...
                %                 'Length', (15 * norm(arrowEnd - arrowStart)) , ...
                %                 'Width', (5 * norm(arrowEnd - arrowStart)), ...
                %                 'BaseAngle', 90, 'TipAngle', 25);
                plot3( [arrowStart(1)  arrowEnd(1)], ...
                    [arrowStart(2)  arrowEnd(2)], ...
                    [arrowStart(3)  arrowEnd(3)], ...
                    'Color', 'k', 'LineWidth', norm(arrowDir) * 0.8);



                arrowBackDir_11 = cosd(arrowTipAngle)*arrowDir(1) - sind(arrowTipAngle)*arrowDir(3);
                arrowBackDir_13 = sind(arrowTipAngle)*arrowDir(1) + cosd(arrowTipAngle)*arrowDir(3);

                arrowBackPos_11 = arrowEnd(1) - arrowBackDir_11*arrowbackLength;
                arrowBackPos_13 = arrowEnd(3) - arrowBackDir_13*arrowbackLength;

                arrowBackDir_21 = cosd(-arrowTipAngle)*arrowDir(1) - sind(-arrowTipAngle)*arrowDir(3);
                arrowBackDir_23 = sind(-arrowTipAngle)*arrowDir(1) + cosd(-arrowTipAngle)*arrowDir(3);

                arrowBackPos_21 = arrowEnd(1) - arrowBackDir_21*arrowbackLength;
                arrowBackPos_23 = arrowEnd(3) - arrowBackDir_23*arrowbackLength;

                plot3( [arrowBackPos_21 arrowEnd(1) arrowBackPos_11 ], ...
                    [arrowEnd(2) arrowEnd(2)  arrowEnd(2) ], ...
                    [arrowBackPos_23 arrowEnd(3) arrowBackPos_13 ], ...
                    'Color', 'k', 'LineWidth', norm(arrowDir) * 0.8);


            end
        end


        %     view(0, 0);
        %     scaleFactor = 0.1;
        %     quiverHolder = quiver3(quiverXq-8,quiverYq - 0.01,quiverZq, ...
        %         quiverU .* scaleFactor,quiverV .* scaleFactor,quiverW .* scaleFactor, '-', ...
        %         'Color', 'k', 'LineWidth', 2);
        %     %     quiverHolder.AutoScale = "off";
        %     quiverHolder.ShowArrowHead = 'on';


        for kk = -1:1:7
            plot3([kk kk], [0 0], [-1.00 -0.90], '-k', 'LineWidth', 1.0)
            plot3([kk kk], [0 0], [1 0.9], '-k', 'LineWidth', 1.0)
        end


        for kk = -0.5:0.5:0.5
            plot3([-2.00 -1.90], [0 0], [kk kk], '-k', 'LineWidth', 1.0)
            plot3([8 7.9], [0 0], [kk kk], '-k', 'LineWidth', 1.0)
        end

        plot3( [-2 8 8 -2 -2], [0 0 0 0 0], [-1 -1 1 1 -1], '-k', 'LineWidth', 0.5 )


        axis equal;

        caxis([minVel maxVel]);
        colormap(COLORMAP_u);


        plot3([0 0], [0 0], [-0.5 0.5], '-k', 'LineWidth', 4)

         %% xplanes

        for xxx = xplaneIntersted

            indexOfPlane = indexOfPlane + 1;

            deltaGrid = 1/80;

            meshGrid_yMin = -1 + deltaGrid/2;
            meshGrid_yMax = 1 - deltaGrid/2;

            meshGrid_zMin = -1 + deltaGrid/2;
            meshGrid_zMax = 1 - deltaGrid/2;

            [yq,zq] = meshgrid( ...
                meshGrid_yMin:deltaGrid:meshGrid_yMax, ...
                meshGrid_zMin:deltaGrid:meshGrid_zMax);


            wq = griddata(openFaomY{pp, indexOfPlane}/D, openFaomZ{pp, indexOfPlane}/D, ...
                openFaomW{pp, indexOfPlane},yq,zq);

            uq = griddata(openFaomY{pp, indexOfPlane}/D, openFaomZ{pp, indexOfPlane}/D, ...
                openFaomU{pp, indexOfPlane},yq,zq);

            C2 = surf(squeeze(yq) .* 0 + xxx, squeeze(yq), ...
                squeeze(zq), squeeze(uq));
            C2.EdgeColor = 'none';

            caxis([minVel maxVel]);
            colormap(COLORMAP_u);


            quiverDelta = 1/4;
            quiverClear1 = 1*quiverDelta;
            quiverClear2 = 1*quiverDelta;


            quiverFactorXplane = quiverFactor;

            [quiverYq,quiverZq] = meshgrid( ...
                (-1 + quiverClear1):(quiverDelta):(1 - quiverClear1), ...
                (-1 + quiverClear2):(quiverDelta):(1 - quiverClear2) ...
                );

            quiverV = griddata(openFaomY{pp, indexOfPlane}/D, openFaomZ{pp, indexOfPlane}/D, ...
                openFaomV{pp, indexOfPlane}/V0 * quiverFactorXplane,quiverYq,quiverZq);
            quiverW = griddata(openFaomY{pp, indexOfPlane}/D, openFaomZ{pp, indexOfPlane}/D, ...
                openFaomW{pp, indexOfPlane}/V0 * quiverFactorXplane,quiverYq,quiverZq);

            quiverV = quiverV ./ sqrt( norm(quiverV + quiverW) );
            quiverW = quiverW ./ sqrt( norm(quiverV + quiverW) );

            quiverXq = quiverYq .* 0.0 + xxx;
            quiverU = quiverV .* 0.0;

            %     quiverW = quiverW + 1 * quiverFactor;


            for ii = 1:size(quiverXq,1)
                for jj = 1:size(quiverXq,2)

                    nowX = quiverXq(ii, jj) - 1/120;
                    nowY = quiverYq(ii, jj);
                    nowZ = quiverZq(ii, jj);

                    arrowStart = [nowX nowY nowZ];
                    arrowDir = [quiverU(ii, jj) quiverV(ii, jj) quiverW(ii, jj)];
                    arrowEnd = arrowStart + arrowScale * arrowDir;

                    plot3( [arrowStart(1)  arrowEnd(1)], ...
                        [arrowStart(2)  arrowEnd(2)], ...
                        [arrowStart(3)  arrowEnd(3)], ...
                        'Color', 'k', 'LineWidth', norm(arrowDir) * 0.8);

                    arrowBackDir_12 = cosd(arrowTipAngle)*arrowDir(2) - sind(arrowTipAngle)*arrowDir(3);
                    arrowBackDir_13 = sind(arrowTipAngle)*arrowDir(2) + cosd(arrowTipAngle)*arrowDir(3);

                    arrowBackPos_12 = arrowEnd(2) - arrowBackDir_12*arrowbackLength;
                    arrowBackPos_13 = arrowEnd(3) - arrowBackDir_13*arrowbackLength;

                    arrowBackDir_22 = cosd(-arrowTipAngle)*arrowDir(2) - sind(-arrowTipAngle)*arrowDir(3);
                    arrowBackDir_23 = sind(-arrowTipAngle)*arrowDir(2) + cosd(-arrowTipAngle)*arrowDir(3);

                    arrowBackPos_22 = arrowEnd(2) - arrowBackDir_22*arrowbackLength;
                    arrowBackPos_23 = arrowEnd(3) - arrowBackDir_23*arrowbackLength;


                    plot3( [arrowEnd(1) arrowEnd(1) arrowEnd(1) ], ...
                        [arrowBackPos_22 arrowEnd(2) arrowBackPos_12 ], ...
                        [arrowBackPos_23 arrowEnd(3) arrowBackPos_13 ], ...
                        'Color', 'k', 'LineWidth', norm(arrowDir) * 0.8);


                end
            end

            degForCir = 0:(2*pi/100):(2*pi);
            zCir = 0.5 * sin(degForCir);
            yCir = 0.5 * cos(degForCir);
            xCir = degForCir.*0 + xxx;

            plot3( xCir, yCir, zCir, 'Color', 'k', ...
                'LineStyle', '-', 'LineWidth', 2);

            for kk = -0.5:0.5:0.5
                plot3([xxx xxx], [kk kk], [-1.0 -0.9], '-k', 'LineWidth', 1.0)
                plot3([xxx xxx], [kk kk], [1     0.9], '-k', 'LineWidth', 1.0)
            end


            for kk = -0.5:0.5:0.5
                plot3([xxx xxx], [-1 -0.9], [kk kk], '-k', 'LineWidth', 1.0)
                plot3([xxx xxx], [ 1  0.9], [kk kk], '-k', 'LineWidth', 1.0)
            end

            plot3([xxx xxx], [0 0], [-1 1], '-k', 'LineWidth', 0.5)
            plot3([xxx xxx xxx xxx xxx], [1 1 -1 -1 1], [-1 1 1 -1 -1], '-k', 'LineWidth', 0.5)

        end
        %%


        c100 = colorbar('eastoutside');



        set(c100, 'FontSize', colorbarTicks_fontSize)

        set(gca, 'FontSize', gca_font_size)
        set(gca, 'TickLabelInterpreter', 'latex')

        c100.Label.String = colobar_string;
        c100.Label.FontSize = colorbar_font_size;
        c100.Label.Interpreter = 'latex';
        c100.Label.Rotation = 90;

        c100.TickLabelInterpreter = 'latex';


        hcolor=text(1,1, colobar_string);
        hcolor.Position = [12.3 -2 4.2];  %[8.2, -1.0];
        hcolor.Interpreter = 'latex';
        hcolor.FontSize = colorbar_font_size;


        xlabel('$x/D$','Interpreter','latex','FontSize',lable_font_size,'fontWeight','bold');
        ylabel('$y/D$','Interpreter','latex','FontSize',lable_font_size,'fontWeight','bold');
        zLabelz = zlabel('$z/D$','Interpreter','latex','FontSize',lable_font_size,'fontWeight','bold', 'Rotation',0);
        zLabelz.Position = [-3.1000    1.0    0.3072];



        annotationInformation = "$\phi_{\Omega} = " + sprintf( '%.2f', indexOfCSVfile * timeStepSize * V0*7*2/D / ( pi) ) + "\pi$";
        h=text(1,1, annotationInformation);
        h.Position = [-1 2 1.7];  %[8.2, -1.0];
        h.Interpreter = 'latex';
        h.FontSize = annotation_font_size;
        

        hTI=text(1,1, "TI $= 4.99\%$");
        hTI.Position = [-1 2 2.1];  %[8.2, -1.0];
        hTI.Interpreter = 'latex';
        hTI.FontSize = annotation_font_size;

        hLu=text(1,1, "$L_u = 31.5$~m");
        hLu.Position = [-1 2 2.5];  %[8.2, -1.0];
        hLu.Interpreter = 'latex';
        hLu.FontSize = annotation_font_size;

        xticks(-2:1:8);
        yticks(-1:0.5:1);
        zticks(-1:0.5:1);

        xtickformat('%,.1f')
        ytickformat('%,.1f')
        ztickformat('%,.1f')


        %     plot([ 0.0,   0.0],  [36/D, 336/D], '-k', 'LineWidth', 2.5);
        %     plot([ 6.0,   6.0],  [36/D, 336/D], '-k', 'LineWidth', 2.5);
        %     plot([ 12.0,  12.0], [36/D, 336/D], '-k', 'LineWidth', 2.5);
        %     plot([ 18.0,  18.0], [36/D, 336/D], '-k', 'LineWidth', 2.5);
        %     plot([ 24.0,  24.0], [36/D, 336/D], '-k', 'LineWidth', 2.5);

        xmin=-2.0;
        xmax=8;
        ymin=-1;
        ymax=1;
        zmin=-1;
        zmax=1;

        xlim([xmin xmax]);
        ylim([ymin ymax]);
        zlim([zmin zmax]);

        % boxLineWidth = 0.8;
        % plot3([-2 -2], [-2 -2], [0 2.5], '-k', 'LineWidth', boxLineWidth)
        % plot3([-2 -2], [-2 2], [2.5 2.5], '-k', 'LineWidth', boxLineWidth)
        % plot3([-2 12], [-2 -2], [2.5 2.5], '-k', 'LineWidth', boxLineWidth)

        % set(gca,'visible','off')

        box on;

        %if c100.TickLabels{7} == '1'
            %c100.TickLabels{7} = '1.0';
        %end

        % camlight;
        % camlight('forntleft')



        if saveFigure == 1
            exportFolder = "./NREL_FXXXXX_5D_000_00025_copy/exportedPNG_Vorticity";

            % Check if the folder exists, create it if it doesn't
            if ~exist(exportFolder, 'dir')
                mkdir(exportFolder);
            end
            exportgraphics(FigMaster, exportFolder + "/3D" + indexOfCSVfileString + ".png");
            % saveas(FigMaster, "./NREL_FXXXXX_5D_000_00025_copy/exportedPNG/yPlane_" + indexOfCSVfileString + ".png");

        end

    end

end


if (closeAllIndex == 1)
    close all;

    clear FigMaster;
    %clear Fig1;
end




%%

if (closeAllIndex == 1)
    close all;

    clear FigMaster;
    clear Fig1;

end


toc









