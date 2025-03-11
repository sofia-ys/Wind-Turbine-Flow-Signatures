clear all;
close all;


frequencyLimit = 11.4 / (126/80 * 2);   % Velocity over waveLength

correlartioThreshhold = 0.1;

plotProbes = 0;

saveFigureProbe = 1;

saveFigureTurb = 0;

saveIndex = 0;

% fileSaving = "NREL_FXXXXX_XD_XX_00";
% 
% fileNmaeProbe = "./0309/"+ fileSaving +"/U_Inflow.txt";  


filenameSavingArray = {};

filenameSavingArray{1} = "NREL_FXXXXX_5D_000_00025_copy";
filenameSavingArray{2} = "NREL_FXXXXX_5D_000_00025_copy";
filenameSavingArray{3} = "NREL_FXXXXX_5D_000_00025_copy";

TI_Array = {};
TI_Array{1} = "TI $= ??.?\%$";
TI_Array{2} = "TI $= ??.?\%$";
TI_Array{3} = "TI $= ??.?\%$";


figureSizeSet = [100 500 1800 430];
FigMaster = figure('Renderer','painters','Position',figureSizeSet);
anno_position = [-0.7, 0.93];
Fig1 = tiledlayout(1,4,"TileSpacing","compact");


mean_TI = {};
sigma_TI = {};

autoCorrelationStore = {};
integralLengthScaleStore = {};
integralLengthScaleMeanStore = {};

integralLengthScaleAreaStore = {};
integralLengthScaleMeanAreaStore = {};

spectrumStore = {};
spectrumMeanStore = {};


for tttt = 1:3

    %% Read Probes


    load("./" + filenameSavingArray{tttt} + "/probedDataMat.mat");


    uToAnalyze = probedDataMat.uToAnalyze;
    vToAnalyze = probedDataMat.vToAnalyze;
    wToAnalyze = probedDataMat.wToAnalyze;

    tToAnalyze = probedDataMat.tToAnalyze;

    deltaT = tToAnalyze(101) - tToAnalyze(100);

    porbesNum = size(uToAnalyze, 2);

    probefocused = 0;

    %%

    freeStream = 11.4;

    U_mean = mean(uToAnalyze, 1).';
    V_mean = mean(vToAnalyze, 1).';
    W_mean = mean(wToAnalyze, 1).';


    U_std = std(uToAnalyze, 1).';
    V_std = std(vToAnalyze, 1).';
    W_std = std(wToAnalyze, 1).';

    % turbulenceIntensityWRONG = 1/3 * (U_std/freeStream + V_std/freeStream  + W_std/freeStream ) * 100;

    turbulenceIntensity = sqrt(1/3 * (U_std.^2 + V_std.^2  + W_std.^2 ) ) / freeStream * 100;

    sigma_TI{tttt} = std(turbulenceIntensity);
    mean_TI{tttt} = mean(turbulenceIntensity);

%     string(mean(U_mean/11.4))) 
    sprintf('U_mean = %.3f', mean(U_mean/11.4) )

    fprintf("\n");
    fprintf("\n");

% sprintf('%.0f', mean(turbulenceIntensity))

    fprintf("TI = " + sprintf('%.2f', mean(turbulenceIntensity)) + "%%"  );

    fprintf("\n");
    fprintf("\n");
 

    fprintf("sigma_u/Vo = " + sprintf('%.3f', mean(U_std/11.4) )  );

    fprintf("\n");
    fprintf("\n");

    fprintf("sigma_v/sigma_u = " + sprintf('%.3f', mean(V_std./U_std) )   );

    fprintf("\n");
    fprintf("\n");


    fprintf("sigma_w/sigma_u = " + sprintf('%.3f', mean(W_std./U_std) )   );

    fprintf("\n");
    fprintf("\n");


    u_focused = uToAnalyze(:, probefocused + 1);
    t_focused = tToAnalyze(:, probefocused + 1);

    U_std_focused = U_std(probefocused + 1);
    V_std_focused = V_std(probefocused + 1);
    W_std_focused = W_std(probefocused + 1);



    %% Autocorrelation


    integralLengthScale = 100 * ones(porbesNum, 1);
    integralLengthScaleArea = 100 * ones(porbesNum, 1);
    r = [];

    for ee = 1:porbesNum
        upp    = uToAnalyze(:, ee ) - mean(uToAnalyze(:, ee ));
        t_focusedShift = tToAnalyze(:, ee) - tToAnalyze(1, ee);
        dt = deltaT;
        nt = numel(t_focusedShift);

%         r = xcorr(upp); 
%         r = r(nt:end); 
%         r = r/r(1);

        rtemp = xcorr(upp); 
        rtemp = rtemp(nt:end);
%         rtemp(1)
        rtemp = rtemp/rtemp(1);
        
        r(:, ee) = rtemp;

        indz = find(r(:, ee) <= correlartioThreshhold);
        indz = indz(1);

        dtplot = 0.0;

        fi = 1/(t_focusedShift(indz)-dtplot);

        integralLengthScale(ee, 1) = t_focusedShift(indz) * freeStream;

        integralLengthScaleArea(ee, 1) = trapz( rtemp(1:indz) ) * deltaT * freeStream;

    end

    autoCorrelationStore{tttt} = r;
    integralLengthScaleStore{tttt} = integralLengthScale;
    integralLengthScaleMeanStore{tttt} = mean(integralLengthScale);

    integralLengthScaleAreaStore{tttt} = integralLengthScaleArea;
    integralLengthScaleMeanAreaStore{tttt} = mean(integralLengthScaleArea);



%     figure();
%     plot(autoCorrelationStore{tttt}(1:1000, :));


    %% Q7. S(f)
    % Compute and plot the 1D spectrum

    %% Plot settings


    lable_font_size = 24;
    title_font_size = 20;
    legend_font_size = 16;
    gca_font_size = 14;

    colorbar_font_size = 24;

    line_Width = 1.8;
    Marker_Size = 20.0;

    % annotation_font_size = 20.0;

    colorIndex{1} = [0.0 0.0 0.8];
    colorIndex{2} = [0.1 0.5 0.1];
    colorIndex{3} = [1.0 0.0 0.0];
    colorIndex{4} = [0.75, 0, 0.75];

    annotation_font_size = 24.0;

    D = 126;
    V0 = 11.4;

    
    S_Array = [];
    for ee = 1:porbesNum

        u_focused = uToAnalyze(:, ee );

        N = length(u_focused);
        f_N = 1/(2*deltaT);
        df = f_N/(N/2);
        f = 0:df:f_N;
        U = fft(u_focused.')/N;                      % Discrete Fourier transform (m/s)
        A = abs([U(1) 2*U(2:N/2) U(N/2+1)]);        % Harm amps (m/s)


        S = 0.5*A.^2/df;  % Energy density (m^2/s^2/Hz=m^2/s)

        S_Array(:, ee) = S;

    end

    spectrumStore{tttt} = S_Array;

    spectrumMeanStore{tttt} = mean(S_Array, 2);

    nexttile(tttt + 1);

    indexlokking = find(f <= frequencyLimit);


    loglog(f(indexlokking),spectrumMeanStore{tttt}(indexlokking));
    set(gca, 'FontSize', gca_font_size)

    xlabel("$f$~[Hz]", 'Interpreter','latex','FontSize',lable_font_size, 'fontWeight','bold');

    if tttt == 1
        ylabel("$S_u(f)$~[m$^2$/s]", 'Interpreter','latex','FontSize',lable_font_size, 'fontWeight','bold');
    end
    hold on;

    nw = 3;
    grid on;



    ftoPlot = f(3:end).';

    indexlokking = find(f <= 2 * frequencyLimit);

    % plot(ftoPlot,3e-7.*ftoPlot.^(-5/3),'k-','Linewidth',2 );

%     KamialSpec = 0.05 * mean(U_std)^2 * (mean(integralLengthScale)/freeStream)^(-2/3) * ftoPlot.^(-5/3);
    KamialSpec = 0.05 * mean(U_std)^2 * (mean(integralLengthScaleArea)/freeStream)^(-2/3) * ftoPlot.^(-5/3);


    plot(ftoPlot(indexlokking), KamialSpec(indexlokking),'k-','Linewidth',2 );

    % plot([frequencyLimit frequencyLimit], [0 1000], 'Linewidth', 0.5 );


    if tttt == 1
        legend('Turbulence spectrum', 'Kaimal IEC', 'Interpreter','Latex', 'fontsize',legend_font_size, 'Location', 'SouthWest');

    end


    xlim([1e-3 1e1]);
    ylim([1e-6 1e2]);


    if tttt == 1
        h=text(1,1, "\textbf{(b)}");
        h.Position = [1.5e-4 2.7e2];  %[-3.5, 1.15]
        h.Interpreter = 'latex';
        h.FontSize = annotation_font_size;
    end

    h3=text(1,1, "TI $= " + sprintf('%.2f', mean(turbulenceIntensity)) + "\%$");
    h3.Position = [2.4e-2 2.7e2];  %[-3.5, 1.15]
    h3.Interpreter = 'latex';
    h3.FontSize = 24;

end




%% Plot settings


nexttile(1);

legendPlot = 1;
legendTransparent = 1;


lable_font_size = 24;
title_font_size = 20;
legend_font_size = 16;
gca_font_size = 14;

colorbar_font_size = 24;

line_Width = 1.8;
Marker_Size = 20.0;

% annotation_font_size = 20.0;

colorIndex{1} = [0.0 0.0 0.8];
colorIndex{2} = [0.1 0.5 0.1];
colorIndex{3} = [1.0 0.0 0.0];
colorIndex{4} = [0.75, 0, 0.75];

annotation_font_size = 24.0;

D = 126;
V0 = 11.4;

%

hold on;

disky = [ (-1:0.01:0.99), (1.00:-0.01:-1.00)];
diskz = sqrt(1 - disky.^2);

diskz(1, 201:end) = -diskz(1, 201:end);


scatter( probedDataMat.probeYLocations / D, probedDataMat.probeZLocations / D, 100,'X', 'k', 'LineWidth', 2);
plot( disky/2, diskz/2,  'LineWidth', 1.2, 'Color', 'k', 'LineStyle','--');

axis equal;
set(gca, 'FontSize', gca_font_size)



ylabel('$z/D$~[-]','Interpreter','latex','FontSize',lable_font_size,'fontWeight','bold');
xlabel('$y/D$~[-]','Interpreter','latex','FontSize',lable_font_size,'fontWeight','bold');
xtickformat('%,.1f')
ytickformat('%,.1f')


xticks(-0.8:0.4:0.8);
yticks(-0.8:0.4:0.8);

ylim([-0.8, 0.80]);
xlim([-0.8, 0.80]);
grid on;
box on;

h2=text(1,1, "\textbf{(a)}");
h2.Position = [-1.1 0.88];  %[-3.5, 1.15]
h2.Interpreter = 'latex';
h2.FontSize = annotation_font_size;


localName = "./" + "paper1_V2_freqMoreProbe";

if  saveFigureTurb == 1
    exportgraphics(FigMaster,localName + ".png")
    exportgraphics(FigMaster,localName + ".eps")
end



saveInformation.KamialSpec = KamialSpec;
saveInformation.f = f;
saveInformation.S = S;
% close all;

% saveSaveSave = "./inflowTurb/" + fileSaving;


disp("mean_TI")
disp("")
disp(mean_TI);

disp("sigma_TI")
disp("")
disp(sigma_TI);

disp("Integral Length Scale")
disp("")
disp(integralLengthScaleMeanStore);

disp("STD Integral Length Scale")
disp("")
disp(std(integralLengthScaleStore{1}));
disp(std(integralLengthScaleStore{2}));
disp(std(integralLengthScaleStore{3}));


disp("Integral Length Scale (Area)")
disp("")
disp(integralLengthScaleMeanAreaStore);


disp("STD Integral Length Scale (Area)")
disp("")
disp(std(integralLengthScaleAreaStore{1}));
disp(std(integralLengthScaleAreaStore{2}));
disp(std(integralLengthScaleAreaStore{3}));




disp("sigma_TI/TI")

for tt = 1:3
    fprintf('%.3f', (sigma_TI{tt}/mean_TI{tt}) );
    fprintf("\n");
end

disp("sigma_Lu/Lu")

for tt = 1:3
    fprintf('%.2f', (std(integralLengthScaleStore{tt})/integralLengthScaleMeanStore{tt}) );
    fprintf("\n");
end


disp("sigma_Lu/Lu (Area)")

for tt = 1:3
    fprintf('%.2f', (std(integralLengthScaleAreaStore{tt})/integralLengthScaleMeanAreaStore{tt}) );
    fprintf("\n");
end


% probedDataMat.uToAnalyze = uToAnalyze
% probedDataMat.vToAnalyze = vToAnalyze
% probedDataMat.wToAnalyze = wToAnalyze
% probedDataMat.tToAnalyze = tToAnalyze
% probedDataMat.probeYLocations = probeLocations.y(probeInterested)
% probedDataMat.probeZLocations = probeLocations.z(probeInterested)
% probedDataMat.probeXLocations = probeLocations.x(probeInterested)
% save("./NREL_FXXXXX_5D_000_00025_copy/probedDataMat.mat", "probedDataMat")



