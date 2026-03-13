%% ================================================================
%  Figure Generation — Paper 2 (Observability)
%  ----------------------------------------------------------------
%  "Equivalence of Kalman and Hewer Observability for Generalised
%   Bilinear Matrix Periodic Systems"
%  Authors: Sri Venkata Durga Sudarsan Madhyannapu &
%           Sravanam Pradheep Kumar
%
%  Generates 4 publication-quality figures at 500 DPI:
%    Fig1_Obs_EigModuli.pdf
%    Fig2_Obs_GramianGrowth.pdf
%    Fig3_Obs_SpectralBound.pdf
%    Fig4_Obs_ScalingStudy.pdf
%
%  All values hardcoded from MATLAB v3 verified output.
%  ================================================================

clear; clc; close all;

% ----------------------------------------------------------------
%  VERIFIED VALUES FROM MATLAB v3 RUN
% ----------------------------------------------------------------

% Monodromy eigenvalues (same system as Paper 1)
eig_moduli_raw = [1.0000, 1.0000, 1.5883, 0.6296, 1.0000, ...
                  1.0000, 1.0000, 1.0000, 1.0000];
[eig_moduli_sorted, ~] = sort(eig_moduli_raw, 'descend');

% Observability Gramian (Table 1) — MATLAB v3 verified
periods    = [1, 2, 3];
lam_min_O  = [1.2807,  1.6946,  2.0129];
lam_max_O  = [8.413e6, 1.461e9, 5.588e10];
kappa_O    = [6.569e6, 8.621e8, 2.776e10];

% Spectral lower bound — MATLAB v3 verified
sigma_min_M = 2.2941e-3;
lam_min_O1  = 1.2807;
bound_i2    = sigma_min_M^(2*(2-1)) * lam_min_O1;   % 6.74e-6
bound_i3    = sigma_min_M^(2*(3-1)) * lam_min_O1;   % 3.55e-11

% Scaling study (Table 2) — MATLAB v3 verified
n_vals     = [2, 3];
kappa_mean = [3.794e9,  1.243e14];
kappa_std  = [1.194e10, 3.874e14];

% ----------------------------------------------------------------
%  GLOBAL STYLE
% ----------------------------------------------------------------
set(0,'DefaultAxesFontSize',   13);
set(0,'DefaultAxesFontName',   'Times New Roman');
set(0,'DefaultTextFontName',   'Times New Roman');
set(0,'DefaultAxesLineWidth',  1.2);
set(0,'DefaultLineLineWidth',  2.0);

clr_blue   = [0.122, 0.471, 0.706];
clr_orange = [1.000, 0.498, 0.055];
clr_green  = [0.173, 0.627, 0.173];
clr_red    = [0.839, 0.153, 0.157];
clr_purple = [0.580, 0.404, 0.741];
clr_gray   = [0.5,   0.5,   0.5  ];

% ================================================================
%  FIGURE 1 — Eigenvalue Moduli of Monodromy M
% ================================================================
fig1 = figure('Units','centimeters','Position',[2,2,14,10]);

b = bar(1:9, eig_moduli_sorted, 0.65);
b.FaceColor = 'flat';
for k = 1:9
    if eig_moduli_sorted(k) > 1.001
        b.CData(k,:) = clr_red;
    elseif abs(eig_moduli_sorted(k) - 1.0) < 0.002
        b.CData(k,:) = clr_orange;
    else
        b.CData(k,:) = clr_blue;
    end
end

hold on;
yline(1.0,'--k','LineWidth',1.5);
text(8.5, 1.04, '$|\lambda|=1$','Interpreter','latex','FontSize',11);

% Annotate unstable eigenvalue
text(1, eig_moduli_sorted(1)+0.07, ...
    sprintf('$%.4f$', eig_moduli_sorted(1)), ...
    'Interpreter','latex','FontSize',10,...
    'HorizontalAlignment','center','Color',clr_red);

% Annotate stable eigenvalue
stab_idx = find(eig_moduli_sorted < 0.999, 1);
if ~isempty(stab_idx)
    text(stab_idx, eig_moduli_sorted(stab_idx)+0.07, ...
        sprintf('$%.4f$', eig_moduli_sorted(stab_idx)), ...
        'Interpreter','latex','FontSize',10,...
        'HorizontalAlignment','center','Color',clr_blue);
end

hold on;
hb1 = bar(NaN,NaN,'FaceColor',clr_red);
hb2 = bar(NaN,NaN,'FaceColor',clr_orange);
hb3 = bar(NaN,NaN,'FaceColor',clr_blue);
hl  = plot(NaN,NaN,'--k','LineWidth',1.5);
legend([hb1,hb2,hb3,hl], ...
    {'$|\lambda|>1$','$|\lambda|=1$','$|\lambda|<1$','Unit circle'}, ...
    'Interpreter','latex','Location','northeast','FontSize',10,'Box','on');

xlabel('Eigenvalue index $j$','Interpreter','latex');
ylabel('$|\lambda_j(\mathcal{M})|$','Interpreter','latex');
title({'Moduli of Monodromy Eigenvalues';
       '$\mathcal{M}=\Phi_{\mathcal{A}}(2\pi,0)$, $n=3$, Observability system'}, ...
    'Interpreter','latex');
xlim([0.5,9.5]); ylim([0,1.85]);
xticks(1:9); grid on; grid minor; box on;

set(fig1,'PaperUnits','centimeters','PaperSize',[14,10],...
    'PaperPosition',[0,0,14,10]);
print(fig1,'Fig1_Obs_EigModuli','-dpdf','-r500');
fprintf('Fig1_Obs_EigModuli.pdf saved.\n');

% ================================================================
%  FIGURE 2 — Observability Gramian: λ_min GROWTH across periods
% ================================================================
fig2 = figure('Units','centimeters','Position',[2,2,14,10]);

yyaxis left
plot(periods, lam_min_O, 'o-', 'Color', clr_green, ...
    'MarkerFaceColor', clr_green, 'MarkerSize', 9);
ylabel('$\lambda_{\min}(\widetilde{O}_i)$','Interpreter','latex');
ylim([1.0, 2.5]);

% Annotate lambda_min values
for i = 1:3
    text(periods(i)+0.05, lam_min_O(i)+0.04, ...
        sprintf('$%.4f$', lam_min_O(i)), ...
        'Interpreter','latex','FontSize',10,'Color',clr_green);
end

yyaxis right
semilogy(periods, lam_max_O, 's--', 'Color', clr_red, ...
    'MarkerFaceColor', clr_red, 'MarkerSize', 9);
ylabel('$\lambda_{\max}(\widetilde{O}_i)$','Interpreter','latex');

for i = 1:3
    text(periods(i)+0.05, lam_max_O(i)*1.6, ...
        sprintf('$%.3e$', lam_max_O(i)), ...
        'Interpreter','latex','FontSize',9,'Color',clr_red);
end

ax2 = gca;
ax2.YAxis(1).Color = clr_green;
ax2.YAxis(2).Color = clr_red;

xlabel('Period index $i$','Interpreter','latex');
title({'Observability Gramian Spectrum vs.\ Period';
       '$\lambda_{\min}$ grows (K-observable), $\lambda_{\max}$ grows geometrically'}, ...
    'Interpreter','latex');
legend({'$\lambda_{\min}(\widetilde{O}_i)$ (increasing — K-obs.\ confirmed)', ...
        '$\lambda_{\max}(\widetilde{O}_i)$ (geometric growth)'}, ...
    'Interpreter','latex','Location','west','FontSize',10,'Box','on');

xticks([1,2,3]); xlim([0.7,3.5]); grid on; box on;

set(fig2,'PaperUnits','centimeters','PaperSize',[14,10],...
    'PaperPosition',[0,0,14,10]);
print(fig2,'Fig2_Obs_GramianGrowth','-dpdf','-r500');
fprintf('Fig2_Obs_GramianGrowth.pdf saved.\n');

% ================================================================
%  FIGURE 3 — Spectral Lower Bound vs Actual λ_min
% ================================================================
fig3 = figure('Units','centimeters','Position',[2,2,13,10]);

% Actual lambda_min values
actual  = lam_min_O;          % [1.2807, 1.6946, 2.0129]

% Theoretical lower bounds
bounds  = [lam_min_O1, bound_i2, bound_i3];   % [1.2807, 6.74e-6, 3.55e-11]

semilogy(periods, actual, 'o-', 'Color', clr_green, ...
    'MarkerFaceColor', clr_green, 'MarkerSize', 9, 'LineWidth', 2.0);
hold on;
semilogy(periods, bounds, 's--', 'Color', clr_red, ...
    'MarkerFaceColor', clr_red, 'MarkerSize', 9, 'LineWidth', 2.0);

% Annotate actual values
for i = 1:3
    text(periods(i)+0.06, actual(i)*1.3, ...
        sprintf('$%.4f$', actual(i)), ...
        'Interpreter','latex','FontSize',10,'Color',clr_green);
end

% Annotate bound values
text(periods(2)+0.06, bounds(2)*3, ...
    sprintf('$%.2e$', bounds(2)), ...
    'Interpreter','latex','FontSize',10,'Color',clr_red);
text(periods(3)+0.06, bounds(3)*5, ...
    sprintf('$%.2e$', bounds(3)), ...
    'Interpreter','latex','FontSize',10,'Color',clr_red);

% Shade the gap to show bound is conservative
for i = 2:3
    plot([periods(i),periods(i)],[bounds(i),actual(i)],...
        ':','Color',clr_gray,'LineWidth',1.5);
end

% Add sigma_min annotation
dim = [0.18, 0.15, 0.35, 0.08];
annotation('textbox', dim, 'String', ...
    '$\sigma_{\min}(\mathcal{M}) \approx 2.29\times10^{-3}$', ...
    'Interpreter','latex','EdgeColor','none','FontSize',10,...
    'Color',clr_red);

xlabel('Period index $i$','Interpreter','latex');
ylabel('Value (log scale)','Interpreter','latex');
title({'Spectral Lower Bound vs.\ Actual $\lambda_{\min}(\widetilde{O}_i)$';
       'Bound~(17): $\lambda_{\min}(\widetilde{O}_i)\geq\sigma_{\min}(\mathcal{M})^{2(i-1)}\lambda_{\min}(\widetilde{O}_1)$'}, ...
    'Interpreter','latex');
legend({'Actual $\lambda_{\min}(\widetilde{O}_i)$', ...
        'Theoretical lower bound'}, ...
    'Interpreter','latex','Location','southwest','FontSize',10,'Box','on');

xticks([1,2,3]); xlim([0.7,3.7]); grid on; box on;

set(fig3,'PaperUnits','centimeters','PaperSize',[13,10],...
    'PaperPosition',[0,0,13,10]);
print(fig3,'Fig3_Obs_SpectralBound','-dpdf','-r500');
fprintf('Fig3_Obs_SpectralBound.pdf saved.\n');

% ================================================================
%  FIGURE 4 — Scaling Study: Condition Number κ(Ot_1) vs n
% ================================================================
fig4 = figure('Units','centimeters','Position',[2,2,12,9]);

errorbar(n_vals, kappa_mean, kappa_std, 'o-', ...
    'Color', clr_purple, 'MarkerFaceColor', clr_purple, ...
    'MarkerSize', 9, 'LineWidth', 2.0, 'CapSize', 10);

hold on;
% Log-linear fit
p   = polyfit(log10(n_vals), log10(kappa_mean), 1);
n_f = linspace(1.8, 3.2, 50);
plot(n_f, 10.^polyval(p,log10(n_f)), '--', ...
    'Color', clr_gray, 'LineWidth', 1.5);

% Annotate
for i = 1:2
    text(n_vals(i)+0.05, kappa_mean(i)*3, ...
        sprintf('$%.3e$', kappa_mean(i)), ...
        'Interpreter','latex','FontSize',10,'Color',clr_purple);
end

set(gca,'YScale','log');
xlabel('State dimension $n$','Interpreter','latex');
ylabel('$\kappa(\widetilde{O}_1)$ (mean $\pm$ std, 10 trials)','Interpreter','latex');
title({'Observability Gramian Conditioning vs.\ Dimension';
       'Large $\kappa$ reflects monodromy structure, not numerical error'}, ...
    'Interpreter','latex');
legend({'Mean $\pm$ std','Log-linear fit'}, ...
    'Interpreter','latex','Location','northwest','FontSize',10,'Box','on');

% n=4 excluded note
text(3.05, kappa_mean(2)*0.05, ...
    {'$n=4$ excluded:';'$256{\times}256$ ODE';'$>5\,\mathrm{GB}$'}, ...
    'Interpreter','latex','FontSize',9,'Color',clr_red,...
    'HorizontalAlignment','left');

xticks([2,3]);
xticklabels({'$n=2$\ ($n^2=4$)','$n=3$\ ($n^2=9$)'});
ax4 = gca;
ax4.XAxis.TickLabelInterpreter = 'latex';
xlim([1.7, 3.5]); grid on; box on;

set(fig4,'PaperUnits','centimeters','PaperSize',[12,9],...
    'PaperPosition',[0,0,12,9]);
print(fig4,'Fig4_Obs_ScalingStudy','-dpdf','-r500');
fprintf('Fig4_Obs_ScalingStudy.pdf saved.\n');

fprintf('\n=== All Paper 2 figures saved at 500 DPI ===\n');
fprintf('Files: Fig1_Obs_EigModuli.pdf\n');
fprintf('       Fig2_Obs_GramianGrowth.pdf\n');
fprintf('       Fig3_Obs_SpectralBound.pdf\n');
fprintf('       Fig4_Obs_ScalingStudy.pdf\n');
