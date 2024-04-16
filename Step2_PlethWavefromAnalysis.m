clc
clear
close all hidden

%% Define file path
cd(file_path) %this is the save path from step 1

%%
files = dir('*.mat');

%%
CurData = [];
for f = 1:numel(files)
    cd(file_path)
    fname = files(f).name;
    load(fname)

    CurData = [CurData; BreathMat];
end

%%
new_ms = fs_target/1000;

%%
[COEFF, score, LATENT, TSQUARED, EXPLAINED, MU] = pca(CurData);
per_expl = 0;
for z_per = 1:numel(EXPLAINED)
    per_expl = per_expl + EXPLAINED(z_per);
    if per_expl >= 90
        %             disp(z_per)
        break
    else
    end
end
Z = linkage(score(:,1:z_per),'ward','euclidean','savememory','on');

%% Marsiske cutoff
max2consider = 30;

PerClust = ([1:max2consider]/max2consider)';
PerCoeff = flipud(Z(:,3))/max(Z(:,3));
PerCoeff = PerCoeff(1:max2consider);
Minus1PerCoeff = 1-PerCoeff;
MarsiskeIndex = abs(PerClust + PerCoeff - 1)';

[Y_Marsiske, I_Marsiske] = max(MarsiskeIndex);

cutoff = Z(end-(I_Marsiske-1),3)+1;

figure(1000)
subplot(2,2,1)
[h,t_tre] = dendrogram(Z,'ColorThreshold',cutoff);
hold on
xlims = get(gca,'xlim');
plot([xlims(1) xlims(2)],[cutoff cutoff],'r--')
xt = get(gca,'xtick');
xt_l = get(gca,'xticklabel');
for h_m = 1:numel(xt)
    xt_ln(h_m) = str2double(xt_l(h_m,:));
end
ylabel('Inner-squared Euclidean distance')


for h_m = 1:numel(h)
    lin_col{h_m} = [num2str(h(h_m).Color(1)) num2str(h(h_m).Color(2)) num2str(h(h_m).Color(3))];
    lin_x(h_m,:) = h(h_m).XData;
    lin_y(h_m,:) = h(h_m).YData;
end

%% set cluster number
if numel(unique(lin_col)) > 1
    breath_num_calc = numel(unique(lin_col))-1;
else
    breath_num_calc = 1;
end
figure(2801)
[h2,t2] = dendrogram(Z,breath_num_calc);
close 2801
breath_num = breath_num_calc;

figure(3050)
subplot(4,1,1)
[H, T] = dendrogram(Z,breath_num);
close 3050

%%
U_T = unique(T);
for t = 1:numel(U_T)

    T_Wave = mean(CurData(T == U_T(t),:),1,'omitnan');
    x_T_wave = [1:numel(T_Wave)]/new_ms;


    figure(1000)
    set(gcf,'pos',[213   229   852   529])
    subplot(2,2,2)
    hold on
    plot(x_T_wave,T_Wave,'linewidth',2)
    xlabel('Time (ms)')
    text(x_T_wave(end),T_Wave(end),num2str(U_T(t)))
end


%%
figure(1000)
subplot(2,3,4)
scatter3(score(:,1),score(:,2),score(:,3),2,T,'filled')
xlabel(['PC1 - ' num2str(round(EXPLAINED(1))) '%'])
ylabel(['PC2 - ' num2str(round(EXPLAINED(2))) '%'])
zlabel(['PC3 - ' num2str(round(EXPLAINED(3))) '%'])
colormap(lines)

subplot(2,3,5)
scatter3(score(:,1),score(:,2),score(:,3),2,T,'filled')
xlabel(['PC1 - ' num2str(round(EXPLAINED(1))) '%'])
ylabel(['PC2 - ' num2str(round(EXPLAINED(2))) '%'])
zlabel(['PC3 - ' num2str(round(EXPLAINED(3))) '%'])
view([0 -90])
colormap(lines)

subplot(2,3,6)
scatter3(score(:,3),score(:,2),score(:,1),2,T,'filled')
xlabel(['PC3 - ' num2str(round(EXPLAINED(3))) '%'])
ylabel(['PC2 - ' num2str(round(EXPLAINED(2))) '%'])
zlabel(['PC1 - ' num2str(round(EXPLAINED(1))) '%'])
view([0 -90])
colormap(lines)


