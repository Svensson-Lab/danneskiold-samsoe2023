%% Work either with Tissue or Single cell type
TXT = "Tissue";
TXT = "Single_cell_type";
if strcmp(TXT, "Tissue")
    HTTP = "g,rnatd,t_RNA_adipose_tissue,t_RNA_adrenal_gland,t_RNA_amygdala,t_RNA_appendix,t_RNA_basal_ganglia,t_RNA_bone_marrow,t_RNA_breast,t_RNA_cerebellum,t_RNA_cerebral_cortex,t_RNA_cervix,t_RNA_choroid_plexus,t_RNA_colon,t_RNA_duodenum,t_RNA_endometrium_1,t_RNA_epididymis,t_RNA_esophagus,t_RNA_fallopian_tube,t_RNA_gallbladder,t_RNA_heart_muscle,t_RNA_hippocampal_formation,t_RNA_hypothalamus,t_RNA_kidney,t_RNA_liver,t_RNA_lung,t_RNA_lymph_node,t_RNA_medulla_oblongata,t_RNA_midbrain,t_RNA_ovary,t_RNA_pancreas,t_RNA_parathyroid_gland,t_RNA_pituitary_gland,t_RNA_placenta,t_RNA_pons,t_RNA_prostate,t_RNA_rectum,t_RNA_retina,t_RNA_salivary_gland,t_RNA_seminal_vesicle,t_RNA_skeletal_muscle,t_RNA_skin_1,t_RNA_small_intestine,t_RNA_smooth_muscle,t_RNA_spinal_cord,t_RNA_spleen,t_RNA_stomach_1,t_RNA_testis,t_RNA_thalamus,t_RNA_thymus,t_RNA_thyroid_gland,t_RNA_tongue,t_RNA_tonsil,t_RNA_urinary_bladder,t_RNA_vagina,t_RNA_white_matter";
    INI = 11;
    LOW = 0.52;
    YMAX1 = 840;
    YSEP1 = 200;
    YSEP2 = 10;
    YMAX2 = 55;
    FS = 26;
else
    HTTP = "g,rnatd,sc_RNA_Adipocytes,sc_RNA_Alveolar_cells_type_1,sc_RNA_Alveolar_cells_type_2,sc_RNA_Astrocytes,sc_RNA_B-cells,sc_RNA_Basal_keratinocytes,sc_RNA_Basal_prostatic_cells,sc_RNA_Basal_respiratory_cells,sc_RNA_Basal_squamous_epithelial_cells,sc_RNA_Bipolar_cells,sc_RNA_Breast_glandular_cells,sc_RNA_Breast_myoepithelial_cells,sc_RNA_Cardiomyocytes,sc_RNA_Cholangiocytes,sc_RNA_Club_cells,sc_RNA_Collecting_duct_cells,sc_RNA_Cone_photoreceptor_cells,sc_RNA_Cytotrophoblasts,sc_RNA_dendritic_cells,sc_RNA_Distal_enterocytes,sc_RNA_Distal_tubular_cells,sc_RNA_Ductal_cells,sc_RNA_Early_spermatids,sc_RNA_Endometrial_ciliated_cells,sc_RNA_Endometrial_stromal_cells,sc_RNA_Endothelial_cells,sc_RNA_Enteroendocrine_cells,sc_RNA_Erythroid_cells,sc_RNA_Excitatory_neurons,sc_RNA_Exocrine_glandular_cells,sc_RNA_Extravillous_trophoblasts,sc_RNA_Fibroblasts,sc_RNA_Gastric_mucus-secreting_cells,sc_RNA_Glandular_and_luminal_cells,sc_RNA_granulocytes,sc_RNA_Granulosa_cells,sc_RNA_Hepatocytes,sc_RNA_Hofbauer_cells,sc_RNA_Horizontal_cells,sc_RNA_Inhibitory_neurons,sc_RNA_Intestinal_goblet_cells,sc_RNA_Ionocytes,sc_RNA_Kupffer_cells,sc_RNA_Langerhans_cells,sc_RNA_Late_spermatids,sc_RNA_Leydig_cells,sc_RNA_Macrophages,sc_RNA_Melanocytes,sc_RNA_Microglial_cells,sc_RNA_monocytes,sc_RNA_Mucus_glandular_cells,sc_RNA_Muller_glia_cells,sc_RNA_NK-cells,sc_RNA_Oligodendrocyte_precursor_cells,sc_RNA_Oligodendrocytes,sc_RNA_Pancreatic_endocrine_cells,sc_RNA_Paneth_cells,sc_RNA_Peritubular_cells,sc_RNA_Plasma_cells,sc_RNA_Prostatic_glandular_cells,sc_RNA_Proximal_enterocytes,sc_RNA_Proximal_tubular_cells,sc_RNA_Respiratory_ciliated_cells,sc_RNA_Rod_photoreceptor_cells,sc_RNA_Salivary_duct_cells,sc_RNA_Schwann_cells,sc_RNA_Serous_glandular_cells,sc_RNA_Sertoli_cells,sc_RNA_Skeletal_myocytes,sc_RNA_Smooth_muscle_cells,sc_RNA_Spermatocytes,sc_RNA_Spermatogonia,sc_RNA_Squamous_epithelial_cells,sc_RNA_Suprabasal_keratinocytes,sc_RNA_Syncytiotrophoblasts,sc_RNA_T-cells,sc_RNA_Theca_cells,sc_RNA_Thymic_epithelial_cells,sc_RNA_Undifferentiated_cells";
    INI = 19;
    LOW = 0.61;
    YMAX1 = 720;
    YSEP1 = 200;
    YSEP2 = 50;
    YMAX2 = 340;
    FS = 23;
end

%% Import data and export it as a MAT-file
tic;
opts = delimitedTextImportOptions("NumVariables", 1);
opts.DataLines = [2, Inf];
opts.VariableNames = "Receptor";
R = readmatrix("input/single_pass_library.txt", opts);
clear opts;

R = R(~(strcmp(R, "CAD23") + strcmp(R, "CSMD1") + strcmp(R, "FAT1") + strcmp(R, "FAT2") + strcmp(R, "FAT3") + strcmp(R, "FAT4") + strcmp(R, "FRAS1") + strcmp(R, "FREM2") + strcmp(R, "LRP1") + strcmp(R, "LRP1B") + strcmp(R, "LRP2") + strcmp(R, "MUC12") + strcmp(R, "MUC16") + strcmp(R, "PKHD1")));

SHORT = [];
for i = 1:numel(R)
    if strcmp(R{i}, "CCSMST1"); r = 'C16orf91';
    elseif strcmp(R{i}, "CIROP"); r = 'LMLN2';
    elseif strcmp(R{i}, "FCGR1BP"); r = 'FCGR1B';
    elseif strcmp(R{i}, "HIDE1"); r = 'C19orf38';
    elseif strcmp(R{i}, "KCT2"); r = 'C5orf15';
    elseif strcmp(R{i}, "KIR2DL2"); r = 'KIR3DL2';
    else; r = R{i};
    end
    try TEMP = webread(['https://www.proteinatlas.org/api/search_download.php?search=' r '&format=json&columns=' HTTP{1} '&compress=no'], weboptions('Timeout', 120));
    catch; disp(i);
    end
    if numel(TEMP)
        TEMP(1).Gene = R{i};
        SHORT = [SHORT; TEMP(1, :)];
    else
        disp(R{i});
    end
end
save(['output/' TXT{1} '_SHORT'], 'SHORT');
toc; 

%% Convert SHORT to LONG format
tic;
load(['output/' TXT{1} '_SHORT']);
LONG = [];
for i = 1:size(SHORT, 1)
    t1 = struct2cell(SHORT(i));
    t2 = t1(3:end);
    f1 = fieldnames(SHORT(i));
    f2 = f1(3:end);
    N = sum(~strcmp(t2, '0.0') & ~cellfun(@isempty,t2));
    t = cell(N, 1);
    u = cell(N, 1);
    j = 1;
    for k = find(~strcmp(t2, '0.0') & ~cellfun(@isempty,t2))'
        u{j} = f2{k}(INI:end-6);
        j = j + 1;
    end
    for k = 1:N
        t{k} = t1{1};
    end
    TEMP = [u t t2(~strcmp(t2, '0.0') & ~cellfun(@isempty,t2))];
    LONG = [LONG; TEMP];
end
save(['output/' TXT{1} '_LONG'], 'LONG');
toc;

%% Convert SHORT to MAT format
tic;
load(['output/' TXT{1} '_SHORT']);
MAT = [];
for i = 1:size(SHORT, 1)
    t1 = struct2cell(SHORT(i));
    t2 = t1(3:end);
    f1 = fieldnames(SHORT(i));
    f2 = f1(3:end);
    N = sum(~strcmp(t2, '0.0') & ~cellfun(@isempty,t2));
    t = zeros(N, 1);
    u = zeros(N, 1);
    j = 1;
    for k = find(~strcmp(t2, '0.0') & ~cellfun(@isempty,t2))'
        u(j) = k;
        j = j + 1;
    end
    for k = 1:N
        t(k) = i;
    end
    TEMP = [u t str2double(t2(~strcmp(t2, '0.0') & ~cellfun(@isempty,t2)))];
    MAT = [MAT; TEMP];
end
save(['output/' TXT{1} '_MAT'], 'MAT');
toc;

%% Make supplementary figure
tic;
load(['output/' TXT{1} '_MAT']);
load(['output/' TXT{1} '_LONG']);

figure(1); clf;
tightsubplot(1, 1, [0 0], [LOW 0.01],  [0.053 0.048])
set(gcf, 'color', 'w', 'Position', [0 0 1670 600]);
MAT = MAT(MAT(:, 3) > 1, :);
G = groupsummary(MAT, MAT(:, 1), "mean");
N = histcounts(MAT(:, 1), 0.5:1:size(G, 1)+0.5);
[~, I] = sort(G(:, 3));
CAT = strrep(strrep(regexprep(LONG(:, 1), '([A-Z])', ' $1'), '1', ' '), '_', ' ');

bar(reordercats(categorical(CAT(I, 1)), CAT(I, 1)), N(I), 'FaceColor', [53 183 121]/255);
ylabel('#Receptors');
ylim([0 YMAX1]);
set(gca, 'FontSize', FS, 'YTick', 0:YSEP1:YMAX1, 'YColor', [72 40 120]/255);

yyaxis right;
bar(reordercats(categorical(CAT(I, 1)), CAT(I, 1)), G(I, 3), 'FaceColor', [72 40 120]/255);
ylim([0 YMAX2]);
ylabel('Mean nTPM');
xtickangle(90);

set(gca, 'FontName', 'sansserif', 'FontSize', FS, 'YTick', 0:YSEP2:YMAX2, 'YColor', [72 40 120]/255);
saveas(gcf, ['output/' TXT{1} '_FIGURE'], 'svg');
toc;