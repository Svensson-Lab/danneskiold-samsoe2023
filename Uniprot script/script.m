%% PART 1: Import accession numbers
tic; clear;
opts = delimitedTextImportOptions("NumVariables", 7);
opts.VariableNames = ["c1", "c2", "c3", "c4", "c5", "c6", "c7"];
opts.SelectedVariableNames = ["c2", "c3"];
opts.DataLines = [468, Inf];
DATA = table2cell(readtable("input/cellphoneDB_interaction_curated_2023-01-19.csv", opts));
DATA = [DATA repmat({'CellPhoneDB'}, size(DATA, 1), 1)];

opts = spreadsheetImportOptions("NumVariables", 26);
opts.Sheet = "HumanGPCRLigands";
opts.VariableNames = ["c1", "c2", "c3", "c4", "c5", "c6", "Uniprot", "c8", "c9", "c10", "c11", "c12", "c13", "c14", "LigandPrUniprot", "c16", "c17", "c18", "c19", "c20", "c21", "c22", "c23", "c24", "c25", "c26"];
temp = table2cell(readtable("input/1-s2.0-S0092867419311262-mmc1.xlsx", opts, "UseExcel", false));
take = cellfun(@(x) ~isnan(str2double(x)), temp);
temp(take) = cellfun(@(x) {str2double(x)}, temp(take));
DATA = [DATA; temp(strcmp(temp(:, 12), 'Peptide'), [7 15]) repmat({'GPCRdb'}, sum(strcmp(temp(:, 12), 'Peptide')), 1)];

DATA(strcmp(DATA(:, 2), ''), :) = []; % removes 5 non-paired accession numbers 
DATA = DATA(~contains(DATA(:, 1), 'HLA'), :); % removes 4 HLAs
tabl = cell2table(cellfun(@num2str, DATA, 'Un', 0 ));
[~, ia] = unique(tabl(:, 1:2), 'rows');
DATA = DATA(ia, :); % removes 117 duplicates

save('output/DATApart1', 'DATA');
DATA = [{'accession', 'accession', 'reference'}; DATA];
writecell(DATA, 'output/DATApart1.xls');
toc; % fast

%% PART 2: Import gene names, single/multi-pass, and secreted or not
tic; clear;
load('output/DATApart1');
data = cell(size(DATA, 1), 6);
for i = 1:size(DATA, 1)*2
    disp(['Did ' num2str(i) ' out of ' num2str(size(DATA, 1)*2)]);
    TEMP = urlread(['https://rest.uniprot.org/uniprotkb/search?query=' DATA{i} '&format=tsv&fields=gene_primary,cc_subcellular_location']);
    gene = TEMP(48:48+strfind(TEMP(48:end), char(9))-2);
    pass = '';
    secr = 0;
    if contains(TEMP, 'Secreted')
        secr = 1;
    end
    if contains(TEMP, 'Single-pass')
        k = strfind(TEMP, "Single-pass");
        l = strfind(TEMP(k:end), ".");
        l2 = 999;
        if contains(TEMP(k:end), "{")
            l2 = strfind(TEMP(k:end), "{");
        end
        pass = strtrim(TEMP(k:k+min(l(1), l2(1))-2));
    elseif contains(TEMP, 'Multi-pass')
        k = strfind(TEMP, "Multi-pass");
        l = strfind(TEMP(k:end), ".");
        l2 = 999;
        if contains(TEMP(k:end), "{")
            l2 = strfind(TEMP(k:end), "{");
        end
        pass = strtrim(TEMP(k:k+min(l(1), l2(1))-2));
    end
    data{i} = gene;
    data{i+size(DATA, 1)*2} = pass;
    data{i+size(DATA, 1)*4} = secr;
end
DATA = [DATA(:, 1:2) data DATA(:, 3)];
save('output/DATApart2', 'DATA');
DATA = [{'accession', 'accession', 'gene_name', 'gene_name', 'receptor_type', 'receptor_type', 'secreted? yes/no', 'secreted? yes/no', 'reference'}; DATA];
writecell(DATA, 'output/DATApart2.xls');
toc; % takes 240 seconds

%% PART 3: Import gene and signal peptide lengths
tic; clear;
load('output/DATApart2');
DATA = DATA(sum(cell2mat(DATA(:, 7:8)), 2) == 1, :); % removes 202 where both or none are secreted
DATA(cell2mat(DATA(:, 7)) == 1, :) = DATA(cell2mat(DATA(:, 7)) == 1, [2 1 4 3 6 5 8 7 9]); % writes the receptor in the first column and the ligand in the second column
DATA = DATA(~cellfun(@isempty, DATA(:, 5)), :); % removes 12 with no receptor info
data = cell(size(DATA, 1), 1);
for i = 1:size(DATA, 1)
    disp(['Did ' num2str(i) ' out of ' num2str(size(DATA, 1))]);
    TEMP = urlread(['https://rest.uniprot.org/uniprotkb/search?query=' DATA{i+size(DATA, 1)} '&format=tsv&fields=length,ft_signal,ft_chain,ft_peptide']); % ft_chain,ft_peptide not used and can be removed
    len = str2double(TEMP(37:37+min(strfind(TEMP(37:end), char(9)))-2));
    tab = strfind(TEMP, char(9));
    sig = 0;
    if tab(5)-tab(4) > 1
        sig = str2double(TEMP(min(strfind(TEMP, '..'))+2:min(strfind(TEMP, '..'))+3));
    end
    data{i} = len-sig;
end
DATA = [DATA(:, 1:5) data DATA(:, 9)];
tabl = cell2table(cellfun(@num2str, DATA, 'Un', 0 ));
[~, ia] = unique(tabl(:, 1:2), 'rows');
DATA = DATA(ia, :); % removes 31 duplicates
save('output/DATApart3', 'DATA');
DATA = [{'receptor_accession', 'ligand_accession', 'receptor_name', 'ligand_name', 'receptor_type', 'ligand_gene_length', 'reference'}; DATA];
writecell(DATA, 'output/DATApart3.xls');
toc; % takes 90 seconds

%% Visualize data
tic; clear;
load('output/DATApart3');
DATA(cell2mat(cellfun(@(x)any(isnan(x)), DATA(:, 6), 'UniformOutput', false)) == 1, :) = []; % removes 1 with no gene length

figure(1); clf;
set(gcf, 'color', 'w', 'Position', [0 0 500 800]);
[~, ia] = unique(DATA(:, 4)); % removes 192 duplicates
violinplot(cell2mat(DATA(ia, 6)), extractBefore(DATA(ia, 5),' '), 'ShowBox', false, 'ShowWhiskers', false, 'ShowMedian', false, 'MarkerSize', 40, 'ViolinColor', [[53 183 121]/255; [72 40 120]/255]);
title('Gene length (AA)');
ylabel('Ligand length (amino acids)');
set(gca, 'FontSize', 20, 'LineWidth', 2, 'YTick', 0:200:800);
ylim([0 800]);

DATA = [{'receptor_accession', 'ligand_accession', 'receptor_name', 'ligand_name', 'receptor_type', 'ligand_gene_length', 'reference'}; DATA(ia, :)];
f = gcf; f.PaperSize = [f.PaperPosition(3) f.PaperPosition(4)];
writecell(DATA, 'output/gene_lengths.xls');
print('output/figure', '-dpdf');
toc;