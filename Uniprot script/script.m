%% Import data
tic; clear;
opts = delimitedTextImportOptions("NumVariables", 7);
opts.VariableNames = ["Var1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7"];
opts.SelectedVariableNames = ["Var2", "Var3"];
opts.DataLines = [468, Inf];
DATA = table2cell(readtable("cellphoneDB_interaction_curated_2023-01-19.csv", opts));
clear opts;

opts = spreadsheetImportOptions("NumVariables", 26);
opts.Sheet = "HumanGPCRLigands";
opts.VariableNames = ["EntryName", "Family", "Ligandtype", "Class", "Gene", "GeneAlternative", "Uniprot", "iupharId", "ligandId", "name", "species", "type", "mass", "LigandPrGene", "LigandPrUniprot", "LigandLength", "LigandPreLength", "SignalPeptide", "PepStartLoc", "PreAA", "AfterAA", "pKi_value", "pKd_value", "pEC50_value", "pIC50_value", "rankorderFlag"];
s2S1 = table2cell(readtable("1-s2.0-S0092867419311262-mmc1.xlsx", opts, "UseExcel", false));
numIdx = cellfun(@(x) ~isnan(str2double(x)), s2S1);
s2S1(numIdx) = cellfun(@(x) {str2double(x)}, s2S1(numIdx));
DATA = [DATA; s2S1(strcmp(s2S1(:, 12), 'Peptide'), [7 15])];
clear numIdx opts

CEL = cell(size(DATA, 1), 22);
for i = 1:size(DATA, 1)*2
    CEL{i} = DATA{i};
    TEMP = urlread(['https://rest.uniprot.org/uniprotkb/search?query=' DATA{i} '&format=tsv&fields=gene_primary']);
    CEL{i+size(DATA, 1)*2} = TEMP(22:end-1);
    TEMP = urlread(['https://rest.uniprot.org/uniprotkb/search?query=' DATA{i} '&format=tsv&fields=cc_subcellular_location']);
    txt1 = '';
    txt2 = 0;
    if contains(TEMP, 'Secreted')
        txt2 = 1;
    end
    if contains(TEMP, 'Single-pass')
        k = strfind(TEMP, "Single-pass");
        l = strfind(TEMP(k:end), ".");
        l2 = 999;
        if contains(TEMP(k:end), "{")
            l2 = strfind(TEMP(k:end), "{");
        end
        txt1 = strtrim(TEMP(k:k+min(l(1), l2(1))-2));
    end
    if contains(TEMP, 'Multi-pass')
        k = strfind(TEMP, "Multi-pass");
        l = strfind(TEMP(k:end), ".");
        l2 = 999;
        if contains(TEMP(k:end), "{")
            l2 = strfind(TEMP(k:end), "{");
        end
        txt1 = strtrim(TEMP(k:k+min(l(1), l2(1))-2));
    end
    CEL{i+size(DATA, 1)*4} = txt1;
    CEL{i+size(DATA, 1)*6} = txt2;
    TEMP = urlread(['https://rest.uniprot.org/uniprotkb/search?query=' DATA{i} '&format=tsv&fields=length']);
    LEN = str2double(TEMP(8:end-1));
    CEL{i+size(DATA, 1)*8} = LEN;
    TEMP = urlread(['https://rest.uniprot.org/uniprotkb/search?query=' DATA{i} '&format=tsv&fields=ft_chain']);
    CEL{i+size(DATA, 1)*10} = TEMP;
    txt3 = TEMP(13:strfind(TEMP, ";")-1);
    CEL{i+size(DATA, 1)*12} = txt3;
    CEL{i+size(DATA, 1)*14} = count(TEMP, 'CHAIN');
    if ~isempty(txt3)
        ol = split(txt3, '..');
        CEL{i+size(DATA, 1)*16} = str2double(cell2mat(ol(2))) - str2double(cell2mat(ol(1)));
    end
    TEMP = urlread(['https://rest.uniprot.org/uniprotkb/search?query=' DATA{i} '&format=tsv&fields=ft_signal']);
    SIG = TEMP(23:strfind(TEMP, ";")-1);
    CEL{i+size(DATA, 1)*18} = SIG;
    if ~isempty(SIG)
        SIG = str2double(SIG(end-1:end));
    else
        SIG = 0;
    end
    CEL{i+size(DATA, 1)*20} = LEN-SIG;
end
tak1 = cell2mat(CEL(:, 7));
tak2 = cell2mat(CEL(:, 8));
CEL(tak1 == 1, 23) = CEL(tak1 == 1, 17);
CEL(tak2 == 1, 23) = CEL(tak2 == 1, 18);
CEL(tak1 == 1, 24) = CEL(tak1 == 1, 21);
CEL(tak2 == 1, 24) = CEL(tak2 == 1, 22);
CEL(tak1 == 0, 25) = CEL(tak1 == 0, 5);
CEL(tak2 == 0, 25) = CEL(tak2 == 0, 6);

CELcopy = CEL;
CELcopy(2:end+1, :) = CELcopy(1:end, :);
CELcopy(1, :) = {'accession', 'accession', 'gene_name', 'gene_name', 'receptor_type', 'receptor_type', 'secreted', 'secreted', 'length', 'length', 'chainFullText', 'chainFullText', 'chain', 'chain', 'chainCOUNT', 'chainCOUNT', 'chainLength', 'chainLength', 'signal', 'signal', 'length-signal', 'length-signal', 'secretedChainLength', 'secretedLength-Signal', 'nonSecretedReceptorType'};
writecell(CELcopy, 'outputJAN27all.xls');
writecell(CELcopy, 'outputJAN27all');
save('resJAN27all', 'CEL');
toc;

%% Visualize data
tic; clear;
load('resJAN27all')
CEL = CEL(sum(cell2mat(CEL(:, 7:8)), 2) == 1, :);
CEL = CEL(~cellfun(@isempty,CEL(:, 23)), :);
CEL = CEL(~cellfun(@isempty,CEL(:, 25)), :);
CEL = CEL(~strcmp(CEL(:, 1), 'HLAC'), :);
tak3 = sum(cell2mat(CEL(:, 15:16)), 2) == 2;

figure(1); clf;
set(gcf, 'color', 'w', 'Position', [0 0 800 800]);
violinplot(log2(cell2mat(CEL(tak3, 23))), extractBefore(CEL(tak3, 25),' '), 'ShowBox', false, 'ShowWhiskers', false, 'ShowMedian', false, 'MarkerSize', 40, 'ViolinColor', [[53 183 121]/255; [72 40 120]/255]);
title('Chain length (AA)');
set(gca, 'FontSize', 20, 'LineWidth', 2, 'YTick', 4:2:11)
ylim([4 11.5]);

figure(2); clf;
set(gcf, 'color', 'w', 'Position', [800 0 800 800]);
violinplot(log2(cell2mat(CEL(tak3, 24))), extractBefore(CEL(tak3, 25),' '), 'ShowBox', false, 'ShowWhiskers', false, 'ShowMedian', false, 'MarkerSize', 40, 'ViolinColor', [[53 183 121]/255; [72 40 120]/255]);
title('Gene length (AA)');
set(gca, 'FontSize', 20, 'LineWidth', 2, 'YTick', 4:2:11)
ylim([4 11.5]);

CELcopy = CEL;
CELcopy(2:end+1, :) = CELcopy(1:end, :);
CELcopy(1, :) = {'accession', 'accession', 'gene_name', 'gene_name', 'receptor_type', 'receptor_type', 'secreted', 'secreted', 'length', 'length', 'chainFullText', 'chainFullText', 'chain', 'chain', 'chainCOUNT', 'chainCOUNT', 'chainLength', 'chainLength', 'signal', 'signal', 'length-signal', 'length-signal', 'secretedChainLength', 'secretedLength-Signal', 'nonSecretedReceptorType'};
writecell(CELcopy, 'outputJAN26plot.xls');
writecell(CELcopy, 'outputJAN26plot');
save('resJAN27plot', 'CEL');
toc;